#ifndef M_PI
#define M_PI 3.1415926
#endif

#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <random>
#include <omp.h>
#include <cstdint>
#include <filesystem>
namespace fs = std::filesystem;


using namespace std;
using namespace chrono;

vector<float> generateSequenceRandom(int count) {
    vector<float> seq(count);
    static std::mt19937 rng(std::random_device{}()); 
    std::uniform_real_distribution<float> dist(0.00f, 1.00f);
    for (int i = 0; i < count; ++i) {
        seq[i] = dist(rng);
    }
    return seq;
}


bool readSimulationData(const string& filepath, vector<float>& xx, vector<float>& yy, vector<float>& theta, int N) {
    std::ifstream file(filepath);
    if (!file) {
        std::cerr << "Failed to open " << filepath << "\n";
        return false;
    }
    
    xx.resize(N);
    yy.resize(N);
    theta.resize(N);
    
    for (int i = 0; i < N; ++i) {
        if (!(file >> xx[i] >> yy[i] >> theta[i])) {
            std::cerr << "Failed to read data for particle " << i << " from " << filepath << "\n";
            return false;
        }
    }
    
    file.close();
    return true;
}


// Define the parameters for the simulation
struct Params {
    float noise;
    float r;  // Interaction radius
    float v;  // Speed of each particle
    int N;    // Number of particles
    float L;  // System size
};

float computeOrder(const vector<float>& theta, const Params& pm) {
    float sum_vx = 0.0f, sum_vy = 0.0f;
    #pragma omp parallel for reduction(+:sum_vx,sum_vy)
    for(int i = 0; i < pm.N; ++i) {
        sum_vx += pm.v * cosf(theta[i]);
        sum_vy += pm.v * sinf(theta[i]);
    }
    return sqrtf(sum_vx*sum_vx + sum_vy*sum_vy)/(pm.N * pm.v);
}

void update(vector<float>& xx, vector<float>& yy, vector<float>& theta, const Params& pm) {
    int N        = pm.N;
    float L      = pm.L;
    float L_half = 0.5f * L;
    float r2     = pm.r * pm.r;

    int cells_per_side = static_cast<int>(floor(L / pm.r));
    if (cells_per_side <= 0) {
        std::cerr << "Error: cells_per_side <= 0. L=" << L << ", r=" << pm.r << std::endl;
        exit(1);
    }

    float cell_size = L / cells_per_side;
    vector<int> head(cells_per_side * cells_per_side, -1);
    vector<int> next(N, -1);

    for (int i = 0; i < N; ++i) {
        int cx = std::min(static_cast<int>(xx[i] / cell_size), cells_per_side - 1);
        int cy = std::min(static_cast<int>(yy[i] / cell_size), cells_per_side - 1);
        int cell = cy * cells_per_side + cx;
        next[i] = head[cell];
        head[cell] = i;
    }
    
    // ───────────────────────────────────────────────────────────────────

    // 2) prepare new‐state buffers
    vector<float> new_xx(N), new_yy(N), new_theta(N);
    // precompute trig to avoid repeated cosf/sinf in the inner loop
    vector<float> c(N), s(N);
    for (int j = 0; j < N; ++j) {
        c[j] = cosf(theta[j]);
        s[j] = sinf(theta[j]);
    }
    //auto noise_seq = generateSequenceConst(N);
    auto noise_seq = generateSequenceRandom(N);


    // 3) parallel update, each i reads only old xx,yy,theta
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        float xi = xx[i], yi = yy[i];
        float sum_cos = 0.0f, sum_sin = 0.0f;

        // find neighbors and accumulate
        // only check the 3×3 neighbor cells
        int cx_i = static_cast<int>(xi / cell_size);
        int cy_i = static_cast<int>(yi / cell_size);
        for (int dxC = -1; dxC <= 1; ++dxC) {
            int ncx = (cx_i + dxC + cells_per_side) % cells_per_side;
            for (int dyC = -1; dyC <= 1; ++dyC) {
                int ncy = (cy_i + dyC + cells_per_side) % cells_per_side;
                int cell = ncy*cells_per_side + ncx;
                for (int j = head[cell]; j != -1; j = next[j]) {
                    float dx = xx[j] - xi;
                    if      (dx >  L_half) dx -= L;
                    else if (dx < -L_half) dx += L;
                    float dy = yy[j] - yi;
                    if      (dy >  L_half) dy -= L;
                    else if (dy < -L_half) dy += L;
                    if (dx*dx + dy*dy < r2) {
                        sum_cos += c[j];
                        sum_sin += s[j];
                    }
                }
            }
        }


        // compute new angle + noise
        float theta_mean    = atan2f(sum_sin, sum_cos);
        float noise_term    = pm.noise * 2.0f * M_PI * (noise_seq[i] - 0.5f);
        float theta_new_val = theta_mean + noise_term;
        new_theta[i]        = theta_new_val;

        // step positions, wrap with if‐branches
        float cx = cosf(theta_new_val), sy = sinf(theta_new_val);
        float x_new = xi + pm.v * cx;
        if      (x_new >= L) x_new -= L;
        else if (x_new <  0) x_new += L;
        new_xx[i] = x_new;

        float y_new = yi + pm.v * sy;
        if      (y_new >= L) y_new -= L;
        else if (y_new <  0) y_new += L;
        new_yy[i] = y_new;
    }

    // 4) write back
    xx    = move(new_xx);
    yy    = move(new_yy);
    theta = move(new_theta);
}


// Main function for the simulation
int main() {
    srand(time(0));

    
    float rho = {RHO_VALUE};
    float speed = {SPEED_VALUE};
    vector<float> noises;
    for (float n = {NOISE_VALUE}f; n <= 1.0f; n += 2.0f) {noises.push_back(n);}
    vector<int> sim_indices = {1};
    int total_time = 999;
    int sample_interval = 1;
    
    fs::create_directories("movie/order_params");
    fs::create_directories("movie/sim_data");

    for (float noise : noises) {
        
        cout << "Running simulation with rho=" << rho << ", s=" << speed << ", noise=" << noise << endl;
        
        for (int sim_idx : sim_indices) {

            // Simulation parameters
            Params pm;
            pm.noise = noise;
            pm.r = 1.00f;
            pm.v = speed;
            pm.N = pow(2,15);
            pm.L = sqrt(pm.N / rho);

            // Initialize particles' positions and angles
            vector<float> xx(pm.N), yy(pm.N), theta(pm.N);

            /*
            auto seq_x = generateSequenceRandom(pm.N);
            auto seq_y = generateSequenceRandom(pm.N);
            auto seq_theta = generateSequenceRandom(pm.N);

            for (int i = 0; i < pm.N; ++i) {
                xx[i] = seq_x[i] * pm.L;
                yy[i] = seq_y[i] * pm.L;
                theta[i] = seq_theta[i] * 2 * M_PI - M_PI;
            }
            */
            
            char initial_filepath[300];
            snprintf(initial_filepath, sizeof(initial_filepath),
                     "sim_data/sim_data_noise_%.2f_sim_%02d.txt",
                     noise, sim_idx);
            if (!readSimulationData(initial_filepath, xx, yy, theta, pm.N)) {
                std::cerr << "Failed to read initial simulation data from " << initial_filepath << "\n";
                return 1;
            }
            cout << "Successfully loaded initial conditions from " << initial_filepath << endl;
            
            char ts_filepath[300];
            snprintf(ts_filepath, sizeof(ts_filepath),
                     "movie/order_params/order_param_noise_%.2f_sim_%02d_timeseries.txt",
                     noise, sim_idx);
            std::ofstream ts_file(ts_filepath);
            if (!ts_file) {
                std::cerr << "Failed to open " << ts_filepath << "\n";
                return 1;
            }
            
            float last_P = 0.0f;
            
            float P = computeOrder(theta, pm);
            ts_file << P << "\n";
            
            
            char sim_filepath[300];
            snprintf(sim_filepath, sizeof(sim_filepath),
                     "movie/sim_data/sim_data_noise_%.2f_sim_%02d_t_%07d.txt",
                     noise, sim_idx, 0);
            std::ofstream sim_file(sim_filepath);
            if (!sim_file) {
                std::cerr << "Failed to open " << sim_filepath << "\n";
                return 1;
            }
            for (int i = 0; i < pm.N; ++i) {
                sim_file << xx[i] << " " << yy[i] << " " << theta[i] << "\n";
            }
            sim_file.close();
            
            for (int t = 1; t <= total_time; ++t) {
                update(xx, yy, theta, pm);
                float P = computeOrder(theta, pm);
                ts_file << P << "\n";
                
                // Save simulation data every sample_interval timesteps
                if (t % sample_interval == 0) {
                    char sim_filepath_t[300];
                    snprintf(sim_filepath_t, sizeof(sim_filepath_t),
                             "movie/sim_data/sim_data_noise_%.2f_sim_%02d_t_%07d.txt",
                             noise, sim_idx, t);
                    std::ofstream sim_file_t(sim_filepath_t);
                    if (!sim_file_t) {
                        std::cerr << "Failed to open " << sim_filepath_t << "\n";
                        return 1;
                    }
                    for (int i = 0; i < pm.N; ++i) {
                        sim_file_t << xx[i] << " " << yy[i] << " " << theta[i] << "\n";
                    }
                    sim_file_t.close();
                    
                    if (t % 100000 == 0) {  // Progress indicator
                        cout << "Completed timestep " << t << " / " << total_time << endl;
                    }
                }
            }
            
            ts_file.close();
            
            cout << "Simulation completed successfully for sim_idx=" << sim_idx << endl;
        }
    }
        return 0;
}
