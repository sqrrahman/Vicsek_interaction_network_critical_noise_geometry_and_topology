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
    float cell_size = L / cells_per_side;
    
    // 1) Move all particles first (parallel)
    vector<float> new_xx(N), new_yy(N), new_theta(N);
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        float xi = xx[i], yi = yy[i];
        
        // step positions first using current angle
        float cx = cosf(theta[i]), sy = sinf(theta[i]);
        float x_new = xi + pm.v * cx;
        if      (x_new >= L) x_new -= L;
        else if (x_new <  0) x_new += L;
        new_xx[i] = x_new;
        
        float y_new = yi + pm.v * sy;
        if      (y_new >= L) y_new -= L;
        else if (y_new <  0) y_new += L;
        new_yy[i] = y_new;
    }
    
    // 2) Rebuild cell list with new positions (sequential)
    vector<int> head(cells_per_side * cells_per_side, -1);
    vector<int> next(N, -1);
    for (int i = 0; i < N; ++i) {
        int cx = std::min(static_cast<int>(new_xx[i] / cell_size), cells_per_side - 1);
        int cy = std::min(static_cast<int>(new_yy[i] / cell_size), cells_per_side - 1);
        int cell = cy * cells_per_side + cx;
        next[i] = head[cell];
        head[cell] = i;
    }
    
    // 3) Update angles based on new positions (parallel)
    auto noise_seq = generateSequenceRandom(N);
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        float sum_cos = 0.0f, sum_sin = 0.0f;
        
        // find neighbors using new positions
        int cx_i = static_cast<int>(new_xx[i] / cell_size);
        int cy_i = static_cast<int>(new_yy[i] / cell_size);
        for (int dxC = -1; dxC <= 1; ++dxC) {
            int ncx = (cx_i + dxC + cells_per_side) % cells_per_side;
            for (int dyC = -1; dyC <= 1; ++dyC) {
                int ncy = (cy_i + dyC + cells_per_side) % cells_per_side;
                int cell = ncy*cells_per_side + ncx;
                
                for (int j = head[cell]; j != -1; j = next[j]) {
                    float dx = new_xx[j] - new_xx[i];
                    if      (dx >  L_half) dx -= L;
                    else if (dx < -L_half) dx += L;
                    
                    float dy = new_yy[j] - new_yy[i];
                    if      (dy >  L_half) dy -= L;
                    else if (dy < -L_half) dy += L;
                    
                    if (dx*dx + dy*dy < r2) {
                        sum_cos += cosf(theta[j]);
                        sum_sin += sinf(theta[j]);
                    }
                }
            }
        }
        
        // compute new angle + noise
        float theta_mean    = atan2f(sum_sin, sum_cos);
        float noise_term    = pm.noise * 2.0f * M_PI * (noise_seq[i] - 0.5f);
        new_theta[i]        = theta_mean + noise_term;
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
    vector<float> noises = {NOISE_VALUES};
    int total_time = pow(10,6);
    
    fs::create_directories("order_params");
    fs::create_directories("sim_data");

    for (float noise : noises) {
        std::ostringstream oss;
        oss << "order_params/order_param_noise_" << std::fixed << std::setprecision(2) << noise << ".txt";
        string noise_file = oss.str();
        std::ofstream noise_out(noise_file, std::ios::app);  // Append mode
        if (!noise_out) {std::cerr << "Failed to open " << noise_file << "\n"; return 1;}
        
        cout << "Running simulation with rho=" << rho << ", s=" << speed << ", noise=" << noise << endl;
        
        for (int sim_idx = 1; sim_idx <= 30; ++sim_idx) {
        
            char ts_filepath[300];
            snprintf(ts_filepath, sizeof(ts_filepath),
                     "order_params/order_param_noise_%.2f_sim_%02d_timeseries.txt",
                     noise, sim_idx);
                     
            char sim_filepath[300];
            snprintf(sim_filepath, sizeof(sim_filepath),
                     "sim_data/sim_data_noise_%.2f_sim_%02d.txt",
                     noise, sim_idx);
            
            // Skip if both files already exist

            if (fs::exists(ts_filepath) && fs::exists(sim_filepath)) {cout << "  Simulation " << sim_idx << " already exists, skipping..." << endl; continue;}

            // Simulation parameters
            Params pm;
            pm.noise = noise;
            pm.r = 1.00f;
            pm.v = speed;
            pm.N = pow(2,15);
            pm.L = sqrt(pm.N / rho);

            // Initialize particles' positions and angles
            vector<float> xx(pm.N), yy(pm.N), theta(pm.N);

            auto seq_x = generateSequenceRandom(pm.N);
            auto seq_y = generateSequenceRandom(pm.N);
            auto seq_theta = generateSequenceRandom(pm.N);

            for (int i = 0; i < pm.N; ++i) {
                xx[i] = seq_x[i] * pm.L;
                yy[i] = seq_y[i] * pm.L;
                theta[i] = seq_theta[i] * 2 * M_PI - M_PI;
            }
            
            
            std::ofstream ts_file(ts_filepath);
            if (!ts_file) {
                std::cerr << "Failed to open " << ts_filepath << "\n";
                return 1;
            }
            
            float last_P = 0.0f;
            
            float P = computeOrder(theta, pm);
            last_P = P;
            ts_file << P << "\n";

            for (int t = 2; t <= total_time; ++t) {
                update(xx, yy, theta, pm);
                float P = computeOrder(theta, pm);
                last_P = P;
                ts_file << P << "\n";
            }
            
            ts_file.close();
            
            noise_out << last_P << "\n";
            
            
            // Save final sim_data for this replicate
            
            std::ofstream sim_file(sim_filepath);
            if (!sim_file) {
                std::cerr << "Failed to open " << sim_filepath << "\n";
                return 1;
            }

            for (int i = 0; i < pm.N; ++i) {
                sim_file << xx[i] << " " << yy[i] << " " << theta[i] << "\n";
            }
            sim_file.close();
        }
        noise_out.close();
    }
    return 0;
}
