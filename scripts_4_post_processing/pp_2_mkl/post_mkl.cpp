#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <omp.h>

namespace fs = std::filesystem;

std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> result(num);
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    double step = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = std::pow(10, log_start + step * i);
    return result;
}

void read_matrix(const std::string& filename,
                 std::vector<double>& mass,
                 std::vector<double>& degree,
                 std::vector<double>& apl) {
    std::ifstream file(filename);
    std::string line;
    double m, k, l;
    while (std::getline(file, line)) {
        if (line[0] == '#') continue;
        std::istringstream iss(line);
        if (iss >> m >> k >> l) {
            mass.push_back(m);
            degree.push_back(k);
            apl.push_back(l);
        }
    }
}

int main() {
    int num_bins = 32;
    auto binedges = logspace(1, 32768, num_bins + 1);
    std::vector<double> binmids(num_bins);
    for (int i = 0; i < num_bins; ++i)
        binmids[i] = std::sqrt(binedges[i] * binedges[i + 1]);
    
    std::string out_dir = "mkl_binned";
    fs::create_directories(out_dir);
    
    #pragma omp parallel for schedule(dynamic)
    for (int set_num = 1; set_num <= 9; ++set_num) {
        std::string folder_root = "../../set_" + std::to_string(set_num) + "/sampled/mkl_data/";
        std::vector<std::vector<double>> all_k(num_bins);
        std::vector<std::vector<double>> all_l(num_bins);
        std::vector<int> cluster_counts(num_bins, 0);
        
        for (const auto& entry : fs::recursive_directory_iterator(folder_root)) {
            if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                std::vector<double> mass, k, l;
                read_matrix(entry.path().string(), mass, k, l);

                for (int i = 0; i < num_bins; ++i) {
                    std::vector<double> k_vals, l_vals;
                    for (size_t j = 0; j < mass.size(); ++j) {
                        if (mass[j] >= binedges[i] && mass[j] < binedges[i + 1]) {
                            k_vals.push_back(k[j]);
                            l_vals.push_back(l[j]);
                        }
                    }
                    cluster_counts[i] += k_vals.size();
                    if (!k_vals.empty()) all_k[i].push_back(std::accumulate(k_vals.begin(), k_vals.end(), 0.0) / k_vals.size());
                    if (!l_vals.empty()) all_l[i].push_back(std::accumulate(l_vals.begin(), l_vals.end(), 0.0) / l_vals.size());
                }
            }
        }
        
        
        
        std::vector<double> avg_k(num_bins, NAN), avg_l(num_bins, NAN);
        std::vector<double> std_k(num_bins, NAN), std_l(num_bins, NAN);
        
        for (int i = 0; i < num_bins; ++i) {
            if (!all_k[i].empty()) avg_k[i] = std::accumulate(all_k[i].begin(), all_k[i].end(), 0.0) / all_k[i].size();
            if (!all_l[i].empty()) avg_l[i] = std::accumulate(all_l[i].begin(), all_l[i].end(), 0.0) / all_l[i].size();
        }
        
        for (int i = 0; i < num_bins; ++i) {
            // a) degree
            size_t nk = all_k[i].size();
            if (nk >= 2) {
                double mean = avg_k[i];
                double sumsq = 0;
                for (double v : all_k[i])
                    sumsq += (v - mean) * (v - mean);
                std_k[i]  = std::sqrt(sumsq / (nk - 1));    // SD
            }
        
            // b) apl
            size_t nl = all_l[i].size();
            if (nl >= 2) {
                double mean = avg_l[i];
                double sumsq = 0;
                for (double v : all_l[i])
                    sumsq += (v - mean) * (v - mean);
                std_l[i]  = std::sqrt(sumsq / (nl - 1));   // SD
            }
        }

        std::ofstream out(out_dir + "/mkl_set_" + std::to_string(set_num) + ".txt");
        out << "mass\tcount\tdegree\tdegree_std\tapl\tapl_std\n";
        
        for (int i = 0; i < num_bins; ++i) {
            if (!all_k[i].empty()) {
                out << binmids[i]              << "\t"
                    << cluster_counts[i]       << "\t"
                    << avg_k[i]                << "\t"
                    << std_k[i]                << "\t"
                    << avg_l[i]                << "\t"
                    << std_l[i]                << "\n";
            }
        }
        out.close();
        std::cout << "Written mkl_set_" << set_num << ".txt\n";
    }
    
    


    return 0;
}