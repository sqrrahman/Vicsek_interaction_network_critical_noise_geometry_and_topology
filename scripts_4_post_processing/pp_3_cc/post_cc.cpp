#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <string>

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
    
    std::string out_dir = "cc_binned";
    fs::create_directories(out_dir);
    
    #pragma omp parallel for schedule(dynamic)
    for (int set_num = 1; set_num <= 9; ++set_num) {
        std::string folder_root = "../../set_" + std::to_string(set_num) + "/sampled/cc_data/";
        std::vector<std::vector<double>> all_cc(num_bins);
        std::vector<int> cluster_counts(num_bins, 0);
        
        for (const auto& entry : fs::recursive_directory_iterator(folder_root)) {
            if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                std::vector<double> mass, cc, dummy;
                read_matrix(entry.path().string(), mass, cc, dummy);


                for (int i = 0; i < num_bins; ++i) {
                    std::vector<double> cc_vals;
                    for (size_t j = 0; j < mass.size(); ++j) {
                        if (mass[j] >= binedges[i] && mass[j] < binedges[i + 1]) {
                            cc_vals.push_back(cc[j]);
                        }
                    }
                    cluster_counts[i] += cc_vals.size();
                    if (!cc_vals.empty()) all_cc[i].push_back(std::accumulate(cc_vals.begin(), cc_vals.end(), 0.0) / cc_vals.size());
                }
            }
        }
        
        
        
        std::vector<double> avg_cc(num_bins, NAN), std_cc(num_bins, NAN);
        
        for (int i = 0; i < num_bins; ++i) {
            if (!all_cc[i].empty()) avg_cc[i] = std::accumulate(all_cc[i].begin(), all_cc[i].end(), 0.0) / all_cc[i].size();
        }
        
        for (int i = 0; i < num_bins; ++i) {
            // a) cc
            size_t ncc = all_cc[i].size();
            if (ncc >= 2) {
                double mean = avg_cc[i];
                double sumsq = 0;
                for (double v : all_cc[i])
                    sumsq += (v - mean) * (v - mean);
                std_cc[i]  = std::sqrt(sumsq / (ncc - 1));    // SD
            }

        }

        std::ofstream out(out_dir + "/cc_set_" + std::to_string(set_num) + ".txt");
        out << "mass\tcount\tcc\tcc_std\n";
        
        for (int i = 0; i < num_bins; ++i) {
            if (!all_cc[i].empty()) {
                out << binmids[i]              << "\t"
                    << cluster_counts[i]       << "\t"
                    << avg_cc[i]                << "\t"
                    << std_cc[i]                << "\n";
            }
        }
        out.close();
        std::cout << "Written cc_set_" << set_num << ".txt\n";
    }
    
    


    return 0;
}