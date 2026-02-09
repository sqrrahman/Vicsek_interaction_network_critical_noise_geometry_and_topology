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
    
    std::string out_dir = "area_binned";
    fs::create_directories(out_dir);
    
    #pragma omp parallel for schedule(dynamic)
    for (int set_num = 1; set_num <= 9; ++set_num) {
        std::string folder_root = "../../set_" + std::to_string(set_num) + "/sampled/area_data/";
        std::vector<std::vector<double>> all_area(num_bins);
        std::vector<int> cluster_counts(num_bins, 0);
        
        for (const auto& entry : fs::recursive_directory_iterator(folder_root)) {
            if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                std::vector<double> mass, area, dummy;
                read_matrix(entry.path().string(), mass, area, dummy);


                for (int i = 0; i < num_bins; ++i) {
                    std::vector<double> area_vals;
                    for (size_t j = 0; j < mass.size(); ++j) {
                        if (mass[j] >= binedges[i] && mass[j] < binedges[i + 1]) {
                            area_vals.push_back(area[j]);
                        }
                    }
                    cluster_counts[i] += area_vals.size();
                    if (!area_vals.empty()) all_area[i].push_back(std::accumulate(area_vals.begin(), area_vals.end(), 0.0) / area_vals.size());
                }
            }
        }
        
        
        
        std::vector<double> avg_area(num_bins, NAN), std_area(num_bins, NAN);
        
        for (int i = 0; i < num_bins; ++i) {
            if (!all_area[i].empty()) avg_area[i] = std::accumulate(all_area[i].begin(), all_area[i].end(), 0.0) / all_area[i].size();
        }
        
        for (int i = 0; i < num_bins; ++i) {
            // a) area
            size_t narea = all_area[i].size();
            if (narea >= 2) {
                double mean = avg_area[i];
                double sumsq = 0;
                for (double v : all_area[i])
                    sumsq += (v - mean) * (v - mean);
                std_area[i]  = std::sqrt(sumsq / (narea - 1));    // SD
            }

        }

        std::ofstream out(out_dir + "/area_set_" + std::to_string(set_num) + ".txt");
        out << "mass\tcount\tarea\tarea_std\n";
        
        for (int i = 0; i < num_bins; ++i) {
            if (!all_area[i].empty()) {
                out << binmids[i]              << "\t"
                    << cluster_counts[i]       << "\t"
                    << avg_area[i]                << "\t"
                    << std_area[i]                << "\n";
            }
        }
        out.close();
        std::cout << "Written area_set_" << set_num << ".txt\n";
    }
    
    


    return 0;
}