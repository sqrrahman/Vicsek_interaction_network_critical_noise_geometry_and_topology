#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <omp.h>

namespace fs = std::filesystem;

// Log-spaced bin edges
std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> result(num);
    double log_start = std::log10(start);
    double log_end   = std::log10(end);
    double step      = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = std::pow(10, log_start + step * i);
    return result;
}

// Histogram count
std::vector<int> histcounts(const std::vector<double>& data,
                            const std::vector<double>& edges) {
    std::vector<int> counts(edges.size() - 1, 0);
    for (double val : data) {
        for (size_t i = 0; i + 1 < edges.size(); ++i) {
            if (val >= edges[i] && val < edges[i + 1]) {
                counts[i]++;
                break;
            }
        }
    }
    return counts;
}

int main() {
    //--- parameters
    int num_bins = 32;
    auto binedges = logspace(1.0, 32768.0, num_bins + 1);
    
    // integer bin widths for density normalization
    std::vector<int> binwidth(num_bins);
    for (int i = 0; i < num_bins; ++i)
        binwidth[i] = static_cast<int>(std::floor(binedges[i + 1]))
                    - static_cast<int>(std::ceil (binedges[i])) + 1;
    
    std::string out_dir = "cmpd_binned";
    fs::create_directories(out_dir);
    
    //--- loop over each set
    #pragma omp parallel for schedule(dynamic)
    for (int set_num = 1; set_num <= 9; ++set_num) {
        std::string folder_root = "../../set_" + std::to_string(set_num) + "/sampled/mkl_data/";
        std::vector<std::vector<int>>   all_counts;
        std::vector<std::vector<double>> all_probs;
        int file_count = 0;

        // read each .txt file and build histograms
        for (const auto& entry : fs::recursive_directory_iterator(folder_root)) {
            if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                std::ifstream in(entry.path());
                std::string line;
                std::vector<double> masses;
                while (std::getline(in, line)) {
                    if (line.empty() || line[0] == '#') continue;
                    std::istringstream iss(line);
                    double m, k, l;
                    if (iss >> m >> k >> l)
                        masses.push_back(m);
                }
                auto counts = histcounts(masses, binedges);
                all_counts.push_back(counts);
                file_count++;
                
                // compute this file’s probability‐density per bin
                double sample_total = std::accumulate(counts.begin(),
                                                      counts.end(), 0.0);
                std::vector<double> probs(num_bins, 0.0);
                if (sample_total > 0) {
                    for (int i = 0; i < num_bins; ++i)
                        probs[i] = (counts[i] / sample_total) / binwidth[i];
                }
                all_probs.push_back(probs);
            }
        }

        if (file_count == 0) {
            std::cerr << "No files found for set_" << set_num << "\n";
            continue;
        }

        // compute bin centers
        std::vector<double> bin_centers(num_bins);
        for (int i = 0; i < num_bins; ++i) {
            double log_mid = (std::log10(binedges[i]) +
                              std::log10(binedges[i + 1])) / 2.0;
            bin_centers[i] = std::pow(10, log_mid);
        }

        //--- compute mean & sample‐std‐dev for each bin
        int M = static_cast<int>(all_probs.size());
        std::vector<double> mean_prob(num_bins, 0.0),
                            std_prob (num_bins, 0.0);

        // mean
        for (int i = 0; i < num_bins; ++i) {
            for (auto& p : all_probs)
                mean_prob[i] += p[i];
            mean_prob[i] /= M;
        }
        // sample standard deviation
        for (int i = 0; i < num_bins; ++i) {
            double sumsq = 0.0;
            for (auto& p : all_probs)
                sumsq += (p[i] - mean_prob[i]) * (p[i] - mean_prob[i]);
            std_prob[i] = std::sqrt(sumsq / (M - 1));
        }

        //--- new: sum clusters across all files
        std::vector<int> total_clusters(num_bins, 0);
        for (auto &counts : all_counts)
            for (int i = 0; i < num_bins; ++i)
                total_clusters[i] += counts[i];

        //--- write results
        std::ofstream out(out_dir + "/cmpd_set_" +
                          std::to_string(set_num) + ".txt");
        out << "mass_center\tnum_clusters\tmean_probability_density\tstd_probability_density\n";
        for (int i = 0; i < num_bins; ++i) {
            if (mean_prob[i] > 0) {
                out << bin_centers[i] << "\t"
                    << total_clusters[i] << "\t"
                    << mean_prob[i] << "\t"
                    << std_prob[i] << "\n";
            }
        }
        out.close();
        std::cout << "Written cmpd_set_" << set_num << ".txt\n";
    }

    return 0;
}
