#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>

// Generate logarithmically spaced numbers
std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> result(num);
    double log_start = std::log10(start);
    double log_end   = std::log10(end);
    double step      = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = std::pow(10, log_start + step * i);
    return result;
}

// Structure to store k/<k> value and count
struct KOverKData {
    double k_over_k;
    double count; // raw count
};

// Process a single file and calculate per-file probability for each bin
void processFile(const std::string& filePath,
                 std::map<std::pair<int,int>, std::map<double,double>>& fileProbMap) {
    std::ifstream fin(filePath);
    if (!fin.is_open()) {
        std::cerr << "Cannot open file: " << filePath << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::vector<int> data;
        int value;
        while (ss >> value) data.push_back(value);

        if (data.size() >= 3) {
            int bin_min = data[0];
            int bin_max = data[1];
            std::vector<int> degree_counts(data.begin() + 2, data.end());

            int total_count = 0;
            double weighted_sum = 0.0;
            for (size_t k = 0; k < degree_counts.size(); ++k) {
                int count = degree_counts[k];
                int degree = k + 1;
                total_count += count;
                weighted_sum += degree * count;
            }

            if (total_count > 0) {
                double k_mean = weighted_sum / total_count;
                double bin_width = 0.1;
                std::map<double,double> prob_map; // k/<k> -> probability

                for (size_t k = 0; k < degree_counts.size(); ++k) {
                    int count = degree_counts[k];
                    if (count > 0) {
                        int degree = k + 1;
                        double k_over_k = degree / k_mean;
                        double rounded_k = std::round(k_over_k / bin_width) * bin_width;
                        prob_map[rounded_k] += count;
                    }
                }

                // Normalize counts to probabilities for this file
                for (auto& pair : prob_map) {
                    pair.second /= total_count;
                    fileProbMap[{bin_min, bin_max}][pair.first] = pair.second;
                }
            }
        }
    }

    fin.close();
}

// Compute mean and std across multiple files
void computeMeanStd(const std::vector<std::map<double,double>>& fileProbs,
                    std::vector<double>& unique_k_over_k,
                    std::vector<double>& mean_probs,
                    std::vector<double>& std_probs) {
    if (fileProbs.empty()) return;

    // Collect all unique k/<k> values
    std::map<double,int> k_indices;
    int idx = 0;
    for (const auto& fmap : fileProbs) {
        for (const auto& p : fmap) {
            if (k_indices.find(p.first) == k_indices.end()) {
                k_indices[p.first] = idx++;
            }
        }
    }

    unique_k_over_k.resize(k_indices.size());
    mean_probs.assign(k_indices.size(), 0.0);
    std_probs.assign(k_indices.size(), 0.0);

    // Fill unique_k_over_k
    for (const auto& pair : k_indices) {
        unique_k_over_k[pair.second] = pair.first;
    }

    // Compute mean
    for (const auto& fmap : fileProbs) {
        for (const auto& pair : fmap) {
            int i = k_indices[pair.first];
            mean_probs[i] += pair.second;
        }
    }
    for (auto& p : mean_probs) p /= fileProbs.size();

    // Compute std
    for (const auto& fmap : fileProbs) {
        for (const auto& pair : fmap) {
            int i = k_indices[pair.first];
            double diff = pair.second - mean_probs[i];
            std_probs[i] += diff * diff;
        }
    }
    for (auto& s : std_probs) s = std::sqrt(s / (fileProbs.size() - 1));
}

// Write output
void writeOutput(const std::string& filename,
                 const std::vector<double>& k_over_k_values,
                 const std::vector<double>& mean_probs,
                 const std::vector<double>& std_probs) {
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        std::cerr << "Cannot create output file: " << filename << std::endl;
        return;
    }

    fout << "k/<k>\tmean_prob\tstd_prob\n";
    for (size_t i = 0; i < k_over_k_values.size(); ++i) {
        fout << std::fixed << std::setprecision(6) << k_over_k_values[i] << "\t"
             << std::scientific << std::setprecision(8) << mean_probs[i] << "\t"
             << std::scientific << std::setprecision(8) << std_probs[i] << "\n";
    }

    fout.close();
}

int main() {
    std::filesystem::create_directories("pdd_cpp");

    int num_bins = 32;
    auto binedges = logspace(1.0, 32768.0, num_bins + 1);

    // Generate bin ranges, skip first 10 bins
    std::vector<std::pair<int,int>> validBins;
    for (int i = 10; i < num_bins; ++i) {
        int bin_min = static_cast<int>(std::ceil(binedges[i]));
        int bin_max = static_cast<int>(std::floor(binedges[i + 1]));
        if (bin_max < bin_min) bin_max = bin_min;
        validBins.push_back({bin_min, bin_max});
    }

    for (int set_num = 1; set_num <= 9; ++set_num) {
        std::cout << "Processing set_" << set_num << "..." << std::endl;

        std::map<std::pair<int,int>, std::vector<std::map<double,double>>> allFileProbs;
        std::ostringstream set_dir;
        set_dir << "../../set_" << set_num << "/sampled/dd_data";

        if (!std::filesystem::exists(set_dir.str())) {
            std::cerr << "Warning: Directory " << set_dir.str() << " does not exist. Skipping." << std::endl;
            continue;
        }

        for (const auto& sim_entry : std::filesystem::directory_iterator(set_dir.str())) {
            if (sim_entry.is_directory() && sim_entry.path().filename().string().substr(0,4) == "sim_") {
                for (const auto& file_entry : std::filesystem::directory_iterator(sim_entry.path())) {
                    if (file_entry.is_regular_file()) {
                        std::string fname = file_entry.path().filename().string();
                        if (fname.substr(0,3) == "dd_" && fname.substr(fname.length()-4) == ".txt") {
                            std::map<std::pair<int,int>, std::map<double,double>> fileProbMap;
                            processFile(file_entry.path().string(), fileProbMap);
                            for (auto& pair : fileProbMap) {
                                allFileProbs[pair.first].push_back(pair.second);
                            }
                        }
                    }
                }
            }
        }

        // Compute mean and std for each bin
        for (const auto& binRange : validBins) {
            int bin_min = binRange.first;
            int bin_max = binRange.second;
            auto it = allFileProbs.find({bin_min, bin_max});
            if (it != allFileProbs.end() && !it->second.empty()) {
                std::vector<double> unique_k_over_k, mean_probs, std_probs;
                computeMeanStd(it->second, unique_k_over_k, mean_probs, std_probs);

                if (!unique_k_over_k.empty()) {
                    std::ostringstream output_file;
                    output_file << "pdd_cpp/dd_set_" << set_num << "_bin" << bin_min << "_" << bin_max << ".txt";
                    writeOutput(output_file.str(), unique_k_over_k, mean_probs, std_probs);
                    std::cout << "Written: " << output_file.str() << " (" << unique_k_over_k.size() << " unique k/<k>)\n";
                }
            }
        }

        std::cout << "Set_" << set_num << " completed\n" << std::endl;
    }

    std::cout << "All processing completed!" << std::endl;
    return 0;
}
