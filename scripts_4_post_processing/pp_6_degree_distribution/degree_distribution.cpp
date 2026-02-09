#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>

// Generate logarithmically spaced numbers (for mass windows)
std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> result(num);
    double log_start = std::log10(start);
    double log_end   = std::log10(end);
    double step      = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = std::pow(10, log_start + step * i);
    return result;
}

// Generate linearly spaced edges: start, start+step, ... , end
std::vector<double> linspace_step(double start, double end, double step) {
    std::vector<double> v;
    for (double x = start; x <= end + 1e-12; x += step) {  // +eps to include end
        v.push_back(x);
    }
    return v;
}

// Process a single file and calculate per-file **probability density** for each bin
// now we pass k/<k> linear-bin edges
void processFile(const std::string& filePath,
                 const std::vector<double>& koverk_edges,
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
                int count  = degree_counts[k];
                int degree = static_cast<int>(k) + 1;
                total_count  += count;
                weighted_sum += degree * count;
            }

            if (total_count > 0) {
                double k_mean = weighted_sum / total_count;

                // k/<k> -> (mass)  (still unnormalized)
                std::map<double,double> prob_map;
                // k/<k> -> bin width
                std::map<double,double> width_map;

                for (size_t k = 0; k < degree_counts.size(); ++k) {
                    int count = degree_counts[k];
                    if (count <= 0) continue;

                    int degree = static_cast<int>(k) + 1;
                    double k_over_k = degree / k_mean;

                    // keep only within the linear range
                    if (k_over_k < koverk_edges.front() || k_over_k > koverk_edges.back())
                        continue;

                    // find bin: upper_bound -> first edge > value
                    auto it = std::upper_bound(koverk_edges.begin(), koverk_edges.end(), k_over_k);
                    int bin_idx = static_cast<int>(it - koverk_edges.begin()) - 1;
                    if (bin_idx < 0) bin_idx = 0;
                    if (bin_idx >= static_cast<int>(koverk_edges.size()) - 1)
                        bin_idx = static_cast<int>(koverk_edges.size()) - 2;

                    double x1 = koverk_edges[bin_idx];
                    double x2 = koverk_edges[bin_idx + 1];
                    // for linear bins, center is just midpoint
                    double bin_center = 0.5 * (x1 + x2);
                    double width = x2 - x1;   // this will be 0.1 everywhere

                    prob_map[bin_center] += count;
                    width_map[bin_center] = width;   // overwrite is fine
                }

                // Normalize counts to probability density for this file
                for (auto& pair : prob_map) {
                    double bin_center = pair.first;
                    double mass       = pair.second;                  // count in that bin
                    double width      = width_map[bin_center];        // 0.1
                    // probability mass:
                    double prob_mass  = mass / static_cast<double>(total_count);
                    // convert to density:
                    double density    = prob_mass / width;

                    // store in the outer structure
                    fileProbMap[{bin_min, bin_max}][bin_center] = density;
                }
            }
        }
    }

    fin.close();
}

// Compute mean and std across multiple files (now over densities)
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

    // Compute std (sample)
    for (const auto& fmap : fileProbs) {
        for (const auto& pair : fmap) {
            int i = k_indices[pair.first];
            double diff = pair.second - mean_probs[i];
            std_probs[i] += diff * diff;
        }
    }
    if (fileProbs.size() > 1) {
        for (auto& s : std_probs) s = std::sqrt(s / (fileProbs.size() - 1));
    } else {
        for (auto& s : std_probs) s = 0.0;
    }
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

    fout << "k/<k>\tmean_density\tstd_density\n";
    for (size_t i = 0; i < k_over_k_values.size(); ++i) {
        fout << std::fixed << std::setprecision(6) << k_over_k_values[i] << "\t"
             << std::scientific << std::setprecision(8) << mean_probs[i] << "\t"
             << std::scientific << std::setprecision(8) << std_probs[i] << "\n";
    }

    fout.close();
}

int main() {
    // output folder (changed name)
    const std::string out_dir = "dd_binned";
    std::filesystem::create_directories(out_dir);

    // mass bins (still log-spaced, but you set 128 here)
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

    // linear bins for k/<k>: 0, 0.1, 0.2, ..., 10.0
    std::vector<double> koverk_edges = linspace_step(0.0, 10.0, 0.1);

    for (int set_num = 1; set_num <= 9; ++set_num) {
        std::cout << "Processing set_" << set_num << "..." << std::endl;

        // (mass_window) -> vector of (k/<k> -> density) maps
        std::map<std::pair<int,int>, std::vector<std::map<double,double>>> allFileProbs;
        std::ostringstream set_dir;
        set_dir << "../../set_" << set_num << "/sampled/dd_data";

        if (!std::filesystem::exists(set_dir.str())) {
            std::cerr << "Warning: Directory " << set_dir.str() << " does not exist. Skipping." << std::endl;
            continue;
        }

        // scan sims
        for (const auto& sim_entry : std::filesystem::directory_iterator(set_dir.str())) {
            if (sim_entry.is_directory() && sim_entry.path().filename().string().substr(0,4) == "sim_") {
                for (const auto& file_entry : std::filesystem::directory_iterator(sim_entry.path())) {
                    if (file_entry.is_regular_file()) {
                        std::string fname = file_entry.path().filename().string();
                        if (fname.substr(0,3) == "dd_" && fname.substr(fname.length()-4) == ".txt") {
                            std::map<std::pair<int,int>, std::map<double,double>> fileProbMap;
                            processFile(file_entry.path().string(), koverk_edges, fileProbMap);
                            for (auto& pair : fileProbMap) {
                                allFileProbs[pair.first].push_back(pair.second);
                            }
                        }
                    }
                }
            }
        }

        // Compute mean and std for each bin (mass window)
        for (const auto& binRange : validBins) {
            int bin_min = binRange.first;
            int bin_max = binRange.second;
            auto it = allFileProbs.find({bin_min, bin_max});
            if (it != allFileProbs.end() && !it->second.empty()) {
                std::vector<double> unique_k_over_k, mean_probs, std_probs;
                computeMeanStd(it->second, unique_k_over_k, mean_probs, std_probs);

                if (!unique_k_over_k.empty()) {
                    std::ostringstream output_file;
                    output_file << out_dir << "/dd_set_" << set_num
                                << "_bin" << bin_min << "_" << bin_max << ".txt";
                    writeOutput(output_file.str(), unique_k_over_k, mean_probs, std_probs);
                    std::cout << "Written: " << output_file.str()
                              << " (" << unique_k_over_k.size() << " k/<k> bins)\n";
                }
            }
        }

        std::cout << "Set_" << set_num << " completed\n" << std::endl;
    }

    std::cout << "All processing completed!" << std::endl;
    return 0;
}
