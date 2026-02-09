#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <map>
#include <iomanip>

// Generate logarithmically spaced values between start and end
std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> result(num);
    double log_start = std::log10(start);
    double log_end   = std::log10(end);
    double step      = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = std::pow(10, log_start + step * i);
    return result;
}

// main function: compute degree distributions binned by cluster size
int main() {
    int N = pow(2,9); // number of agents
    int num_bins = 32;
    int dt = 100;
    
    // Process all simulations (sim_01 to sim_5 or 30)
    for (int sim_num = 1; sim_num <= 5; ++sim_num) {
        for (int time_step = 0; time_step <= 2000; time_step += dt) {
            
            // Format file paths
            std::ostringstream sim_str, time_str;
            sim_str << "sim_" << std::setfill('0') << std::setw(2) << sim_num;
            time_str << "t_" << std::setfill('0') << std::setw(7) << time_step;
            
            std::string input_file = "neighbor_data/" + sim_str.str() + "/neighbor_data_" + time_str.str() + ".txt";
            std::string output_dir = "dd_data/" + sim_str.str();
            std::filesystem::create_directories(output_dir);
            std::string output_file = output_dir + "/dd_" + time_str.str() + ".txt";
            
            std::cout << "Processing: " << sim_str.str() << " " << time_str.str() << std::endl;
            
            // Read neighbor data
            std::vector<std::vector<int>> neighbors(N);
            std::vector<int> degrees(N, 0);
            
            std::ifstream fin(input_file);
            if (!fin.is_open()) {
                std::cerr << "Cannot open input file: " << input_file << "\n";
                continue; // Skip this file and continue with next
            }
            
            std::string line;
            while (std::getline(fin, line)) {
                std::istringstream ss(line);
                int idx, k;
                ss >> idx >> k;
                degrees[idx] = k;
                neighbors[idx].resize(k);
                for (int i = 0; i < k; ++i) ss >> neighbors[idx][i];
            }
            fin.close();
            
            // Find clusters using BFS
            std::vector<bool> visited(N, false);
            std::vector<std::vector<int>> clusters;
            
            for (int i = 0; i < N; ++i) {
                if (!visited[i]) {
                    std::vector<int> cluster;
                    std::queue<int> q;
                    q.push(i);
                    visited[i] = true;
                    
                    while (!q.empty()) {
                        int node = q.front(); q.pop();
                        cluster.push_back(node);
                        for (int nei : neighbors[node]) {
                            if (!visited[nei]) {
                                visited[nei] = true;
                                q.push(nei);
                            }
                        }
                    }
                    if (cluster.size() >= 3) clusters.push_back(cluster);
                }
            }
            
            // Bin edges
            auto binedges = logspace(1.0, 32768.0, num_bins + 1);
            
            // Prepare bins for degree values
            std::vector<std::vector<int>> bin_degrees(num_bins);
            
            // Assign clusters to bins and collect degrees
            for (auto &cluster : clusters) {
                int size = cluster.size();
                // Find bin index
                auto it = std::upper_bound(binedges.begin(), binedges.end(), size);
                int bin_idx = std::distance(binedges.begin(), it) - 1;
                if (bin_idx >= 0 && bin_idx < num_bins) {
                    for (int node : cluster) {
                        bin_degrees[bin_idx].push_back(degrees[node]);
                    }
                }
            }
            
            // Write output with bin ranges and degree counts
            std::ofstream fout(output_file);
            if (!fout.is_open()) {
                std::cerr << "Cannot open output file: " << output_file << "\n";
                continue; // Skip this file and continue with next
            }
            
            for (int i = 0; i < num_bins; ++i) {
                if (!bin_degrees[i].empty()) {
                    int bin_min = static_cast<int>(std::ceil(binedges[i]));
                    int bin_max = static_cast<int>(std::floor(binedges[i + 1]));
                    if (bin_max < bin_min) bin_max = bin_min;
                    fout << bin_min << " " << bin_max << " ";
                    
                    // Count occurrences of each degree
                    int max_degree = *std::max_element(bin_degrees[i].begin(), bin_degrees[i].end());
                    std::vector<int> counts(max_degree + 1, 0); // index 0 unused
                    for (int d : bin_degrees[i]) counts[d]++;
                    
                    for (int d = 1; d <= max_degree; ++d) {
                        fout << counts[d];
                        if (d != max_degree) fout << " ";
                    }
                    fout << "\n";
                }
            }
            fout.close();
            
            std::cout << "Completed: " << sim_str.str() << " " << time_str.str() << std::endl;
        }
    }
    
    std::cout << "All processing completed!" << std::endl;
    return 0;
}