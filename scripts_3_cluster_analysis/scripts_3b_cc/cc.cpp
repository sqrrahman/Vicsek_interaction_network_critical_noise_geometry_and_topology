#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <queue>
#include <iomanip>
#include <filesystem>
#include <omp.h>


namespace fs = std::filesystem;

int main() {
    std::vector<int> time_steps;
    for (int t = 0; t <= 2000; t += 100) time_steps.push_back(t);

    for (int sim = 1; sim <= 5; ++sim) {
        std::ostringstream sim_folder_ss, output_folder_ss;
        sim_folder_ss << "neighbor_data/sim_" << std::setw(2) << std::setfill('0') << sim << "/";
        output_folder_ss << "cc_data/sim_" << std::setw(2) << std::setfill('0') << sim << "/";
        std::string sim_folder = sim_folder_ss.str();
        std::string output_folder = output_folder_ss.str();
        fs::create_directories(output_folder);

        #pragma omp parallel for
        for (int idx = 0; idx < time_steps.size(); ++idx) {
            int t = time_steps[idx];
            std::ostringstream input_file_ss, output_file_ss;
            input_file_ss << sim_folder << "neighbor_data_t_" << std::setw(7) << std::setfill('0') << t << ".txt";
            output_file_ss << output_folder << "cc_data_t_" << std::setw(7) << std::setfill('0') << t << ".txt";
            std::string input_file = input_file_ss.str();
            std::string output_file = output_file_ss.str();

            std::ifstream fin(input_file);
            if (!fin.is_open()) {
                std::cerr << "Cannot open input file: " << input_file << "\n";
                continue; // skip to next file
            }

            // Read neighbor data
            std::vector<std::set<int>> neighbors;
            int N = 0;
            std::string line;
            while (std::getline(fin, line)) {
                std::istringstream ss(line);
                int node, k;
                ss >> node >> k;
                if (node >= neighbors.size()) neighbors.resize(node + 1);
                for (int i = 0; i < k; ++i) {
                    int nb;
                    ss >> nb;
                    neighbors[node].insert(nb);
                }
                N = std::max(N, node + 1);
            }
            fin.close();

            // Cluster detection (BFS)
            std::vector<int> cluster_id(N, -1);
            int current_cluster = 0;
            for (int i = 0; i < N; ++i) {
                if (cluster_id[i] != -1) continue;
                std::queue<int> q;
                q.push(i);
                cluster_id[i] = current_cluster;

                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    for (int v : neighbors[u]) {
                        if (cluster_id[v] == -1) {
                            cluster_id[v] = current_cluster;
                            q.push(v);
                        }
                    }
                }
                current_cluster++;
            }

            // Compute clustering coefficient for each cluster
            struct ClusterInfo {
                int size = 0;
                double cc_sum = 0.0;
                double parent_value = 0.0; // placeholder
            };
            std::vector<ClusterInfo> clusters(current_cluster);

            for (int i = 0; i < N; ++i) {
                int cid = cluster_id[i];
                const auto &nb = neighbors[i];
                int k = nb.size();
                double Ci = 0.0;

                if (k >= 2) {
                    int ni = 0;
                    for (auto it1 = nb.begin(); it1 != nb.end(); ++it1) {
                        for (auto it2 = std::next(it1); it2 != nb.end(); ++it2) {
                            if (neighbors[*it1].count(*it2)) ni++;
                        }
                    }
                    Ci = 2.0 * ni / (k * (k - 1));
                }

                clusters[cid].cc_sum += Ci;
                clusters[cid].size++;
            }

            // Write output
            std::ofstream fout(output_file);
            if (!fout.is_open()) {
                std::cerr << "Cannot open output file: " << output_file << "\n";
                continue;
            }
            fout << std::fixed << std::setprecision(6);
            for (auto &c : clusters) {
                double cluster_cc = c.size > 0 ? c.cc_sum / c.size : 0.0;
                fout << c.size << " " << cluster_cc << " " << c.parent_value << "\n";
            }
            fout.close();

            std::cout << "Processed: " << input_file << " -> " << output_file << "\n";
        }
    }

    return 0;
}
