#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
namespace fs = std::filesystem;

struct Agent {
    double x, y, angle;
};

struct Edge {
    int from, to;
};

struct ClusterData {
    vector<int> agents;
    vector<Edge> edges;
    double avg_degree;
    int cluster_size;
    double centroid_x, centroid_y;
    double box_size;
};

vector<int> find_cluster_sizes(const vector<int>& parents) {
    unordered_map<int, int> parent_count;
    for (int parent : parents) {
        parent_count[parent]++;
    }
    
    vector<int> sizes;
    for (const auto& pair : parent_count) {
        sizes.push_back(pair.second);
    }
    
    sort(sizes.rbegin(), sizes.rend()); // Sort descending
    return sizes;
}

vector<int> get_cluster_agents(const vector<int>& parents, int target_parent) {
    vector<int> agents;
    for (int i = 0; i < parents.size(); i++) {
        if (parents[i] == target_parent) {
            agents.push_back(i);
        }
    }
    return agents;
}

bool is_valid_cluster(const vector<int>& cluster_agents, const vector<Agent>& agents) {
    for (int agent_idx : cluster_agents) {
        if (agents[agent_idx].x < 1.0 || agents[agent_idx].y < 1.0) {
            return false;
        }
    }
    return true;
}

vector<Edge> get_cluster_edges(const vector<int>& cluster_agents, 
                               const vector<vector<int>>& neighbor_lists) {
    vector<Edge> edges;
    unordered_map<int, bool> in_cluster;
    
    for (int agent : cluster_agents) {
        in_cluster[agent] = true;
    }
    
    for (int agent : cluster_agents) {
        for (int neighbor : neighbor_lists[agent]) {
            if (in_cluster.count(neighbor) && agent < neighbor) { // Avoid duplicates
                edges.push_back({agent, neighbor});
            }
        }
    }
    
    return edges;
}

ClusterData process_files(const string& sim_file, const string& p_file, 
                         const string& neighbor_file) {
    ClusterData result;
    
    // Load simulation data
    vector<Agent> agents(32768);
    ifstream sim_stream(sim_file);
    if (!sim_stream.is_open()) {
        cerr << "Cannot open: " << sim_file << endl;
        return result;
    }
    
    for (int i = 0; i < 32768; i++) {
        sim_stream >> agents[i].x >> agents[i].y >> agents[i].angle;
    }
    sim_stream.close();
    
    // Load parent data
    vector<int> parents(32768);
    ifstream p_stream(p_file);
    if (!p_stream.is_open()) {
        cerr << "Cannot open: " << p_file << endl;
        return result;
    }
    
    for (int i = 0; i < 32768; i++) {
        p_stream >> parents[i];
    }
    p_stream.close();
    
    // Load neighbor data
    vector<vector<int>> neighbor_lists(32768);
    vector<int> neighbor_counts(32768);
    
    ifstream neighbor_stream(neighbor_file);
    if (!neighbor_stream.is_open()) {
        cerr << "Cannot open: " << neighbor_file << endl;
        return result;
    }
    
    string line;
    while (getline(neighbor_stream, line) && neighbor_lists.size() <= 32768) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        vector<int> nums;
        int num;
        while (iss >> num) {
            nums.push_back(num);
        }
        
        if (nums.size() >= 2) {
            int agent_idx = nums[0];
            int num_neighbors = nums[1];
            neighbor_counts[agent_idx] = num_neighbors;
            
            if (num_neighbors > 0 && nums.size() >= 2 + num_neighbors) {
                for (int i = 2; i < 2 + num_neighbors; i++) {
                    neighbor_lists[agent_idx].push_back(nums[i]);
                }
            }
        }
    }
    neighbor_stream.close();
    
    // Find clusters by parent
    unordered_map<int, vector<int>> clusters;
    for (int i = 0; i < 32768; i++) {
        clusters[parents[i]].push_back(i);
    }
    
    // Sort clusters by size
    vector<pair<int, int>> cluster_sizes; // (size, parent_id)
    for (const auto& pair : clusters) {
        cluster_sizes.push_back({pair.second.size(), pair.first});
    }
    sort(cluster_sizes.rbegin(), cluster_sizes.rend());
    
    // Find largest valid cluster
    for (const auto& cluster_info : cluster_sizes) {
        int parent_id = cluster_info.second;
        vector<int> cluster_agents = clusters[parent_id];
        
        if (is_valid_cluster(cluster_agents, agents)) {
            result.agents = cluster_agents;
            result.cluster_size = cluster_agents.size();
            
            // Calculate centroid
            double sum_x = 0, sum_y = 0;
            for (int agent_idx : cluster_agents) {
                sum_x += agents[agent_idx].x;
                sum_y += agents[agent_idx].y;
            }
            result.centroid_x = sum_x / cluster_agents.size();
            result.centroid_y = sum_y / cluster_agents.size();
            
            // Calculate box size
            double min_x = agents[cluster_agents[0]].x, max_x = min_x;
            double min_y = agents[cluster_agents[0]].y, max_y = min_y;
            
            for (int agent_idx : cluster_agents) {
                min_x = min(min_x, agents[agent_idx].x);
                max_x = max(max_x, agents[agent_idx].x);
                min_y = min(min_y, agents[agent_idx].y);
                max_y = max(max_y, agents[agent_idx].y);
            }
            
            double x_range = max_x - min_x;
            double y_range = max_y - min_y;
            double max_range = max(x_range, y_range);
            
            vector<double> box_sizes = {25, 50, 100, 200, 400};
            result.box_size = 400; // Default
            for (double bs : box_sizes) {
                if (max_range <= bs * 0.8) {
                    result.box_size = bs;
                    break;
                }
            }
            
            // Get edges
            result.edges = get_cluster_edges(cluster_agents, neighbor_lists);
            
            // Calculate average degree
            double total_degree = 0;
            for (int agent_idx : cluster_agents) {
                total_degree += neighbor_counts[agent_idx];
            }
            result.avg_degree = total_degree / cluster_agents.size();
            
            break; // Found largest valid cluster
        }
    }
    
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <set_number>" << endl;
        return 1;
    }
    
    int set_idx = stoi(argv[1]);
    string base_path = "../";
    string output_dir = "./cluster_data_set_" + to_string(set_idx);
    
    // Create output directory
    fs::create_directories(output_dir);
    
    cout << "Processing set_" << set_idx << "..." << endl;
    
    for (int sim_idx = 1; sim_idx <= 30; sim_idx++) {
        for (int t_val = 0; t_val <= 0; t_val += 50000) {
            // Construct file paths
            string sim_file = base_path + "set_" + to_string(set_idx) + 
                             "/sampled/sim_data/sim_" + 
                             (sim_idx < 10 ? "0" : "") + to_string(sim_idx) + 
                             "/sim_data_t_" + 
                             string(7 - to_string(t_val).length(), '0') + to_string(t_val) + ".txt";
            
            string p_file = base_path + "set_" + to_string(set_idx) + 
                           "/sampled/p_data/sim_" + 
                           (sim_idx < 10 ? "0" : "") + to_string(sim_idx) + 
                           "/p_vector_t_" + 
                           string(7 - to_string(t_val).length(), '0') + to_string(t_val) + ".txt";
            
            string neighbor_file = base_path + "set_" + to_string(set_idx) + 
                                  "/sampled/neighbor_data/sim_" + 
                                  (sim_idx < 10 ? "0" : "") + to_string(sim_idx) + 
                                  "/neighbor_data_t_" + 
                                  string(7 - to_string(t_val).length(), '0') + to_string(t_val) + ".txt";
            
            // Check if files exist
            if (!fs::exists(sim_file) || !fs::exists(p_file) || !fs::exists(neighbor_file)) {
                continue;
            }
            
            ClusterData cluster = process_files(sim_file, p_file, neighbor_file);
            
            if (cluster.agents.empty()) {
                continue;
            }
            
            // Write cluster data
            string output_file = output_dir + "/cluster_set_" + to_string(set_idx) + 
                               "_sim_" + (sim_idx < 10 ? "0" : "") + to_string(sim_idx) + 
                               "_t_" + string(7 - to_string(t_val).length(), '0') + to_string(t_val) + ".txt";
            
            ofstream out(output_file);
            out << cluster.cluster_size << " " << fixed << setprecision(3) 
                << cluster.avg_degree << " " << cluster.centroid_x << " " 
                << cluster.centroid_y << " " << cluster.box_size << endl;
            
            // Write agent coordinates
            ifstream sim_in(sim_file);
            vector<double> all_x(32768), all_y(32768);
            for (int i = 0; i < 32768; i++) {
                double angle;
                sim_in >> all_x[i] >> all_y[i] >> angle;
            }
            sim_in.close();
            
            for (int agent_idx : cluster.agents) {
                out << all_x[agent_idx] << " " << all_y[agent_idx] << endl;
            }
            
            // Write edges
            out << cluster.edges.size() << endl;
            for (const Edge& edge : cluster.edges) {
                out << all_x[edge.from] << " " << all_y[edge.from] << " "
                    << all_x[edge.to] << " " << all_y[edge.to] << endl;
            }
            
            out.close();
            
            cout << "Processed: sim_" << (sim_idx < 10 ? "0" : "") << sim_idx 
                 << "_t_" << string(7 - to_string(t_val).length(), '0') << t_val 
                 << " (N=" << cluster.cluster_size << ")" << endl;
        }
    }
    
    cout << "Processing complete for set_" << set_idx << "!" << endl;
    return 0;
}