#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <omp.h>

// Geometry point used for pixel stamps
struct Point {
    int x, y;
    Point(int x, int y) : x(x), y(y) {}
    bool operator==(const Point& other) const { return x == other.x && y == other.y; }
};

struct PointHash {
    size_t operator()(const Point& p) const {
        return std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1);
    }
};

struct Agent { double x, y, theta; };

struct ClusterResult {
    int size;
    double area;
    int parent;
};

// Precomputed circle pattern
namespace {
constexpr double kPixelSize = 0.05;
constexpr double kRadius    = 1.0;

const std::vector<Point>& get_circle_pattern() {
    static std::vector<Point> pat = []{
        std::vector<Point> p;
        const double r_pix    = kRadius / kPixelSize;
        const int    r_int    = static_cast<int>(std::ceil(r_pix));
        const double r_pix_sq = r_pix * r_pix;
        p.reserve((2*r_int + 1) * (2*r_int + 1));
        for (int dy = -r_int; dy <= r_int; ++dy) {
            for (int dx = -r_int; dx <= r_int; ++dx) {
                if (dx*dx + dy*dy <= r_pix_sq) p.emplace_back(dx, dy);
            }
        }
        return p;
    }();
    return pat;
}
} // namespace

// Area calculation by rasterizing disks around agent coords
double calculateArea(const std::vector<std::pair<double, double>>& coords) {
    if (coords.empty()) return 0.0;

    const double margin = kRadius + kPixelSize;

    auto [min_x, max_x] = std::minmax_element(
        coords.begin(), coords.end(),
        [](const auto& a, const auto& b){ return a.first < b.first; });

    auto [min_y, max_y] = std::minmax_element(
        coords.begin(), coords.end(),
        [](const auto& a, const auto& b){ return a.second < b.second; });

    const double xmin = min_x->first  - margin;
    const double ymin = min_y->second - margin;

    const auto& circle_pattern = get_circle_pattern();

    std::unordered_set<Point, PointHash> black_pixels;
    black_pixels.reserve(coords.size() * 12); // heuristic to reduce rehashes

    for (const auto& coord : coords) {
        int cx = static_cast<int>(std::round((coord.first  - xmin) / kPixelSize));
        int cy = static_cast<int>(std::round((coord.second - ymin) / kPixelSize));
        for (const auto& offset : circle_pattern) {
            int px = cx + offset.x;
            int py = cy + offset.y;
            if (px >= 0 && py >= 0) {
                black_pixels.insert(Point(px, py));
            }
        }
    }

    return black_pixels.size() * kPixelSize * kPixelSize;
}

// IO helpers
std::vector<int> readPVector(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<int> p_vector;
    int value;
    while (file >> value) p_vector.push_back(value);
    return p_vector;
}

std::vector<Agent> readSimData(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Agent> agents;
    double x, y, theta;
    while (file >> x >> y >> theta) agents.push_back({x, y, theta});
    return agents;
}

std::unordered_map<int, std::vector<int>> groupByParent(const std::vector<int>& p_vector) {
    std::unordered_map<int, std::vector<int>> clusters;
    clusters.reserve(p_vector.size());
    for (int i = 0; i < static_cast<int>(p_vector.size()); ++i) {
        clusters[p_vector[i]].push_back(i);
    }
    return clusters;
}

// Process one (sim, t) pair
void processFile(const std::string& p_folder,
                 const std::string& sim_folder,
                 const std::string& output_folder,
                 int sim, int t)
{
    // Build paths
    const std::string sim_str = (sim < 10 ? "0" : "") + std::to_string(sim);
    const std::string t_str   = std::string(7 - std::to_string(t).length(), '0') + std::to_string(t);

    const std::string p_filename =
        p_folder + "/sim_" + sim_str + "/p_vector_t_" + t_str + ".txt";
    const std::string sim_filename =
        sim_folder + "/sim_" + sim_str + "/sim_data_t_" + t_str + ".txt";
    const std::string output_dir =
        output_folder + "/sim_" + sim_str;
    const std::string output_filename =
        output_dir + "/area_data_t_" + t_str + ".txt";

    if (!std::filesystem::exists(p_filename) || !std::filesystem::exists(sim_filename)) {
        std::cerr << "Could not find matching files for sim " << sim << " t " << t << std::endl;
        return;
    }

    // Make sure the per-sim output dir exists
    std::filesystem::create_directories(output_dir);

    // Read & validate
    std::vector<int>    p_vector = readPVector(p_filename);
    std::vector<Agent>  agents   = readSimData(sim_filename);

    if (p_vector.empty() || agents.empty() || p_vector.size() != agents.size()) {
        std::cerr << "Error reading files or size mismatch for sim " << sim << " t " << t << std::endl;
        return;
    }

    // Group and compute areas for ALL clusters
    auto clusters = groupByParent(p_vector);

    std::vector<ClusterResult> results;
    results.reserve(clusters.size());
    
    for (const auto& [parent, cluster_agents] : clusters) {
        if (cluster_agents.size() <= 1) continue; // skip size 1 or empty
        std::vector<std::pair<double,double>> coords;
        coords.reserve(cluster_agents.size());
        for (int agent_id : cluster_agents) {
            coords.emplace_back(agents[agent_id].x, agents[agent_id].y);
        }
        double area = calculateArea(coords);
        results.push_back({static_cast<int>(cluster_agents.size()), area, parent});
    }

    // Write: <size> <area> <parent>
    std::ofstream out(output_filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot write to " << output_filename << std::endl;
        return;
    }
    for (const auto& r : results) {
        out << r.size << ' ' << r.area << ' ' << r.parent << '\n';
    }
    out.close();

    std::cout << "Processed: " << output_filename << " with " << results.size() << " clusters" << std::endl;
}

// Main function
int main() {
    const std::string p_folder = "p_data";
    const std::string sim_folder       = "sim_data";
    const std::string output_folder    = "area_data";

    std::filesystem::create_directories(output_folder);

    std::vector<int> sim_values;
    for (int i = 1; i <= 30; ++i) sim_values.push_back(i);

    std::vector<int> t_values;
    for (int t = 0; t <= 2000000; t += 10000) t_values.push_back(t);

    for (int sim : sim_values) {
        std::filesystem::create_directories(output_folder + "/sim_" + (sim < 10 ? "0" : "") + std::to_string(sim));
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(t_values.size()); ++i) {
            int t = t_values[i];
            processFile(p_folder, sim_folder, output_folder, sim, t);
        }
    }

    std::cout << "All processing complete." << std::endl;
    return 0;
}
