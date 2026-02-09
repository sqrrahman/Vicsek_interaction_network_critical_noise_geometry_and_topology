#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <omp.h>
#include <string>

namespace fs = std::filesystem;

// log-spaced edges
std::vector<double> logspace(double start, double end, int num) {
    std::vector<double> result(num);
    double log_start = std::log10(start);
    double log_end   = std::log10(end);
    double step      = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = std::pow(10.0, log_start + step * i);
    return result;
}

// Read data: mass, L, nPerim, nInside, nLeaf, degP, degI, degL, Area, parent
void read_matrix(const std::string& filename,
                 std::vector<double>& mass,
                 std::vector<double>& L,
                 std::vector<double>& nPerim,
                 std::vector<double>& nInside,
                 std::vector<double>& nLeaf,
                 std::vector<double>& degP,
                 std::vector<double>& degI,
                 std::vector<double>& degL,
                 std::vector<double>& Area,
                 std::vector<double>& parent)
{
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double m, l, nP, nI, nL, dP, dI, dL, A, p;
        if (iss >> m >> l >> nP >> nI >> nL >> dP >> dI >> dL >> A >> p) {
            mass.push_back(m);
            L.push_back(l);
            nPerim.push_back(nP);
            nInside.push_back(nI);
            nLeaf.push_back(nL);
            degP.push_back(dP);
            degI.push_back(dI);
            degL.push_back(dL);
            Area.push_back(A);
            parent.push_back(p);
        }
    }
}

int main() {
    const int num_bins = 32;
    const auto binedges = logspace(1.0, 32768.0, num_bins + 1);
    std::vector<double> binmids(num_bins);
    for (int i = 0; i < num_bins; ++i)
        binmids[i] = std::sqrt(binedges[i] * binedges[i + 1]);

    const std::string out_dir = "perimeter_binned";
    fs::create_directories(out_dir);

    #pragma omp parallel for schedule(dynamic)
    for (int set_num = 1; set_num <= 9; ++set_num) {
        const std::string folder_root =
            "../../set_" + std::to_string(set_num) + "/sampled/perimeter_data/";

        // per-bin collections of per-file means
        std::vector<std::vector<double>> all_L(num_bins), all_A(num_bins), all_nP(num_bins), all_nI(num_bins), all_nL(num_bins);
        std::vector<std::vector<double>> all_degP(num_bins), all_degI(num_bins), all_degL(num_bins);
        std::vector<std::vector<double>> all_perim_over_area(num_bins);  // Added perimeter/area ratio
        std::vector<int> cluster_counts(num_bins, 0);

        for (const auto& entry : fs::recursive_directory_iterator(folder_root)) {
            if (!(entry.is_regular_file() && entry.path().extension() == ".txt")) continue;

            std::vector<double> mass, L, nP, nI, nL, dP, dI, dL, Area, parent;
            read_matrix(entry.path().string(), mass, L, nP, nI, nL, dP, dI, dL, Area, parent);

            for (int i = 0; i < num_bins; ++i) {
                std::vector<double> vL, vA, vnP, vnI, vnL, vdP, vdI, vdL, vPerimOverArea;
                for (size_t j = 0; j < mass.size(); ++j) {
                    if (mass[j] >= binedges[i] && mass[j] < binedges[i + 1]) {
                        vL.push_back(L[j]);
                        vA.push_back(Area[j]);
                        vnP.push_back(nP[j]);
                        vnI.push_back(nI[j]);
                        vnL.push_back(nL[j]);
                        vdP.push_back(dP[j]);
                        vdI.push_back(dI[j]);
                        vdL.push_back(dL[j]);
                        // Calculate perimeter/area ratio (avoid division by zero)
                        if (Area[j] > 0.0) {
                            vPerimOverArea.push_back(L[j] / Area[j]);
                        }
                    }
                }
                cluster_counts[i] += static_cast<int>(vL.size());
                auto mean_of = [](const std::vector<double>& v){
                    return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
                };
                if (!vL.empty())         { all_L[i].push_back(mean_of(vL));   }
                if (!vA.empty())         { all_A[i].push_back(mean_of(vA));    }
                if (!vnP.empty())        { all_nP[i].push_back(mean_of(vnP)); }
                if (!vnI.empty())        { all_nI[i].push_back(mean_of(vnI)); }
                if (!vnL.empty())        { all_nL[i].push_back(mean_of(vnL)); }
                if (!vdP.empty())        { all_degP[i].push_back(mean_of(vdP)); }
                if (!vdI.empty())        { all_degI[i].push_back(mean_of(vdI)); }
                if (!vdL.empty())        { all_degL[i].push_back(mean_of(vdL)); }
                if (!vPerimOverArea.empty()) { all_perim_over_area[i].push_back(mean_of(vPerimOverArea)); }
            }
        }

        // per-bin across-file means and population std (÷ n)
        auto set_mean = [](const std::vector<std::vector<double>>& all, std::vector<double>& avg){
            for (size_t i = 0; i < all.size(); ++i)
                if (!all[i].empty())
                    avg[i] = std::accumulate(all[i].begin(), all[i].end(), 0.0) / static_cast<double>(all[i].size());
        };
        auto set_std_pop = [](const std::vector<std::vector<double>>& all,
                              const std::vector<double>& avg,
                              std::vector<double>& sd){
            for (size_t i = 0; i < all.size(); ++i) {
                size_t n = all[i].size();
                if (n >= 1) {
                    double s2 = 0.0;
                    for (double v : all[i]) s2 += (v - avg[i]) * (v - avg[i]);
                    sd[i] = std::sqrt(s2 / static_cast<double>(n));
                }
            }
        };

        std::vector<double>
            avg_L(num_bins, NAN),   std_L(num_bins, NAN),
            avg_A(num_bins, NAN),   std_A(num_bins, NAN),
            avg_nP(num_bins, NAN),  std_nP(num_bins, NAN),
            avg_nI(num_bins, NAN),  std_nI(num_bins, NAN),
            avg_nL(num_bins, NAN),  std_nL(num_bins, NAN),
            avg_dP(num_bins, NAN),  std_dP(num_bins, NAN),
            avg_dI(num_bins, NAN),  std_dI(num_bins, NAN),
            avg_dL(num_bins, NAN),  std_dL(num_bins, NAN),
            avg_perim_over_area(num_bins, NAN), std_perim_over_area(num_bins, NAN);

        set_mean(all_L,   avg_L);   set_std_pop(all_L,   avg_L,   std_L);
        set_mean(all_A,   avg_A);   set_std_pop(all_A,   avg_A,   std_A);
        set_mean(all_nP,  avg_nP);  set_std_pop(all_nP,  avg_nP,  std_nP);
        set_mean(all_nI,  avg_nI);  set_std_pop(all_nI,  avg_nI,  std_nI);
        set_mean(all_nL,  avg_nL);  set_std_pop(all_nL,  avg_nL,  std_nL);
        set_mean(all_degP,avg_dP);  set_std_pop(all_degP,avg_dP,  std_dP);
        set_mean(all_degI,avg_dI);  set_std_pop(all_degI,avg_dI,  std_dI);
        set_mean(all_degL,avg_dL);  set_std_pop(all_degL,avg_dL,  std_dL);
        set_mean(all_perim_over_area, avg_perim_over_area); 
        set_std_pop(all_perim_over_area, avg_perim_over_area, std_perim_over_area);

        std::ofstream out(out_dir + "/perimeter_set_" + std::to_string(set_num) + ".txt");
        out << "mass\tcount\tL\tL_std\tnumPerimeterPoints\tnumPerimeterPoints_std\tnumInsidePoints\tnumInsidePoints_std\tnumLeafPoints\tnumLeafPoints_std\t"
               "avg_degree_perimeter\tavg_degree_perimeter_std\tavg_degree_inside\tavg_degree_inside_std\tavg_degree_leaf\tavg_degree_leaf_std\tarea\tarea_std\t"
               "perimeter_over_area\tperimeter_over_area_std\n";

        for (int i = 0; i < num_bins; ++i) {
            if (!std::isnan(avg_L[i])) {
                out << binmids[i]        << "\t"
                    << cluster_counts[i] << "\t"
                    << avg_L[i]  << "\t" << std_L[i]  << "\t"
                    << avg_nP[i] << "\t" << std_nP[i] << "\t"
                    << avg_nI[i] << "\t" << std_nI[i] << "\t"
                    << avg_nL[i] << "\t" << std_nL[i] << "\t"
                    << avg_dP[i] << "\t" << std_dP[i] << "\t"
                    << avg_dI[i] << "\t" << std_dI[i] << "\t"
                    << avg_dL[i] << "\t" << std_dL[i] << "\t"
                    << avg_A[i]  << "\t" << std_A[i] << "\t"
                    << avg_perim_over_area[i] << "\t" << std_perim_over_area[i] << "\n";
            }
        }
        out.close();
        std::cout << "Written perimeter_set_" << set_num << ".txt\n";
    }

    return 0;
}