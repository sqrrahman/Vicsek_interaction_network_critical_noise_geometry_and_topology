#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <iomanip>
#include <omp.h>

using namespace std;

// Function to read neighbor data and calculate degree distribution
map<int, int> calculateDegreeDistribution(const string& filename) {
    map<int, int> degreeCount;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return degreeCount;
    }
    
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string first_col, second_col;
        
        // Read first column (agent ID) and second column (degree)
        if (iss >> first_col >> second_col) {
            int degree = stoi(second_col);
            degreeCount[degree]++;
        }
    }
    
    file.close();
    return degreeCount;
}

// Function to write degree distribution to output stream
void writeDegreeDistribution(ofstream& outFile, const map<int, int>& degreeCount, int maxDegree) {
    for (int deg = 0; deg <= maxDegree; deg++) {
        auto it = degreeCount.find(deg);
        int count = (it != degreeCount.end()) ? it->second : 0;
        outFile << count;
        if (deg < maxDegree) outFile << " ";
    }
    outFile << endl;
}

int main() {
    const int NUM_SETS = 9;
    const int NUM_SIMS = 5;
    const int TIME_START = 0;
    const int TIME_END = 2000;
    const int TIME_STEP = 100;
    const int NUM_TIME_POINTS = (TIME_END - TIME_START) / TIME_STEP + 1; // 201
    
    cout << "Processing " << NUM_SETS << " sets with " << NUM_SIMS << " simulations each..." << endl;
    cout << "Time points: " << NUM_TIME_POINTS << " (from " << TIME_START << " to " << TIME_END << " step " << TIME_STEP << ")" << endl;
    
    // Process each set
    for (int set = 1; set <= NUM_SETS; set++) {
        cout << "Processing set " << set << "..." << endl;
        
        // Create output file for this set
        string outputFilename = "set_" + to_string(set) + "_pk_over_k.txt";
        ofstream outFile(outputFilename);
        
        if (!outFile.is_open()) {
            cerr << "Error: Cannot create output file " << outputFilename << endl;
            continue;
        }
        
        // Process each simulation in this set
        for (int sim = 1; sim <= NUM_SIMS; sim++) {
            cout << "  Processing sim " << setfill('0') << setw(2) << sim << "..." << endl;
            
            // Store results for all time points of this simulation
            vector<map<int, int>> simResults(NUM_TIME_POINTS);
            int maxDegreeFound = 0;
            
            // Parallelize across time points
            #pragma omp parallel for reduction(max:maxDegreeFound)
            for (int t_idx = 0; t_idx < NUM_TIME_POINTS; t_idx++) {
                int t = TIME_START + t_idx * TIME_STEP;
                
                // Construct filename
                string filename = "../../set_" + to_string(set) + 
                                "/sampled/neighbor_data/sim_" + 
                                (sim < 10 ? "0" : "") + to_string(sim) + 
                                "/neighbor_data_t_" + 
                                string(7 - to_string(t).length(), '0') + to_string(t) + ".txt";
                
                // Calculate degree distribution for this file
                map<int, int> degreeCount = calculateDegreeDistribution(filename);
                simResults[t_idx] = degreeCount;
                
                // Find maximum degree in this file
                if (!degreeCount.empty()) {
                    int localMaxDegree = degreeCount.rbegin()->first;
                    if (localMaxDegree > maxDegreeFound) {
                        maxDegreeFound = localMaxDegree;
                    }
                }
            }
            
            // Write results for this simulation to output file
            for (int t_idx = 0; t_idx < NUM_TIME_POINTS; t_idx++) {
                writeDegreeDistribution(outFile, simResults[t_idx], maxDegreeFound);
            }
        }
        
        outFile.close();
        cout << "Completed set " << set << " -> " << outputFilename << endl;
    }
    
    cout << "All processing completed!" << endl;
    return 0;
}