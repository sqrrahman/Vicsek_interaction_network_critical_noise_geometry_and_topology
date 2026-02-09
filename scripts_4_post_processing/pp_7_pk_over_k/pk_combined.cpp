#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cmath>

using namespace std;

int main() {
    const int NUM_SETS = 9;
    const int NUM_FILES_PER_SET = 105; // 5 sims * 21 time points
    
    // Create output directory if it doesn't exist
    filesystem::create_directories("pk_calc");
    
    cout << "Processing degree count files to calculate mean probabilities..." << endl;
    
    // Process each set
    for (int set = 1; set <= NUM_SETS; set++) {
        cout << "Processing set " << set << "..." << endl;
        
        // Read input file for this set
        string inputFilename = "set_" + to_string(set) + "_pk_over_k.txt";
        ifstream inFile(inputFilename);
        
        if (!inFile.is_open()) {
            cerr << "Error: Cannot open input file " << inputFilename << endl;
            continue;
        }
        
        // Read all rows to determine maximum degree across all files
        vector<vector<int>> allCounts;
        string line;
        int maxDegree = 0;
        
        while (getline(inFile, line)) {
            vector<int> counts;
            istringstream iss(line);
            string count_str;
            
            while (iss >> count_str) {
                counts.push_back(stoi(count_str));
            }
            
            if (!counts.empty()) {
                maxDegree = max(maxDegree, (int)counts.size() - 1);
                allCounts.push_back(counts);
            }
        }
        inFile.close();
        
        if (allCounts.size() != NUM_FILES_PER_SET) {
            cerr << "Warning: Expected " << NUM_FILES_PER_SET << " rows but found " << allCounts.size() << " in " << inputFilename << endl;
        }
        
        // Calculate probabilities for each file and accumulate sums
        vector<double> probabilitySum(maxDegree + 1, 0.0);
        int validFiles = 0;
        
        for (const auto& counts : allCounts) {
            // Calculate total agents in this file
            int totalAgents = 0;
            for (int count : counts) {
                totalAgents += count;
            }
            
            if (totalAgents > 0) {
                // Calculate probabilities and add to sum
                for (int deg = 0; deg < (int)counts.size(); deg++) {
                    double probability = (double)counts[deg] / totalAgents;
                    probabilitySum[deg] += probability;
                }
                validFiles++;
            }
        }
        
        // Calculate mean probabilities
        vector<double> meanProbabilities(maxDegree + 1);
        for (int deg = 0; deg <= maxDegree; deg++) {
            meanProbabilities[deg] = (validFiles > 0) ? probabilitySum[deg] / validFiles : 0.0;
        }
        
        // Calculate variance for standard deviation
        vector<double> varianceSum(maxDegree + 1, 0.0);
        
        // Second pass through data to calculate variance
        for (const auto& counts : allCounts) {
            int totalAgents = 0;
            for (int count : counts) {
                totalAgents += count;
            }
            
            if (totalAgents > 0) {
                for (int deg = 0; deg < (int)counts.size(); deg++) {
                    double probability = (double)counts[deg] / totalAgents;
                    double diff = probability - meanProbabilities[deg];
                    varianceSum[deg] += diff * diff;
                }
            }
        }
        
        // Calculate standard deviations
        vector<double> stdDevs(maxDegree + 1);
        for (int deg = 0; deg <= maxDegree; deg++) {
            stdDevs[deg] = (validFiles > 1) ? sqrt(varianceSum[deg] / (validFiles - 1)) : 0.0;
        }
        
        // Write output file
        string outputFilename = "pk_calc/set_" + to_string(set) + "_pk_calc.txt";
        ofstream outFile(outputFilename);
        
        if (!outFile.is_open()) {
            cerr << "Error: Cannot create output file " << outputFilename << endl;
            continue;
        }
        
        // Write mean probabilities and standard deviations (two columns per line)
        for (int deg = 0; deg <= maxDegree; deg++) {
            outFile << fixed << setprecision(8) << meanProbabilities[deg] 
                    << " " << stdDevs[deg] << endl;
        }
        
        outFile.close();
        
        cout << "  Processed " << validFiles << " files, max degree: " << maxDegree << endl;
        cout << "  Output: " << outputFilename << endl;
    }
    
    cout << "All processing completed!" << endl;
    cout << "Output files created in pk_calc/ directory" << endl;
    
    return 0;
}