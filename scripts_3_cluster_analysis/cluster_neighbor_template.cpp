// Neighbor Analysis

#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

// Analyze one snapshot
void analyze_snapshot(const string &infile,
                      double r, double L,
                      const string &outfile)
{
    const double r2 = r*r;

    // 1) load positions
    vector<pair<double,double>> pos;
    {
        ifstream in(infile);
        if(!in) { 
            cerr<<"Error opening "<<infile<<"\n"; 
            return; 
        }
        double x,y,th;
        while(in>>x>>y>>th)
            pos.emplace_back(x,y);
    }
    int N = pos.size();
    if(N==0) {
        cerr<<"No particles in "<<infile<<"\n";
        return;
    }

    // 2) build periodic cell-list
    int ncell = max(1,int(floor(L/r)));
    double cell_size = L/ncell;
    vector<int> head(ncell*ncell, -1), nxt(N, -1);
    for(int i=0;i<N;++i){
        int cx = (static_cast<int>(pos[i].first  / cell_size) + ncell) % ncell;
        int cy = (static_cast<int>(pos[i].second / cell_size) + ncell) % ncell;
        assert(cx >= 0 && cx < ncell);
        assert(cy >= 0 && cy < ncell);
        int cid = cy*ncell + cx;
        nxt[i] = head[cid];
        head[cid] = i;
    }
    
    // 3) Find neighbors for each particle
    vector<vector<int>> neighbors(N);
    
    for(int i=0;i<N;++i){
        auto [xi, yi] = pos[i];
        int cx = (static_cast<int>(xi/cell_size) + ncell) % ncell;
        int cy = (static_cast<int>(yi/cell_size) + ncell) % ncell;
        assert(cx >= 0 && cx < ncell);
        assert(cy >= 0 && cy < ncell);

        for(int dx=-1; dx<=1; ++dx){
            int ncx = (cx+dx + ncell)%ncell;
            for(int dy=-1; dy<=1; ++dy){
                int ncy = (cy+dy + ncell)%ncell;
                int cid = ncy*ncell + ncx;
                for(int j=head[cid]; j!=-1; j=nxt[j]){
                    if(j==i) continue;  // skip self
                    auto [xj,yj] = pos[j];
                    double ddx = fabs(xi-xj); if(ddx>0.5*L) ddx = L-ddx;
                    double ddy = fabs(yi-yj); if(ddy>0.5*L) ddy = L-ddy;
                    if(ddx*ddx + ddy*ddy <= r2){
                        neighbors[i].push_back(j);
                    }
                }
            }
        }
        
        // Sort neighbors in ascending order
        sort(neighbors[i].begin(), neighbors[i].end());
    }

    // 4) Write output
    ofstream out(outfile);
    if(!out){
        cerr<<"Error opening "<<outfile<<"\n"; 
        return;
    }
    
    for(int i=0; i<N; ++i){
        out << i << " " << neighbors[i].size();
        for(int neighbor : neighbors[i]){
            out << " " << neighbor;
        }
        out << "\n";
    }
}

// Main: loop t = 0:10,000:2,000,000 ────────────────────────────────────────────
int main(){
    float rho = {RHO_VALUE};
    float speed = {SPEED_VALUE};
    vector<int> sim_indices;
    for (int i = 1; i <= 30; ++i) sim_indices.push_back(i);
    const int    t0 = 0;
    const int    t1 = 2000000;
    const int    dt = 10000;
    
    const double r = 1.0;         // Interaction radius
    const double L = sqrt(pow(2,15) / rho);     // Box size
    
    system("mkdir -p neighbor_data");
    
    const string prefix_in  = "sim_data/";
    const string neighbor_prefix = "neighbor_data/";
    
    for (int sim_idx : sim_indices) {
        ostringstream sim_folder;
        sim_folder << setw(2) << setfill('0') << sim_idx;
        system((string("mkdir -p ") + neighbor_prefix + "sim_" + sim_folder.str()).c_str());
        
        #pragma omp parallel for schedule(dynamic)
        for(int t = t0; t <= t1; t += dt){
            char infile[300];
            snprintf(infile, sizeof(infile),
                "%ssim_%02d/sim_data_t_%07d.txt",
                prefix_in.c_str(), sim_idx, t);
            
            char neighbor_outfile[300];
            snprintf(neighbor_outfile, sizeof(neighbor_outfile),
                "%ssim_%02d/neighbor_data_t_%07d.txt",
                neighbor_prefix.c_str(), sim_idx, t);
            
            cout << "Processing t=" << t << " -> " << infile << " …\n";
            analyze_snapshot(string(infile), r, L, string(neighbor_outfile));
        }
    }
    return 0;
}