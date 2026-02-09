#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
#include <queue>


// Union-Find
struct UnionFind {
    vector<int> p, r;
    UnionFind(int n): p(n), r(n,0) { iota(p.begin(), p.end(), 0); }
    int  find(int x){ return p[x]==x?x:p[x]=find(p[x]); }
    void unite(int a,int b){
        a=find(a); b=find(b);
        if(a==b) return;

        if(r[a]<r[b]) p[a]=b;
        else if(r[b]<r[a]) p[b]=a;
        else { p[b]=a; r[a]++; }
    }
};

// Analyze one snapshot
// Build clusters and compute size, average degree, and average path length
void analyze_snapshot(const string &infile,
                      double r, double L,
                      const string &outfile,
                      const string &p_outfile)
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
    
    // 3) Collect edges
    vector<pair<int,int>> edges;
    
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
                    if(j<=i) continue;  // avoid double count
                    auto [xj,yj] = pos[j];
                    double ddx = fabs(xi-xj); if(ddx>0.5*L) ddx = L-ddx;
                    double ddy = fabs(yi-yj); if(ddy>0.5*L) ddy = L-ddy;
                    if(ddx*ddx + ddy*ddy <= r2){
                        edges.emplace_back(i, j);
                    }
                }
            }
        }
    }

    // 4) Process edges serially (thread-safe)
    UnionFind uf(N);
    vector<int> deg(N,0);
    vector<vector<int>> adj(N);
    
    for(const auto& edge : edges) {
        int i = edge.first, j = edge.second;
        uf.unite(i, j);
        deg[i]++;
        deg[j]++;
        adj[i].push_back(j);
        adj[j].push_back(i);
    }

    // 5) aggregate cluster stats
    unordered_map<int,int>    csize;
    unordered_map<int,int>    cdegsum;
    unordered_map<int, vector<int>> cluster_nodes;

    csize.reserve(N);
    cdegsum.reserve(N);
    cluster_nodes.reserve(N);

    for(int i=0;i<N;++i){
        int root = uf.find(i);
        csize[root]    += 1;
        cdegsum[root]  += deg[i];
        cluster_nodes[root].push_back(i);
    }

    unordered_map<int, double> avg_path_length;
    avg_path_length.reserve(N);

    for (const auto& kv : cluster_nodes) {
        int root = kv.first;
        const vector<int>& nodes = kv.second;
        int sz = nodes.size();

        if (sz < 2) {
            avg_path_length[root] = 0;
            continue;
        }        

        // Average path length via BFS
        double total_path = 0;
        int num_pairs = 0;

        for (int u : nodes) {
            queue<int> q;
            unordered_map<int, int> dist;
            dist[u] = 0;
            q.push(u);

            while (!q.empty()) {
                int curr = q.front(); q.pop();
                for (int v : adj[curr]) {
                    if (dist.find(v) == dist.end()) {
                        dist[v] = dist[curr] + 1;
                        q.push(v);
                        total_path += dist[v];
                        num_pairs++;
                    }
                }
            }
        }
        avg_path_length[root] = (num_pairs > 0) ? total_path / num_pairs : 0;
    }


    // 6) write out
    ofstream out(outfile);
    if(!out){
        cerr<<"Error opening "<<outfile<<"\n"; 
        return;
    }
    
    for(auto &kv : csize){
        int root = kv.first;
        int sz   = kv.second;
        double avg_deg = double(cdegsum[root]) / sz;
        out << sz << " "
            << avg_deg << " "
            << avg_path_length[root] << " "
            << root << "\n";
    }

    
    ofstream p_out(p_outfile);
    if(!p_out){
        cerr << "Error opening " << p_outfile << "\n";
        return;
    }
    
    for (int parent : uf.p) {
        p_out << parent << "\n";
    }
}

// Main: loop over snapshot times t = 0:10,000:2,000,000
int main(){


    float rho = {RHO_VALUE};
    float speed = {SPEED_VALUE};
    vector<int> sim_indices;
    for (int i = {START_IDX}; i <= {END_IDX}; ++i) sim_indices.push_back(i);
    const int t0 = 0;
    const int dt = 100;
    const int n_samples = 20;
    const int t1 = n_samples * dt;
        
    
    const double r = 1.0;         // Interaction radius
    const double L = sqrt(pow(2,9) / rho);     // Box size
    
    system("mkdir -p p_data mkl_data");
    
    const string prefix_in  = "sim_data/";
    const string p_prefix   = "p_data/";
    const string mkl_prefix = "mkl_data/";
    
    for (int sim_idx : sim_indices) {
        ostringstream sim_folder;
        sim_folder << setw(2) << setfill('0') << sim_idx;
        system((string("mkdir -p ") + p_prefix   + "sim_" + sim_folder.str()).c_str());
        system((string("mkdir -p ") + mkl_prefix + "sim_" + sim_folder.str()).c_str());
        
        #pragma omp parallel for schedule(dynamic)
        for(int t = t0; t <= t1; t += dt){
            char infile[300];
            snprintf(infile, sizeof(infile),
                "%ssim_%02d/sim_data_t_%07d.txt",
                prefix_in.c_str(), sim_idx, t);
            
            char mkl_outfile[300];
            snprintf(mkl_outfile, sizeof(mkl_outfile),
                "%ssim_%02d/mkl_data_t_%07d.txt",
                mkl_prefix.c_str(), sim_idx, t);
                     
            char p_outfile[300];
            snprintf(p_outfile, sizeof(p_outfile),
                "%ssim_%02d/p_vector_t_%07d.txt",
                p_prefix.c_str(), sim_idx, t);
            
            cout << "Processing t=" << t << " -> " << infile << " …\n";
            analyze_snapshot(string(infile), r, L, string(mkl_outfile), string(p_outfile));
        }
    }
    return 0;
}

