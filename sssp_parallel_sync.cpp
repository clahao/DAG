//
// Created by moyu on 2023/5/31.
//
#include "Graph.cpp"
#include "timer.h"
#include "timer.cpp"
#include "parallel.h"

#define SIMD_SIZE 16

void sssp_sync(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    Timer timer;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    Bitmap bm(n);
    for(auto i:start){
        bm.set(i-start_id);
    }


    Bitmap bm_part(partition);
    queue<int> q;
    for(auto i:start){
        bm_part.set(get_partition(i-start_id, num));
        q.push(get_partition(i-start_id, num));
    }
    while (!q.empty()) {
        int s = q.size();
//        set<int> need_update;
        vector<bool> need_update(num_node,false);
        queue<int> need_process_matrix;
        while(s--){

            int part = q.front();
            q.pop();
            bm_part.reset(part);

            vector<int> vertexs = vertex_n_partion(part, n, num);
            vector<bool> bm_cnt(partition,false);

//            timer.Start();
            vector<vector<int>> adjacency_matrix(num,vector<int>(partition*num,INF));
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                if(bm.get(u) == 1){
                    vector<int> neigh = csr.neighbors(u);
                    vector<int> neighbors_w = csr.neighbors_weights(u);

                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int vweight = neighbors_w[j];
                        adjacency_matrix[i][v] = vweight;
                        int part_v = get_partition(v,num);
                        if(bm_cnt[part_v] == false){
                            need_process_matrix.push(part_v);
                            bm_cnt[part_v] = true;
                        }
                    }
                    bm.reset(u);
                }
            }
//            process_time += timer.Finish();
            while(!need_process_matrix.empty()){
                timer.Start();
                alignas(64) int state[num];
                __m512i state_vector;
                for(int i = 0; i < num; i ++){
                    int u = vertexs[i];
                    state[i] = dist[u];
                }
                state_vector = _mm512_load_epi32(state);
                process_time += timer.Finish();
                int part_v = need_process_matrix.front();
                need_process_matrix.pop();
                int start_node = part_v*num;
                for(int i = start_node; i < start_node + num; i ++){
                    alignas(64) int matrix_col[num];
                    alignas(64) int result[num];
                    timer.Start();
                    for(int j = 0; j < num; j ++){
                        matrix_col[j] = adjacency_matrix[j][i];
                    }
                    __m512i matrix_col_vector = _mm512_load_epi32(matrix_col);
                    process_time += timer.Finish();
//                    timer.Start();
                    __m512i result_vector = _mm512_add_epi32(state_vector,matrix_col_vector);
//                    process_time += timer.Finish();
                    timer.Start();
                    _mm512_store_epi32(result,result_vector);
                    process_time += timer.Finish();
//                    timer.Start();
                    for(int j = 0; j < num; j ++){
                        dist_old[i] = min(dist_old[i],result[j]);
                    }
//                    process_time += timer.Finish();
                    if(dist_old[i] < dist[i]){
                        need_update[i] = true;
                    }
                }
            }
        }
//        timer.Start();
        for(int nu = 0; nu < num_node; nu ++){
            if(need_update[nu]){
                if(dist[nu] > dist_old[nu]){
                    dist[nu] = min(dist[nu], dist_old[nu]);
                    dist_old[nu] = dist[nu];
                    int part = get_partition(nu, num);
                    if(bm_part.get(part) == 0){
                        q.push(part);
                        bm_part.set(part);
                    }
                    bm.set(nu);
                }
            }

        }
//        process_time += timer.Finish();
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
//#pragma omp parallel for num_threads(1)
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_sync_plus(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=16){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    Timer timer;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    Bitmap bm(n);
    queue<int> q;
    for(auto i:start){
        bm.set(i-start_id);
        q.push(i-start_id);
    }

    while (!q.empty()) {
        int s = q.size();
        unordered_set<int> need_update;
//        vector<bool> need_update(num_node,false);
        queue<int> need_process_matrix;
        while(s--){

            int node = q.front();
            q.pop();
            bm.reset(node);

            vector<int> neigh = csr.neighbors(node);
            vector<int> neighbors_w = csr.neighbors_weights(node);

            __m512i state_vector;
            alignas(64) int state[SIMD_SIZE];
            for(int i = 0; i < SIMD_SIZE; i ++){
                state[i] = dist[node];
            }
            state_vector = _mm512_load_epi32(state);
            for(int i = 0; i < neigh.size(); i += SIMD_SIZE){
                alignas(64) int matrix_col[SIMD_SIZE];
                alignas(64) int result[SIMD_SIZE];
                for(int j = i; j < i + SIMD_SIZE && j < neigh.size(); j ++){
                    matrix_col[j-i] = neighbors_w[j];
                }
                __m512i matrix_col_vector = _mm512_load_epi32(matrix_col);
                __m512i result_vector = _mm512_add_epi32(state_vector,matrix_col_vector);
                _mm512_store_epi32(result,result_vector);
                for(int j = i; j < i + SIMD_SIZE && j < neigh.size(); j ++){
                    dist_old[neigh[j]] = min(dist_old[neigh[j]],result[j-i]);
                    if(dist_old[neigh[j]] < dist[neigh[j]]){
//                        need_update[neigh[j]] = true;
                        need_update.insert(neigh[j]);
                    }
                }
            }
        }
//        for(int nu = 0; nu < num_node; nu ++){
//            if(need_update[nu]){
//                if(dist[nu] > dist_old[nu]){
//                    dist[nu] = min(dist[nu], dist_old[nu]);
//                    dist_old[nu] = dist[nu];
//                    if(bm.get(nu) == 0){
//                        q.push(nu);
//                        bm.set(nu);
//                    }
//                }
//            }
//        }
        for(auto nu : need_update){
            if(dist[nu] > dist_old[nu]){
                dist[nu] = min(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                if(bm.get(nu) == 0){
                    q.push(nu);
                    bm.set(nu);
                }
            }
        }
    }
//#pragma omp parallel for num_threads(1)
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

int main(int argc, char ** argv){
    string filename(argv[1]);
    int num = stoi(argv[2]);
    Graph<OutEdgeWeighted> graph(filename, true);
    int source_node = 0;    //以重排序前的0顶点为源顶点
    if(argc == 4){
        source_node = stoi(argv[3]);
    }
    CSR csr;
    csr.n = graph.num_nodes;
    csr.row_ptr.assign(csr.n+1, 0);
    csr.col_idx.assign(graph.num_edges, 0);
    csr.weights.assign(graph.num_edges, 0);

    vector<int> global_dist = vector<int>(graph.num_nodes,INF);
    global_dist[source_node] = 0;

    for(int i = 0; i <= graph.num_nodes; i ++)
        csr.row_ptr[i] = int(graph.offset[i]);
    for(int i = 0; i < graph.num_edges; i ++)
        csr.col_idx[i] = int(graph.edgeList[i].end);
    for(int i = 0; i < graph.num_edges; i ++)
        csr.weights[i] = int(graph.edgeList[i].w8);

    long long int cal_times = 0;
    set<int> start;
    start.insert(source_node);
    float process_time = 0;
    Timer timer;
    timer.Start();
//    sssp_sync(cal_times,process_time,csr,start,global_dist,0,graph.num_nodes,num);
    sssp_sync_plus(cal_times,process_time,csr,start,global_dist,0,graph.num_nodes,num);
    float runtime = timer.Finish();
    cout << "Calculating finished in " << process_time/1000 << " (s).\n";
    cout << "Processing finished in " << runtime/1000 << " (s).\n";
    cout<<"calculation times: "<<cal_times<<endl;
    int update_num = 0;
    for(int i = 0; i < graph.num_nodes; i ++){
        if(global_dist[i] != INF)
            update_num ++;
    }
    cout<<"updated vertex num: "<<update_num<<endl;
    for(int i = 0; i < min(100,(int)graph.num_nodes); i ++){
        cout<<global_dist[i]<<" ";
    }
    cout<<endl;


    return 0;

}