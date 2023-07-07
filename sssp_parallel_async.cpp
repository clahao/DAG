//
// Created by moyu on 2023/6/8.
//
#include "Graph.cpp"
#include "timer.h"
#include "timer.cpp"
#include "parallel.h"



void sssp_async_plus(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=16){
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
                        cal_times ++;
                        dist[neigh[j]] = dist_old[neigh[j]];
                        if(bm.get(neigh[j]) == 0){
                            q.push(neigh[j]);
                            bm.set(neigh[j]);
                        }
                    }
                }
            }
        }
    }
//#pragma omp parallel for num_threads(1)
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_async_split(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=16){
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

        queue<int> need_process_matrix;
        while(s--){
//            unordered_set<int> need_update;
            int node = q.front();
            q.pop();
            bm.reset(node);
//            timer.Start();
//            vector<int> neigh = csr.neighbors(node);
//            vector<int> neighbors_w = csr.neighbors_weights(node);
//            process_time += timer.Finish();
            int outdegree = csr.row_ptr[node+1] - csr.row_ptr[node];
            int start_node = csr.row_ptr[node];
            if(outdegree > num){
                __m512i state_vector;
                alignas(64) int state[SIMD_SIZE];
                for(int i = 0; i < SIMD_SIZE; i ++){
                    state[i] = dist[node];
                }
                state_vector = _mm512_load_epi32(state);
                for(int i = 0; i < outdegree; i += SIMD_SIZE){
                    if(i+SIMD_SIZE < outdegree){
                        alignas(64) int matrix_col[SIMD_SIZE];
                        alignas(64) int result[SIMD_SIZE];

                        for(int j = i; j < min(i+SIMD_SIZE ,outdegree); j ++){
                            matrix_col[j-i] = csr.weights[start_node+j];
                        }
                        __m512i matrix_col_vector = _mm512_load_epi32(matrix_col);
//                    process_time += timer.Finish();
//                    timer.Start();
                        __m512i result_vector = _mm512_add_epi32(state_vector,matrix_col_vector);
//                    process_time += timer.Finish();
                        _mm512_store_epi32(result,result_vector);

                        for(int j = i; j < min(i+SIMD_SIZE ,outdegree); j ++){
                            int node_index = csr.col_idx[start_node+j];
                            dist_old[node_index] = min(dist_old[node_index],result[j-i]);
                            if(dist_old[node_index] < dist[node_index]){
//                            cal_times ++;
                                dist[node_index] = dist_old[node_index];
                                if(bm.get(node_index) == 0){
                                    q.push(node_index);
                                    bm.set(node_index);
                                }
                            }
                        }
                    }
                    else{
                        for(int j = i; j < min(i+SIMD_SIZE ,outdegree); j ++){
                            int v = csr.col_idx[start_node+j];
                            int vweight = csr.weights[start_node+j];
                            dist_old[v] = min(dist_old[v], vweight + dist[node]);
                            if(dist_old[v] < dist[v]){
//                        need_update.insert(v);
                                dist[v] = dist_old[v];
                                if(bm.get(v) == 0){
                                    q.push(v);
                                    bm.set(v);
                                }

                            }
                        }
                    }
                }
            }
            else{
                for(int j=0;j<outdegree;j++){
                    int v = csr.col_idx[start_node+j];
                    int vweight = csr.weights[start_node+j];
                    dist_old[v] = min(dist_old[v], vweight + dist[node]);
                    if(dist_old[v] < dist[v]){
//                        need_update.insert(v);
                        dist[v] = dist_old[v];
                        if(bm.get(v) == 0){
                            q.push(v);
                            bm.set(v);
                        }

                    }
                }
//                for(auto nu : need_update){
//                    if(dist[nu] > dist_old[nu]){
//                        dist[nu] = min(dist[nu], dist_old[nu]);
//                        dist_old[nu] = dist[nu];
//                        int part = get_partition(nu, num);
//                        if(bm.get(nu) == 0) {
//                            q.push(nu);
//                            bm.set(nu);
//                        }
//                    }
//                }
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
    sssp_async_split(cal_times,process_time,csr,start,global_dist,0,graph.num_nodes,num);
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