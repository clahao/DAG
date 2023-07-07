//
// Created by moyu on 2023/6/2.
//
#include "scc_reorder.cpp"
#include "timer.h"
#include "timer.cpp"
using namespace std;

struct cmp_sssp{
    bool operator()(pair<int,int> &a, pair<int,int> &b){
        if(a.second != b.second)
            return a.second > b.second;
        else
            return a.first < b.first;
    }
};

void sssp_DAG_sync_SIMD(long long int &cal_times, float &process_time, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    Timer timer;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
    queue<int> q;
    for(auto i:start){
        bm_part.set(get_partition(node_map[i-start_id], num));
        q.push(get_partition(node_map[i-start_id], num));
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

            timer.Start();
            vector<vector<int>> adjacency_matrix(num,vector<int>(partition*num,INF));
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                if(bm.get(u) == 1){
                    u = node_inv_map[u];
                    vector<int> neigh = csr.neighbors(u);
                    vector<int> neighbors_w = csr.neighbors_weights(u);

                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int vweight = neighbors_w[j];
                        adjacency_matrix[i][node_map[v]] = vweight;
                        int part_v = get_partition(node_map[v],num);
                        if(bm_cnt[part_v] == false){
                            need_process_matrix.push(part_v);
                            bm_cnt[part_v] = true;
                        }
                    }
                    bm.reset(node_map[u]);
                }
            }
//            process_time += timer.Finish();
            while(!need_process_matrix.empty()){
                timer.Start();
                alignas(64) int state[num];
                __m512i state_vector;
                for(int i = 0; i < num; i ++){
                    if(i < vertexs.size()){
                        int u = vertexs[i];
                        u = node_inv_map[u];
                        state[i] = dist[u];
                    }
                    else{
                        state[i] = INF;
                    }
                }
                state_vector = _mm512_load_epi32(state);
                process_time += timer.Finish();
                int part_v = need_process_matrix.front();
                need_process_matrix.pop();
                int start_node = part_v*num;
                for(int i = start_node; i < start_node + num && i < num_node; i ++){
                    int u = node_inv_map[i];
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
                        dist_old[u] = min(dist_old[u],result[j]);
                    }
//                    process_time += timer.Finish();
                    if(dist_old[u] < dist[u]){
                        need_update[u] = true;
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
                    int part = get_partition(node_map[nu], num);
                    if(bm_part.get(part) == 0){
                        q.push(part);
                        bm_part.set(part);
                    }
                    bm.set(node_map[nu]);
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

void sssp_DAG_async_SIMD(long long int &cal_times, float &process_time, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    Timer timer;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
    queue<int> q;
    for(auto i:start){
        bm_part.set(get_partition(node_map[i-start_id], num));
        q.push(get_partition(node_map[i-start_id], num));
    }
    while (!q.empty()) {
        int s = q.size();
        set<int> need_update;
//        vector<bool> need_update(num_node,false);
        queue<int> need_process_matrix;

        int part = q.front();
        q.pop();
        bm_part.reset(part);

        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false);

//        timer.Start();
        vector<vector<int>> adjacency_matrix(num,vector<int>(partition*num,INF));
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);

                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int vweight = neighbors_w[j];
                    adjacency_matrix[i][node_map[v]] = vweight;
                    int part_v = get_partition(node_map[v],num);
                    if(bm_cnt[part_v] == false){
                        need_process_matrix.push(part_v);
                        bm_cnt[part_v] = true;
                    }
                }
                bm.reset(node_map[u]);
            }
        }
//        process_time += timer.Finish();
        while(!need_process_matrix.empty()){
            timer.Start();
            alignas(64) int state[num];
            __m512i state_vector;
            for(int i = 0; i < num; i ++){
                if(i < vertexs.size()){
                    int u = vertexs[i];
                    u = node_inv_map[u];
                    state[i] = dist[u];
                }
                else{
                    state[i] = INF;
                }
            }
            state_vector = _mm512_load_epi32(state);
            process_time += timer.Finish();
            int part_v = need_process_matrix.front();
            need_process_matrix.pop();
            int start_node = part_v*num;
            for(int i = start_node; i < start_node + num && i < num_node; i ++){
                int u = node_inv_map[i];
                alignas(64) int matrix_col[num];
                alignas(64) int result[num];
                timer.Start();
                for(int j = 0; j < num; j ++){
                    matrix_col[j] = adjacency_matrix[j][i];
                }
                __m512i matrix_col_vector = _mm512_load_epi32(matrix_col);
                process_time += timer.Finish();
//                timer.Start();
                __m512i result_vector = _mm512_add_epi32(state_vector,matrix_col_vector);
//                process_time += timer.Finish();
                timer.Start();
                _mm512_store_epi32(result,result_vector);
                process_time += timer.Finish();
//                timer.Start();
                for(int j = 0; j < num; j ++){
                    dist_old[u] = min(dist_old[u],result[j]);
                }
//                process_time += timer.Finish();
                if(dist_old[u] < dist[u]){
                    need_update.insert(u);
//                    need_update[u] = true;
                }
            }
//        timer.Start();
            for(auto nu : need_update){
//            for(int nu = 0; nu < num_node; nu ++){
//                if(need_update[nu]){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part = get_partition(node_map[nu], num);
                        if(bm_part.get(part) == 0){
                            q.push(part);
                            bm_part.set(part);
                        }
                        bm.set(node_map[nu]);
                    }
//                }

            }
//        process_time += timer.Finish();
        }


//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
//#pragma omp parallel for num_threads(1)
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_DAG_SIMD(long long int &cal_times, float &process_time, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
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
                    if(i < vertexs.size()){
                        int u = vertexs[i];
                        state[i] = dist[u];
                    }
                    else{
                        state[i] = INF;
                    }
                }
                state_vector = _mm512_load_epi32(state);
                process_time += timer.Finish();
                int part_v = need_process_matrix.front();
                need_process_matrix.pop();
                int start_node = part_v*num;
                for(int i = start_node; i < start_node + num && i < num_node; i ++){
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

void sssp_diag_iter_priority_SIMD(long long int &cal_times, long long int &useless_cal_times, float &process_time, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    Timer timer;
    unordered_map<int,int> map;//first:part second:priority
    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    priority_queue<pair<int,int>,vector<pair<int,int>>,cmp_sssp> priorityQueue;
    int part_index = node_index / num;  //需要进行块内迭代的块的数量

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
//    queue<int> q;
    for(auto i:start){
        int part = get_partition(node_map[i-start_id], num);
        bm_part.set(part);
        priorityQueue.push(make_pair(part,global_dist[i]));
    }

    queue<int> need_process_matrix;

    while (!priorityQueue.empty()) {
        pair<int,int> top = priorityQueue.top();
        priorityQueue.pop();
        int part = top.first;
        if(bm_part.get(part) == 0)
            continue;

        bm_part.reset(part);
        set<int> need_update;
//        vector<bool> need_update(num_node,false);
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0

//        timer.Start();
        vector<vector<int>> adjacency_matrix(num,vector<int>(partition*num,INF));
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
//            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);

                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int vweight = neighbors_w[j];
                    adjacency_matrix[i][node_map[v]] = vweight;
                    int part_v = get_partition(node_map[v],num);
                    if(bm_cnt[part_v] == false && bm.get(node_map[u]) == 1){
                        need_process_matrix.push(part_v);
                        bm_cnt[part_v] = true;
                    }
                }
//                bm.reset(node_map[u]);
//            }
        }
//        process_time += timer.Finish();

//        //迭代处理16*16块内的状态传递
        bool finished = false;
        bool is_add = false;
        bool is_itr = false;
        int active_num = 0;
        for(int i=0;i<vertexs.size();i++) {
            int u = vertexs[i];
            if (bm.get(u) == 1) {
                active_num ++;
                bm.reset(u);
            }
        }
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                int start_node = part*num;
                alignas(64) int state[num];
                __m512i state_vector;
                timer.Start();
                for(int i = 0; i < num; i ++){
                    if(i < vertexs.size()){
                        int u = vertexs[i];
                        u = node_inv_map[u];
                        state[i] = dist[u];
                    }
                    else{
                        state[i] = INF;
                    }
                }
                state_vector = _mm512_load_epi32(state);
                process_time += timer.Finish();
                for(int i = start_node; i < start_node + num && i < num_node; i ++){
                    cal_times ++;
                    int u = i;
                    if(i < dist.size()){
                        u = node_inv_map[i];
                    }
                    alignas(64) int matrix_col[num];
                    alignas(64) int result[num];
                    timer.Start();
                    bool is_useless = true;
                    for(int j = 0; j < num; j ++){
                        matrix_col[j] = adjacency_matrix[j][i];
                        if(matrix_col[j] != INF)
                            is_useless = false;
                    }
                    if(is_useless){
//                        continue;
                        useless_cal_times ++;
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
                        dist_old[u] = min(dist_old[u],result[j]);
                    }
//                    process_time += timer.Finish();
                    if(dist_old[u] < dist[u]){
//                        need_update[i] = true;
                        need_update.insert(u);
                        finished = false;
                        active_num ++;
                    }
                }
                for(auto nu : need_update){
                    if(dist_old[nu] < dist[nu]){
                        dist[nu] = dist_old[nu];
                        bm.set(node_map[nu]);

                        vector<int> neigh = csr.neighbors(nu);
                        vector<int> neighbors_w = csr.neighbors_weights(nu);

                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int vweight = neighbors_w[j];
                            int part_v = get_partition(node_map[v],num);
                            if(bm_cnt[part_v] == false){
                                need_process_matrix.push(part_v);
                                bm_cnt[part_v] = true;
                            }
                        }
                    }
                }
            }
        }

        while(!need_process_matrix.empty()){
            need_update.clear();
            timer.Start();
            alignas(64) int state[num];
            __m512i state_vector;
            for(int i = 0; i < num; i ++){
                if(i < vertexs.size()){
                    int u = vertexs[i];
                    u = node_inv_map[u];
                    state[i] = dist[u];
                }
                else{
                    state[i] = INF;
                }

            }
            state_vector = _mm512_load_epi32(state);
            process_time += timer.Finish();
            int part_v = need_process_matrix.front();
            need_process_matrix.pop();
            int start_node = part_v*num;
            for(int i = start_node; i < start_node + num && i < num_node; i ++){
                cal_times ++;
                int u = i;
                if(i < dist.size()){
                    u = node_inv_map[i];
                }
                alignas(64) int matrix_col[num];
                alignas(64) int result[num];
                timer.Start();
                bool is_useless = true;
                for(int j = 0; j < num; j ++){
                    matrix_col[j] = adjacency_matrix[j][i];
                    if(matrix_col[j] != INF)
                        is_useless = false;
                }
                if(is_useless){
//                    continue;
                    useless_cal_times ++;
                }
                __m512i matrix_col_vector = _mm512_load_epi32(matrix_col);
                process_time += timer.Finish();
//                timer.Start();
                __m512i result_vector = _mm512_add_epi32(state_vector,matrix_col_vector);
//                process_time += timer.Finish();
                timer.Start();
                _mm512_store_epi32(result,result_vector);
                process_time += timer.Finish();
//                timer.Start();
                for(int j = 0; j < num; j ++){
                    dist_old[u] = min(dist_old[u],result[j]);
                }
//                process_time += timer.Finish();
                if(dist_old[u] < dist[u]){
//                    need_update[i] = true;
                    need_update.insert(u);
                }
            }
            for(auto nu : need_update){
                if(dist[nu] > dist_old[nu]){
                    dist[nu] = min(dist[nu], dist_old[nu]);
                    dist_old[nu] = dist[nu];
                    int part_nu = get_partition(node_map[nu], num);
                    if(bm_part.get(part_nu) == 0){
                        bm_part.set(part_nu);
                    }
                    if(bm.get(node_map[nu]) == 0){

                        priorityQueue.push(make_pair(part_nu,dist[nu]));
                        bm.set(node_map[nu]);
                    }
                }
            }
        }

    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_diag_iter_priority_SIMD_plus(long long int &cal_times, long long int &useless_cal_times, float &process_time, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    Timer timer;
    unordered_map<int,int> map;//first:part second:priority
    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    priority_queue<pair<int,int>,vector<pair<int,int>>,cmp_sssp> priorityQueue;
//    queue<int> Q;
    int part_index = node_index / num;  //需要进行块内迭代的块的数量

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
//    queue<int> q;
    for(auto i:start){
        int part = get_partition(node_map[i-start_id], num);
//        if(bm_part.get(part) == 0){
//            Q.push(part);
//        }
        bm_part.set(part);
        priorityQueue.push(make_pair(part,global_dist[i]));

    }

    queue<int> need_process_matrix;

    while (!priorityQueue.empty()) {
//    while (!Q.empty()) {
//        timer.Start();
        pair<int,int> top = priorityQueue.top();
        priorityQueue.pop();
        int part = top.first;
//        int part = Q.front();
//        Q.pop();
        if(bm_part.get(part) == 0)
            continue;

        bm_part.reset(part);
        set<int> need_update;
//        vector<bool> need_update(num_node,false);
        vector<int> vertexs = vertex_n_partion(part, n, num);
//        vector<bool> bm_cnt(partition,false); //模拟去0

        //迭代处理16*16块内的状态传递
        bool finished = false;
        bool is_add = false;
        bool is_itr = false;
        int active_num = 0;
        for(int i=0;i<vertexs.size();i++) {
            int u = vertexs[i];
            if (bm.get(u) == 1) {
                active_num ++;
            }
        }

//        timer.Start();
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
                            if(part_v == part){
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
                                    cal_times++;
                                    dist[v] = dist_old[v];
                                    bm.set(node_map[v]);
                                    finished = false;
                                    active_num ++;
                                }
                            }

                        }
                    }
                }
            }
        }
//        process_time += timer.Finish();
//        timer.Start();
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];

            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
//                timer.Start();
                alignas(64) int state[num];
                __m512i state_vector;
                for(int j = 0; j < SIMD_SIZE; j ++){
                    state[j] = dist[u];
                }
                state_vector = _mm512_load_epi32(state);
//                process_time += timer.Finish();
                for(int j = 0; j < neigh.size(); j += SIMD_SIZE){
//                    timer.Start();
                    alignas(64) int matrix_col[SIMD_SIZE];
                    alignas(64) int result[SIMD_SIZE];
                    for(int k = j; k < j + SIMD_SIZE && k < neigh.size(); k ++){
                        matrix_col[k-j] = neighbors_w[k];
                    }
                    __m512i matrix_col_vector = _mm512_load_epi32(matrix_col);
//                    process_time += timer.Finish();
                    __m512i result_vector = _mm512_add_epi32(state_vector,matrix_col_vector);
//                    timer.Start();
                    _mm512_store_epi32(result,result_vector);
//                    process_time += timer.Finish();
//                    timer.Start();
                    for(int k = j; k < j + SIMD_SIZE && k < neigh.size(); k ++){
                        dist_old[neigh[k]] = min(dist_old[neigh[k]],result[k-j]);
                        if(dist_old[neigh[k]] < dist[neigh[k]]){
                            cal_times++;
                            dist[neigh[k]] = dist_old[neigh[k]];
                            int part_nu = get_partition(node_map[neigh[k]], num);
                            if(bm_part.get(part_nu) == 0){
                                bm_part.set(part_nu);
//                                Q.push(part_nu);
                            }
                            if(bm.get(node_map[neigh[k]]) == 0){
                                priorityQueue.push(make_pair(part_nu,dist[neigh[k]]));

                                bm.set(node_map[neigh[k]]);
                            }
                        }
                    }
//                    process_time += timer.Finish();
                }
                bm.reset(node_map[u]);
            }

        }
//        process_time += timer.Finish();
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

int main(int argc, char ** argv) {
    string filename(argv[1]);
    int num = stoi(argv[2]);
    Timer timer;
    Timer pre_timer;
    Timer process_timer;
    float process_time = 0;
    Graph<OutEdgeWeighted> graph(filename, true);
    int source_node = 0;    //以重排序前的0顶点为源顶点
    if(argc == 4){
        source_node = stoi(argv[3]);
    }

    vector<int> offset = vector<int>(graph.num_nodes+1,0);
    vector<int> neighbor = vector<int>(graph.num_edges,0);
    vector<int> weight = vector<int>(graph.num_edges,0);

    vector<int> scc = vector<int>(graph.num_nodes,0);
    vector<int> scc_node_num = vector<int>(graph.num_nodes,0);
    vector<int> offset_scc = vector<int>(graph.num_nodes+1,0);
    vector<int> neighbor_scc = vector<int>(graph.num_edges,0);
    vector<int> weight_scc = vector<int>(graph.num_edges,0);
    vector<int> split_scc = vector<int>(graph.num_nodes,0);

    vector<int> global_dist = vector<int>(graph.num_nodes,INF);

    for(int i = 0; i <= graph.num_nodes; i ++)
        offset[i] = int(graph.offset[i]);
    for(int i = 0; i < graph.num_edges; i ++)
        neighbor[i] = int(graph.edgeList[i].end);
    for(int i = 0; i < graph.num_edges; i ++)
        weight[i] = int(graph.edgeList[i].w8);

    timer.Start();
    int scc_num;
    int max_scc_id = 0;
    pre_timer.Start();
    vector<int> map = top_scc(source_node,max_scc_id,graph.num_nodes,graph.num_edges,offset,neighbor,weight,scc_num,scc,scc_node_num,offset_scc,neighbor_scc,weight_scc,split_scc);
    float pre_runtime = pre_timer.Finish();
    float processing_time = 0;
    int node_start_index = 0;     //当前连通分量第一个顶点的id
    long long int cal_times = 0;
    long long int useless_cal_times = 0;
    long long int cal_times_community = 0;
    long long int cal_times_between_communities = 0;
    vector<set<int>> start_node_list(scc_num,set<int>());       //活跃顶点集合
    vector<set<int>> neighbor_node_list(scc_num,set<int>());    //连通分量的邻居顶点集合
    global_dist[map[source_node]] = 0;
    start_node_list[scc[map[source_node]]].insert(map[source_node]);

    cout << "Processing start"<<endl;
    for(int i = 0; i < scc_num; i ++){
        if(start_node_list[i].empty()){
            node_start_index += scc_node_num[i];
            continue;
        }

        if(scc_node_num[i] > 1){
            //构造当前连通分量的csr和csc
            CSR csr;
            CSR csc;
            csr.n = scc_node_num[i];
            csc.n = scc_node_num[i];
            int num_edges = 0;      //当前连通分量内部的边数
            for(int j = 0; j < csr.n; j++){
                int node_id = node_start_index + j;
                num_edges += split_scc[node_id] - offset_scc[node_id];
            }
            csr.row_ptr.assign(csr.n+1, 0);
            csr.col_idx.assign(num_edges, 0);
            csr.weights.assign(num_edges, 0);
            csc.row_ptr.assign(csc.n+1, 0);
            csc.col_idx.assign(num_edges, 0);
            csc.weights.assign(num_edges, 0);
            vector<int> indegree(csc.n,0);
            for(int j = 0; j < csr.n; j ++){
                int node_id = node_start_index + j;
                int outdegree = split_scc[node_id] - offset_scc[node_id];
                csr.row_ptr[j+1] = csr.row_ptr[j] + outdegree;
                for(int k = 0; k < outdegree; k ++){
                    csr.col_idx[csr.row_ptr[j]+k] = neighbor_scc[offset_scc[node_id]+k] - node_start_index;
                    csr.weights[csr.row_ptr[j]+k] = weight_scc[offset_scc[node_id]+k];
                    indegree[neighbor_scc[offset_scc[node_id]+k] - node_start_index] ++;
                }
            }

            vector<int> csc_index(csc.n,0);
            for(int j = 0; j < csc.n; j ++){
                csc.row_ptr[j+1] = csc.row_ptr[j] + indegree[j];
            }
            for(int j = 0; j < csc.n; j ++){
                int node_id = node_start_index + j;
                int outdegree = split_scc[node_id] - offset_scc[node_id];
                for(int k = 0; k < outdegree; k ++){
                    int index = neighbor_scc[offset_scc[node_id]+k] - node_start_index;     //计算邻居顶点在csc对应的位置
                    csc.col_idx[csc.row_ptr[index] + csc_index[index]] = j;
                    csc.weights[csc.row_ptr[index] + csc_index[index]] = weight_scc[offset_scc[node_id]+k];
                    csc_index[index] ++;
                }
            }
            //顶点重排序，稠密化邻接数组
            int row = 0;
            int col = 0;
            vector<bool> row_visit(csr.n, false);
            vector<bool> col_visit(csr.n, false);
            vector<int> row_map(csr.n,0);   //重排序前 -> 重排序后
            vector<int> col_map(csr.n,0);
            vector<int> row_inv_map(csr.n,0);        //重排序后 -> 重排序前
            vector<int> col_inv_map(csr.n,0);
            set<int> row_set;
            set<int> col_set;
            //比较最大强连通分量中社区内元素稠密度和重排序前的稠密度
            if(i == max_scc_id){
                cout<<"vertex num in largest scc: "<<scc_node_num[i]<<endl;
            }
            //16-size-blocK顶点重排序
            pre_timer.Start();
            vector<vector<int>> map_degree = reorder_degree(csr,csc,num);
            pre_runtime += pre_timer.Finish();
            vector<int> node_degree_inv_map = map_degree[0];        //重排序后的顶点id -> 原顶点id
            vector<int> node_degree_map = map_degree[1];            //原顶点id -> 重排序后的顶点id
            int node_index = map_degree[2][0];
            //当前连通分量内部的状态传递
            process_timer.Start();
//            sssp_diag_iter_priority_SIMD(cal_times,useless_cal_times,process_time,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_dist,node_start_index,scc_node_num[i],num);
            sssp_diag_iter_priority_SIMD_plus(cal_times,useless_cal_times,process_time,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_dist,node_start_index,scc_node_num[i],num);
//            sssp_DAG_sync_SIMD(cal_times,process_time,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_dist,node_start_index,scc_node_num[i],num);
//            sssp_DAG_async_SIMD(cal_times,process_time,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_dist,node_start_index,scc_node_num[i],num);
//            sssp_DAG_SIMD(cal_times,process_time,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_dist,node_start_index,scc_node_num[i],num);
            processing_time += process_timer.Finish();
        }

//        if(i != max_scc_id)
//            cal_times = temp;
//        cout<<i<<" "<<start_node_list[i].size()<<endl;
        //连通分量之间的状态传递
//        long long int cal_size = 0;   //连通分量间计算的矩阵维度
        process_timer.Start();
        long long int part_num = 0;
//        int row_part_num = (csr.n+num-1)/num;
        vector<bool> visit((graph.num_nodes+num-1)/num,false);
        for(int j = 0; j < scc_node_num[i]; j ++){
            int node_id = node_start_index + j;
            if(j != 0 && (j/num != (j-1)/num)){
                visit.assign(visit.size(),0);
            }
//            if(global_dist[node_id] == INF)
//                continue;
            int outdegree = offset_scc[node_id+1] - split_scc[node_id];
            for(int k = 0; k < outdegree; k ++){
                int neighbor_id = neighbor_scc[split_scc[node_id]+k];
                if(visit[neighbor_id/num] == false){
                    visit[neighbor_id/num] = true;
                    part_num ++;
                }
                if(global_dist[node_id] + weight_scc[split_scc[node_id]+k] < global_dist[neighbor_id]){
                    global_dist[neighbor_id] = global_dist[node_id] + weight_scc[split_scc[node_id]+k];
                    start_node_list[scc[neighbor_id]].insert(neighbor_id);
                    cal_times ++;
                }
//                if(neighbor_node_list[i].count(neighbor_id) == 0){
//                    neighbor_node_list[i].insert(neighbor_id);
//                    cal_size += scc_node_num[scc[neighbor_id]];
//                }
            }
        }
//        cal_times += scc_node_num[i] * cal_size;

//        cal_times += part_num * num * num;
        node_start_index += scc_node_num[i];
        processing_time += process_timer.Finish();
    }
    float runtime = timer.Finish();
    cout << "Total runtime: " << runtime/1000 << " (s).\n";
    cout << "Pre-processing finished in " << pre_runtime/1000 << " (s).\n";
    cout << "Processing finished in " << processing_time/1000 << " (s).\n";
    cout << "Total processing time: " << (pre_runtime+processing_time)/1000 << " (s).\n";
    cout << "Calculating finished in " << process_time/1000 << " (s).\n";
    cout<<"calculation times: "<<cal_times<<endl;
    cout<<"the ratio of useless calculation times: "<<double(useless_cal_times)/cal_times<<endl;
//    cout<<"Proportion of calculation outside the community: "<<double(cal_times_between_communities)/cal_times_community<<endl;
    int update_num = 0;
    for(int i = 0; i < graph.num_nodes; i ++){
        if(global_dist[i] != INF)
            update_num ++;
    }
    cout<<"updated vertex num: "<<update_num<<endl;
    for(int i = 0; i < min(100,(int)graph.num_nodes); i ++){
        cout<<global_dist[map[i]]<<" ";
    }
    cout<<endl;


    return 0;
}