#ifndef SSSP_PARTITION
#define SSSP_PARTITION
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <queue>
#include <limits>
#include <cassert>
#include <iomanip>
#include "Graph.h"
#include "Graph.cpp"


// 记录迭代轮次
void sssp_reoeder(long long int &cal_times, const CSR& csr, set<int> start, vector<int> row_map, vector<int> col_map, vector<int> row_inv_map, vector<int> col_inv_map,vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    Bitmap bm(n);
    for(auto i:start){
        bm.set(col_map[i-start_id]);
    }


    Bitmap bm_part(partition);
    queue<int> q;
    for(auto i:start){
        bm_part.set(get_partition(col_map[i-start_id], num));
        q.push(get_partition(col_map[i-start_id], num));
    }
    while (!q.empty()) {
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = col_inv_map[u];//原本的id
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(row_map[v],num);
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(col_map[u]);
            }
        }
        for(auto nu : need_update){
            if(dist[nu] > dist_old[nu]){
                dist[nu] = min(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                int part = get_partition(col_map[nu], num);
                if(bm_part.get(part) == 0){
                    q.push(part);
                    bm_part.set(part);
                }
                bm.set(col_map[nu]);
            }
        }
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_rabbit(long long int &cal_times, long long int &cal_times_community, long long int &cal_times_between_communities, const CSR& csr, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
//                        cal_times_community += num * num;
//                        if(community[u] != community[v]){
//                            cal_times ++;
//                            cal_times_between_communities += num * num;
//                        }
//                        else{
//                            cal_times += num * num;
//                        }

                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
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
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_rabbit_push_pull_col(long long int &cal_times, long long int &cal_times_community, long long int &cal_times_between_communities, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        queue<int> q_copy;
        Bitmap bm_update(partition);
        int q_size = q.size();
        for(int itr = 0; itr < q_size; itr ++){
            int part = q.front();
            q_copy.push(part);
            q.pop();
            vector<int> vertexs = vertex_n_partion(part, n, num);
            vector<bool> bm_cnt(partition,false); //模拟去0
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                u = node_inv_map[u];
                vector<int> neigh = csc.neighbors(u);
                vector<int> neighbors_w = csc.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
//                    if(bm_part.get(part_v) == 1){//
                    if(bm.get(node_map[v]) == 1){//
                        if(bm_cnt[part_v] == false){
//                            cal_times += num * num;
                            cal_times_community += num * num;
                            bm_cnt[part_v] = true;
                        }
                        int vweight = neighbors_w[j];
                        dist_old[u] = min(dist_old[u], vweight + dist[v]);
                        if(dist_old[u] < dist[u]){
                            dist[u] = dist_old[u];
                            bm.set(node_map[u]);
//                            if(part_v >= part)
                                bm_update.set(part);
                        }
                    }
                }

            }
        }
        q = q_copy;
        Bitmap bm_part_copy = bm_part;
        for(int itr = 0; itr < q_size; itr ++){
            int part = q.front();
            q.pop();
            bm_part.reset(part);
            set<int> need_update;
            vector<int> vertexs = vertex_n_partion(part, n, num);
            vector<bool> bm_cnt(partition,false); //模拟去0
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                if(bm.get(u) == 1){
                    u = node_inv_map[u];
                    vector<int> neigh = csr.neighbors(u);
                    vector<int> neighbors_w = csr.neighbors_weights(u);
                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int part_v = get_partition(node_map[v],num);
                        if(bm_update.get(part) == 1 || bm_part_copy.get(part_v) == 0){
                            if(bm_cnt[part_v] == false){
                                cal_times += num * num;
                                bm_cnt[part_v] = true;
                            }
                            int vweight = neighbors_w[j];
                            dist_old[v] = min(dist_old[v], vweight + dist[u]);
                            if(dist_old[v] < dist[v]){
                                need_update.insert(v);
                            }
                        }
                    }
                    bm.reset(node_map[u]);
                }
            }
            for(auto nu : need_update){
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

    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_rabbit_push_pull_row(long long int &cal_times, long long int &cal_times_community, long long int &cal_times_between_communities, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
//先处理活跃块之间的边，再更新其他边，在更新其他边的过程中，如果i中有新的顶点被j激活，本次迭代不处理
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        if(q.size() > 0){
            queue<int> q_copy = q;
//            queue<int> update_part;
            Bitmap bm_part_copy = bm_part;
            Bitmap bm_update(partition);
            Bitmap bm_node_update(n);
            int q_size = q.size();
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                bm_part_copy.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1) {
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for (int j = 0; j < neigh.size(); j++) {
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v], num);
                            if (bm_part.get(part_v) == 1) {//
//                            if(bm.get(node_map[v]) == 1){//
                                if (bm_cnt[part_v] == false) {
                                    cal_times += num * num;
//                                    cal_times_community += num * num;
                                    bm_cnt[part_v] = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if (dist_old[v] < dist[v]) {
                                    need_update.insert(v);
                                }
                            }
                        }
                    }
                }
                for(auto nu : need_update){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part_nu = get_partition(node_map[nu], num);
                        if(bm_update.get(part_nu) == 0 && bm_part_copy.get(part_nu) == 0){//&& part_nu <= part){////
                            bm_update.set(part_nu);     //当前块i被之后的块j更新
                            q_copy.push(part_nu);
//                            update_part.push(part_nu);
                        }
                        if(bm_part_copy.get(part_nu) == 0)
                            bm_node_update.set(node_map[nu]);
                        bm.set(node_map[nu]);
                    }
                }
            }
            bm_part_copy = bm_part;
            q = q_copy;
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                if(bm_update.get(part) == 0)    //bm_update为1表示当前块i已经再次被插入队列中，因此bm_part不需要重置，防止一个块被多次插入队列
                    bm_part.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1 && bm_node_update.get(u) == 0){               //important!!!
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
//                        if(bm_update.get(part) == 1 || bm_part_copy.get(part_v) == 0){      //当前块所在行被更新或者块所对应的列还未处理过
                            if(bm_part_copy.get(part_v) == 0){      //当前块所对应的列还未处理过
                                if(bm_cnt[part_v] == false){
                                    cal_times += num * num;
                                    bm_cnt[part_v] = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
                                    need_update.insert(v);
                                }
                            }
                        }
                        if(bm_node_update.get(node_map[u]) == 0)
                            bm.reset(node_map[u]);
                    }
                }
                for(auto nu : need_update){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part_nu = get_partition(node_map[nu], num);
                        if(bm_part.get(part_nu) == 0) {
                            q.push(part_nu);
                            bm_part.set(part_nu);
                        }
//                   else if(bm_part_copy.get(part_nu) == 1 && bm_update.get(part_nu) == 0){////
//                        bm_update.set(part_nu);
////                        q.push(part);
//                    }
                        bm.set(node_map[nu]);
                    }
                }
            }
        }
        else{
            int q_size = q.size();
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                bm_part.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
                            if(bm_cnt[part_v] == false){
                                cal_times += num * num;
                                bm_cnt[part_v] = true;
                            }
                            int vweight = neighbors_w[j];
                            dist_old[v] = min(dist_old[v], vweight + dist[u]);
                            if(dist_old[v] < dist[v]){
                                need_update.insert(v);
                            }
                        }
                        bm.reset(node_map[u]);
                    }
                }
                for(auto nu : need_update){
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
        }

    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_rabbit_push_pull_16size(long long int &cal_times, long long int &cal_times_community, long long int &cal_times_between_communities, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        if(q.size() > 0){
            queue<int> q_copy = q;
//            queue<int> update_part;
            Bitmap bm_part_copy = bm_part;
            Bitmap bm_update(partition);
            Bitmap bm_node_update(n);
            int q_size = q.size();
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                bm_part_copy.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
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
                if(active_num > 4){//} && part < part_index){
                    is_itr = true;
                    while(!finished && active_num > 0){
                        active_num = 0;
                        finished = true;
                        is_add = false;
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
                                        if(!is_add){
                                            cal_times += num * num;
                                            is_add = true;
                                        }
                                        int vweight = neighbors_w[j];
                                        dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                        if(dist_old[v] < dist[v]){
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
                //处理活跃块之间的边
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1) {
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for (int j = 0; j < neigh.size(); j++) {
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v], num);
                            if(is_itr && part_v == part)
                                continue;
                            if (bm_part.get(part_v) == 1) {//
//                            if(bm.get(node_map[v]) == 1){//
                                if (bm_cnt[part_v] == false) {
                                    cal_times += num * num;
//                                    cal_times_community += num * num;
                                    bm_cnt[part_v] = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if (dist_old[v] < dist[v]) {
                                    need_update.insert(v);
                                }
                            }
                        }
                    }
                }
                for(auto nu : need_update){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part_nu = get_partition(node_map[nu], num);
                        if(bm_update.get(part_nu) == 0 && bm_part_copy.get(part_nu) == 0){//&& part_nu <= part){////
                            bm_update.set(part_nu);     //当前块i被之后的块j更新
                            q_copy.push(part_nu);
//                            update_part.push(part_nu);
                        }
                        if(bm_part_copy.get(part_nu) == 0)
                            bm_node_update.set(node_map[nu]);
                        bm.set(node_map[nu]);
                    }
                }
            }
            bm_part_copy = bm_part;
            q = q_copy;
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                if(bm_update.get(part) == 0)        //bm_update为1表示当前块i已经再次被插入队列中，因此bm_part不需要重置，防止一个块被多次插入队列
                    bm_part.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1 && bm_node_update.get(u) == 0){               //important!!!
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
//                        if(bm_update.get(part) == 1 || bm_part_copy.get(part_v) == 0){      //当前块所在行被更新或者块所对应的列还未处理过
                            if(bm_part_copy.get(part_v) == 0){      //当前块所对应的列还未处理过
                                if(bm_cnt[part_v] == false){
                                    cal_times += num * num;
                                    bm_cnt[part_v] = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
                                    need_update.insert(v);
                                }
                            }
                        }
                        if(bm_node_update.get(node_map[u]) == 0)
                            bm.reset(node_map[u]);
                    }
                }
                for(auto nu : need_update){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part_nu = get_partition(node_map[nu], num);
                        if(bm_part.get(part_nu) == 0) {
                            q.push(part_nu);
                            bm_part.set(part_nu);
                        }
//                   else if(bm_part_copy.get(part_nu) == 1 && bm_update.get(part_nu) == 0){////
//                        bm_update.set(part_nu);
////                        q.push(part);
//                    }
                        bm.set(node_map[nu]);
                    }
                }
            }
        }
        else{
            int q_size = q.size();
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                bm_part.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
                            if(bm_cnt[part_v] == false){
                                cal_times += num * num;
                                bm_cnt[part_v] = true;
                            }
                            int vweight = neighbors_w[j];
                            dist_old[v] = min(dist_old[v], vweight + dist[u]);
                            if(dist_old[v] < dist[v]){
                                need_update.insert(v);
                            }
                        }
                        bm.reset(node_map[u]);
                    }
                }
                for(auto nu : need_update){
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
        }

    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}


void sssp_rabbit_diag_iter(long long int &cal_times, long long int &cal_times_community, long long int &cal_times_between_communities, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
        //迭代处理16*16块内的状态传递
        bool finished = false;
        bool is_add = false;
        int active_num = 0;
        for(int i=0;i<vertexs.size();i++) {
            int u = vertexs[i];
            if (bm.get(u) == 1) {
                active_num ++;
            }
        }
        if(active_num > 8){
            while(!finished && active_num > 8){
                active_num = 0;
                finished = true;
                is_add = false;
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
                                if(!is_add){
                                    cal_times += num * num;
                                    is_add = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
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


        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
//                        cal_times_community += num * num;
//                        if(community[u] != community[v]){
//                            cal_times ++;
//                            cal_times_between_communities += num * num;
//                        }
//                        else{
//                            cal_times += num * num;
//                        }

                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
            if(dist[nu] > dist_old[nu]){
                dist[nu] = min(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                int part_nu = get_partition(node_map[nu], num);
                if(bm_part.get(part_nu) == 0){
                    q.push(part_nu);
                    bm_part.set(part_nu);
                }
                bm.set(node_map[nu]);
            }
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}


void sssp_diag_iter(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    int part_index = node_index / num;  //需要进行块内迭代的块的数量

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
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
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
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                is_add = false;
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
                                if(!is_add){
                                    cal_times += num * num;
                                    is_add = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
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


        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(is_itr && part_v == part)
                        continue;
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
            if(dist[nu] > dist_old[nu]){
                dist[nu] = min(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                int part_nu = get_partition(node_map[nu], num);
                if(bm_part.get(part_nu) == 0){
                    q.push(part_nu);
                    bm_part.set(part_nu);
                }
                bm.set(node_map[nu]);
            }
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

struct cmp{
    bool operator()(pair<int,int> &a, pair<int,int> &b){
        if(a.second != b.second)
            return a.second > b.second;
        else
            return a.first < b.first;
    }
};

void sssp_diag_iter_priority(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    unordered_map<int,int> map;//first:part second:priority
    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    priority_queue<pair<int,int>,vector<pair<int,int>>,cmp> priorityQueue;
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


//    bool is_finished = false;
    while (!priorityQueue.empty()) {
        pair<int,int> top = priorityQueue.top();
        priorityQueue.pop();
        int part = top.first;
        if(bm_part.get(part) == 0)
            continue;

        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
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
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                is_add = false;
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
                                if(!is_add){
                                    cal_times += num * num;
                                    is_add = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
                                    need_update.insert(v);
//                                    dist[v] = dist_old[v];
//                                    bm.set(node_map[v]);
                                    finished = false;
                                    active_num ++;
//                                    cal_times ++;
                                }
                            }

                        }
                    }
                }
                for(auto nu : need_update){
                    if(dist_old[nu] < dist[nu]){
                        dist[nu] = dist_old[nu];
                        bm.set(node_map[nu]);
                    }
                }
            }
        }

        need_update.clear();

        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(is_itr && part_v == part)
                        continue;
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
//                        cal_times ++;
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
            if(dist[nu] > dist_old[nu]){
                dist[nu] = min(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                int part_nu = get_partition(node_map[nu], num);
                if(bm_part.get(part_nu) == 0){
//                    q.push(part_nu);
                    bm_part.set(part_nu);
                }
                if(bm.get(node_map[nu]) == 0){
                    priorityQueue.push(make_pair(part_nu,dist[nu]));
                    bm.set(node_map[nu]);
                }
            }
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_rabbit_push_pull_16size_priority(long long int &cal_times, long long int &cal_times_community, long long int &cal_times_between_communities, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> community, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> part_priority(partition,INF);   //每个块的优先级
    priority_queue<pair<int,int>,vector<pair<int,int>>,cmp> priorityQueue;

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
    queue<int> q;
    for(auto i:start){
        int part = get_partition(node_map[i-start_id], num);
        bm_part.set(part);
        q.push(part);
        if(part_priority[part] < global_dist[i]){
            part_priority[part] = global_dist[i];
        }
    }
    while (!q.empty()) {
        if(q.size() > 0){
            queue<int> q_copy = q;
//            queue<int> update_part;
            Bitmap bm_part_copy = bm_part;
            Bitmap bm_update(partition);
            Bitmap bm_node_update(n);
            int q_size = q.size();
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                priorityQueue.push(make_pair(part,part_priority[part]));
                part_priority[part] = INF;
            }
            for(int itr = 0; itr < q_size; itr ++){
//                int part = q.front();
//                q.pop();
                int part = priorityQueue.top().first;
                priorityQueue.pop();
                bm_part_copy.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
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
                if(active_num > 4){//} && part < part_index){
                    is_itr = true;
                    while(!finished && active_num > 0){
                        active_num = 0;
                        finished = true;
                        is_add = false;
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
                                        if(!is_add){
                                            cal_times += num * num;
                                            is_add = true;
                                        }
                                        int vweight = neighbors_w[j];
                                        dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                        if(dist_old[v] < dist[v]){
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
                //处理活跃块之间的边
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1) {
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for (int j = 0; j < neigh.size(); j++) {
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v], num);
                            if(is_itr && part_v == part)
                                continue;
                            if (bm_part.get(part_v) == 1) {//
//                            if(bm.get(node_map[v]) == 1){//
                                if (bm_cnt[part_v] == false) {
                                    cal_times += num * num;
//                                    cal_times_community += num * num;
                                    bm_cnt[part_v] = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if (dist_old[v] < dist[v]) {
                                    need_update.insert(v);
                                }
                            }
                        }
                    }
                }
                for(auto nu : need_update){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part_nu = get_partition(node_map[nu], num);
                        if(bm_update.get(part_nu) == 0 && bm_part_copy.get(part_nu) == 0){//&& part_nu <= part){////
                            bm_update.set(part_nu);     //当前块i被之后的块j更新
                            q_copy.push(part_nu);
                            part_priority[part_nu] = min(part_priority[part_nu],dist[nu]);  //更新块优先级
//                            update_part.push(part_nu);
                        }
                        if(bm_part_copy.get(part_nu) == 0)
                            bm_node_update.set(node_map[nu]);
                        bm.set(node_map[nu]);
                    }
                }
            }
            bm_part_copy = bm_part;
            q = q_copy;
            for(int itr = 0; itr < q_size; itr ++){
                int part = q.front();
                q.pop();
                if(bm_update.get(part) == 0)        //bm_update为1表示当前块i已经再次被插入队列中，因此bm_part不需要重置，防止一个块被多次插入队列
                    bm_part.reset(part);
                set<int> need_update;
                vector<int> vertexs = vertex_n_partion(part, n, num);
                vector<bool> bm_cnt(partition,false); //模拟去0
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1 && bm_node_update.get(u) == 0){               //important!!!
                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
//                        if(bm_update.get(part) == 1 || bm_part_copy.get(part_v) == 0){      //当前块所在行被更新或者块所对应的列还未处理过
                            if(bm_part_copy.get(part_v) == 0){      //当前块所对应的列还未处理过
                                if(bm_cnt[part_v] == false){
                                    cal_times += num * num;
                                    bm_cnt[part_v] = true;
                                }
                                int vweight = neighbors_w[j];
                                dist_old[v] = min(dist_old[v], vweight + dist[u]);
                                if(dist_old[v] < dist[v]){
                                    need_update.insert(v);
                                }
                            }
                        }
                        if(bm_node_update.get(node_map[u]) == 0)
                            bm.reset(node_map[u]);
                    }
                }
                for(auto nu : need_update){
                    if(dist[nu] > dist_old[nu]){
                        dist[nu] = min(dist[nu], dist_old[nu]);
                        dist_old[nu] = dist[nu];
                        int part_nu = get_partition(node_map[nu], num);
                        if(bm_part.get(part_nu) == 0) {
                            q.push(part_nu);
                            bm_part.set(part_nu);
                            part_priority[part_nu] = min(part_priority[part_nu],dist[nu]);  //更新块优先级
                        }
//                   else if(bm_part_copy.get(part_nu) == 1 && bm_update.get(part_nu) == 0){////
//                        bm_update.set(part_nu);
////                        q.push(part);
//                    }
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

void sssp(long long int &cal_times, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num, int &sum0, int &sum1){
    int num0 = 0;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
        num0 = 0;
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                num0++;
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(v,num);
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(u);
            }
        }
        if(num0 == 1){
            sum0++;
        }
        else{
            sum1++;
        }
        for(auto nu : need_update){
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
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_DAG(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        set<int> need_update;
        while(s--){
            int part = q.front();
            q.pop();
            bm_part.reset(part);

            vector<int> vertexs = vertex_n_partion(part, n, num);
            vector<bool> bm_cnt(partition,false); //模拟去0
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                if(bm.get(u) == 1){
                    vector<int> neigh = csr.neighbors(u);
                    vector<int> neighbors_w = csr.neighbors_weights(u);
                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int part_v = get_partition(v,num);
                        if(bm_cnt[part_v] == false){
                            cal_times += num * num;
                            bm_cnt[part_v] = true;
                        }
                        int vweight = neighbors_w[j];
                        dist_old[v] = min(dist_old[v], vweight + dist[u]);
                        if(dist_old[v] < dist[v]){
                            need_update.insert(v);
                        }
                    }
                    bm.reset(u);
                }
            }
        }

        for(auto nu : need_update){
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
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }
}

void sssp_sync(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

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
        while(s--){
            int part = q.front();
            q.pop();
            bm_part.reset(part);

            vector<int> vertexs = vertex_n_partion(part, n, num);
            vector<bool> bm_cnt(partition,false); //模拟去0
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                if(bm.get(u) == 1){
                    u = node_inv_map[u];
                    vector<int> neigh = csr.neighbors(u);
                    vector<int> neighbors_w = csr.neighbors_weights(u);
                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int part_v = get_partition(v,num);
                        if(bm_cnt[part_v] == false){
                            cal_times += num * num;
                            bm_cnt[part_v] = true;
                        }
                        int vweight = neighbors_w[j];
                        dist_old[v] = min(dist_old[v], vweight + dist[u]);
                        if(dist_old[v] < dist[v]){
                            need_update.insert(v);
                        }
                    }
                    bm.reset(node_map[u]);
                }
            }
        }

        for(auto nu : need_update){
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
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }
}

void sssp_async(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
    queue<int> q;
    for(auto i:start){
        int part = get_partition(node_map[i-start_id], num);
        if(bm_part.get(part) == 0){
            bm_part.set(part);
            q.push(part);
        }
    }
    while (!q.empty()) {
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
        vector<bool> bm_cnt(partition,false); //模拟去0
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                vector<int> neigh = csr.neighbors(u);
                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(v,num);
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    int vweight = neighbors_w[j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
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
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }
}

//int main(int argc, char **argv){
//    string filename(argv[1]);
//    int num = stoi(argv[2]);
//    Graph<OutEdgeWeighted> graph(filename, true);
//
//    CSR csr;
//    csr.n = graph.num_nodes;
//    csr.row_ptr.assign(csr.n+1, 0);
//    csr.col_idx.assign(graph.num_edges, 0);
//    csr.weights.assign(graph.num_edges, 0);
//    for (int i = 0; i < csr.row_ptr.size(); i++) {
//        csr.row_ptr[i] = graph.offset[i];
//    }
//    for (int i = 0; i < csr.col_idx.size(); i++) {
//        csr.col_idx[i] = graph.edgeList[i].end;
//    }
//    for (int i = 0; i < csr.weights.size(); i++) {
//        csr.weights[i] = graph.edgeList[i].w8;
//    }
//
//    long long int iterations = 0;
//    set<int> start;
//    start.insert(0);
//    vector<int> value(graph.num_nodes,INF);
//    value[0] = 0;
//    sssp(iterations,csr,start,value,0,graph.num_nodes,num);
//    cout<< "iterations " << left << setw(2) <<iterations << " :" ;
//    for(int i = 0; i < min(100,csr.n); i ++)
//        cout<<value[i]<<" ";
//    cout<<endl;
//#ifdef test_value
//
//    for(int i=0;i<min(100,csr.n);i++){
//      cout << left << setw(2)<< ans[i] << " ";
//    }
//    cout << endl;
//#endif
//    return 0;
//}
#endif