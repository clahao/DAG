//
// Created by moyu on 2023/5/15.
//
#ifndef SSWP_PARTITION
#define SSWP_PARTITION
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
struct cmp_sswp{
    bool operator()(pair<int,int> &a, pair<int,int> &b){
        if(a.second != b.second)
            return a.second < b.second;
        else
            return a.first < b.first;
    }
};

void sswp_diag_iter_priority(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    unordered_map<int,int> map;//first:part second:priority
    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    priority_queue<pair<int,int>,vector<pair<int,int>>,cmp_sswp> priorityQueue;
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

//    vector<bool> is_inserted(partition, false);
//    for(auto i:start){
//        int part = get_partition(node_map[i-start_id], num);
//        if(!is_inserted[part]){
//            is_inserted[part] = true;
//            vector<int> vertexs = vertex_n_partion(part, n, num);
//            double priority = 0;
//            int active_num = 0;
//            for(int j=0;j<vertexs.size();j++) {
//                int u = vertexs[j];
//                if (bm.get(u) == 1) {
//                    active_num ++;
//                    priority += dist[node_inv_map[u]];
//                }
//            }
//            priorityQueue.push(make_pair(part,priority/active_num));
//        }
//    }

//    bool is_finished = false;
    while (!priorityQueue.empty()) {
//        int part = q.front();
//        q.pop();

//        is_finished = true;
//        int priority = 10000;
//        int part;
//        for(auto it = map.begin(); it != map.end(); it ++){
//            if((*it).second >= 0 && (*it).second < priority ){
//                priority = (*it).second;
//                part = (*it).first;
//                is_finished = false;
//            }
//        }
        pair<int,int> top = priorityQueue.top();
        priorityQueue.pop();
        int part = top.first;
        if(bm_part.get(part) == 0)
            continue;
//        if(is_finished)
//            break;
//        map.erase(part);
//        map[part] = 0;
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
                                dist_old[v] = max(dist_old[v], min(vweight,dist[u]));
                                if(dist_old[v] > dist[v]){
                                    need_update.insert(v);
//                                    dist[v] = dist_old[v];
//                                    bm.set(node_map[v]);
                                    finished = false;
                                    active_num ++;
                                }
                            }

                        }
                    }
                }
                for(auto nu : need_update){
                    if(dist_old[nu] > dist[nu]){
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
                    dist_old[v] = max(dist_old[v], min(vweight,dist[u]));
                    if(dist_old[v] > dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
            if(dist[nu] < dist_old[nu]){
                dist[nu] = max(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                int part_nu = get_partition(node_map[nu], num);
                if(bm_part.get(part_nu) == 0){
//                    q.push(part_nu);
                    bm_part.set(part_nu);
                }
                if(bm.get(node_map[nu]) == 0){
//                    vector<int> vertexs = vertex_n_partion(part_nu, n, num);
//                    double priority = 0;
//                    int node_num = 0;
//                    for(int j=0;j<vertexs.size();j++) {
//                        int u = vertexs[j];
//                        if (bm.get(u) == 1) {
//                            node_num ++;
//                            priority += dist[node_inv_map[u]];
//                        }
//                    }
//                    priorityQueue.push(make_pair(part_nu,priority/node_num));
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

void sswp_DAG(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
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
                        dist_old[v] = max(dist_old[v], min(vweight,dist[u]));
                        if(dist_old[v] > dist[v]){
                            need_update.insert(v);
                        }
                    }
                    bm.reset(u);
                }
            }
        }

        for(auto nu : need_update){
            if(dist[nu] < dist_old[nu]){
                dist[nu] = max(dist[nu], dist_old[nu]);
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

void sswp_sync(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
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
                        dist_old[v] = max(dist_old[v], min(vweight,dist[u]));
                        if(dist_old[v] > dist[v]){
                            need_update.insert(v);
                        }
                    }
                    bm.reset(node_map[u]);
                }
            }
        }

        for(auto nu : need_update){
            if(dist[nu] < dist_old[nu]){
                dist[nu] = max(dist[nu], dist_old[nu]);
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

void sswp_async(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<int> &global_dist, int start_id, int num_node ,int num){
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
                    dist_old[v] = max(dist_old[v], min(vweight,dist[u]));
                    if(dist_old[v] > dist[v]){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
            }
        }
        for(auto nu : need_update){
            if(dist[nu] < dist_old[nu]){
                dist[nu] = max(dist[nu], dist_old[nu]);
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

#endif