//
// Created by moyu on 2023/5/16.
//
#ifndef PR_PARTITION
#define PR_PARTITION
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
struct cmp_pr{
    bool operator()(pair<int,double> &a, pair<int,double> &b){
        if(a.second != b.second)
            return a.second < b.second;
        else
            return a.first < b.first;
    }
};

void pr_diag_iter(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<double> &global_delta, vector<double> &global_value, vector<int> outDegree, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    unordered_map<int,int> map;//first:part second:priority
    vector<double> delta(global_delta.begin()+start_id,global_delta.begin()+start_id+num_node);
    vector<double> delta_old = vector<double>(num_node,0);
    vector<double> value(global_value.begin()+start_id,global_value.begin()+start_id+num_node);
    int part_index = node_index / num;  //需要进行块内迭代的块的数量

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
        vector<double> delta_temp(num,0);
        vector<double> value_temp(num,0);
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            int start_node = 0;
            for(int i = 0; i < vertexs.size(); i ++){
                int u = vertexs[i];
                if(i == 0)
                    start_node = u;
                u = node_inv_map[u];
                delta_temp[i] = delta[u];
            }
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                is_add = false;
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        value_temp[i] += delta_temp[i];
                        vector<int> neigh = csr.neighbors(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
                            if(part_v == part){
                                if(!is_add){
                                    cal_times += num * num;
                                    is_add = true;
                                }
                                delta_old[v] += delta_temp[i] * damping_factor / outDegree[u];
                                if(delta_old[v] > threshold){
                                    need_update.insert(v);
//                                    delta[v] += delta_old[v];
//                                    delta_old[v] = 0;
//                                    bm.set(node_map[v]);
                                    finished = false;
                                    active_num ++;
                                }
                            }

                        }
                        delta_temp[i] = 0;
                    }

                }
                for(auto nu : need_update){
                    if(delta_old[nu] > threshold){
                        delta_temp[node_map[nu]-start_node] += delta_old[nu];
                        delta_old[nu] = 0;
                        bm.set(node_map[nu]);
                    }
                }
            }
        }

        if(is_itr){
            for(int i = 0; i < vertexs.size(); i ++){
                int u = vertexs[i];
                u = node_inv_map[u];
                delta[u] = value_temp[i];
            }
        }

        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                value[u] += delta[u];
                vector<int> neigh = csr.neighbors(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(is_itr && part_v == part)
                        continue;
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    delta_old[v] += delta[u] * damping_factor / outDegree[u];
                    if((delta_old[v] + delta[v] > threshold && part != part_v) || (delta_old[v] > threshold && part == part_v)){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
                delta[u] = 0;
            }

        }
        for(auto nu : need_update){
//            if(delta_old[nu] + delta[nu] > threshold){
            delta[nu] += delta_old[nu];
            delta_old[nu] = 0;
            int part_nu = get_partition(node_map[nu], num);
            if(bm_part.get(part_nu) == 0){
                q.push(part_nu);
                bm_part.set(part_nu);
            }
            if(bm.get(node_map[nu]) == 0){
                bm.set(node_map[nu]);
            }
//            }
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_value[start_id+i] = value[i] + delta_old[i];
    }

}

void pr_diag_iter_priority(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<double> &global_delta, vector<double> &global_value, vector<int> outDegree, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    unordered_map<int,int> map;//first:part second:priority
    vector<double> delta(global_delta.begin()+start_id,global_delta.begin()+start_id+num_node);
    vector<double> delta_old = vector<double>(num_node,0);
    vector<double> value(global_value.begin()+start_id,global_value.begin()+start_id+num_node);
    priority_queue<pair<int,double>,vector<pair<int,double>>,cmp_pr> priorityQueue;
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
        priorityQueue.push(make_pair(part,global_delta[i]));
    }


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
        vector<double> delta_temp(num,0);
        vector<double> value_temp(num,0);
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            int start_node = 0;
            for(int i = 0; i < vertexs.size(); i ++){
                int u = vertexs[i];
                if(i == 0)
                    start_node = u;
                u = node_inv_map[u];
                delta_temp[i] = delta[u];
            }
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                is_add = false;
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        value_temp[i] += delta_temp[i];
                        vector<int> neigh = csr.neighbors(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
                            if(part_v == part){
                                if(!is_add){
                                    cal_times += num * num;
                                    is_add = true;
                                }
                                delta_old[v] += delta_temp[i] * damping_factor / outDegree[u];
                                if(delta_old[v] > threshold){
                                    need_update.insert(v);
//                                    delta[v] += delta_old[v];
//                                    delta_old[v] = 0;
//                                    bm.set(node_map[v]);
                                    finished = false;
                                    active_num ++;
                                }
                            }

                        }
                        delta_temp[i] = 0;
                    }

                }
                for(auto nu : need_update){
                    if(delta_old[nu] > threshold){
                        delta_temp[node_map[nu]-start_node] += delta_old[nu];
                        delta_old[nu] = 0;
                        bm.set(node_map[nu]);
                    }
                }
            }
        }

        if(is_itr){
            for(int i = 0; i < vertexs.size(); i ++){
                int u = vertexs[i];
                u = node_inv_map[u];
                delta[u] = value_temp[i];
            }
        }

        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                value[u] += delta[u];
                vector<int> neigh = csr.neighbors(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(is_itr && part_v == part)
                        continue;
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    delta_old[v] += delta[u] * damping_factor / outDegree[u];
                    if((delta_old[v] + delta[v] > threshold && part != part_v) || (delta_old[v] > threshold && part == part_v)){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
                delta[u] = 0;
            }

        }
        for(auto nu : need_update){
//            if(delta_old[nu] + delta[nu] > threshold){
                delta[nu] += delta_old[nu];
                delta_old[nu] = 0;
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
                    priorityQueue.push(make_pair(part_nu,delta[nu]));
                    bm.set(node_map[nu]);
                }
//            }
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_value[start_id+i] = value[i] + delta_old[i];
    }

}

void pr_diag_iter_priority_plus(long long int &cal_times, int node_index, const CSR& csr, const CSR& csc, set<int> start, vector<int> node_map, vector<int> node_inv_map, vector<double> &global_delta, vector<double> &global_value, vector<int> outDegree, int start_id, int num_node ,int num){
//  int n = csr.n;
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    unordered_map<int,int> map;//first:part second:priority
    vector<double> delta(global_delta.begin()+start_id,global_delta.begin()+start_id+num_node);
    vector<double> delta_old = vector<double>(num_node,0);
    vector<double> value(global_value.begin()+start_id,global_value.begin()+start_id+num_node);
    priority_queue<pair<int,double>,vector<pair<int,double>>,cmp_pr> priorityQueue;
    int part_index = node_index / num;  //需要进行块内迭代的块的数量

    Bitmap bm(n);
    for(auto i:start){
        bm.set(node_map[i-start_id]);
    }


    Bitmap bm_part(partition);
    set<int> part_start;
//    queue<int> q;
    for(auto i:start){
        int part = get_partition(node_map[i-start_id], num);
        if(bm_part.get(part) == 0){
            bm_part.set(part);
            part_start.insert(part);
        }
//        priorityQueue.push(make_pair(part,global_delta[i]));
    }

    for(auto i:part_start){
        vector<int> vertexs = vertex_n_partion(i, n, num);
        double priority = 0;
        for(int j=0;j<vertexs.size();j++) {
            int u = vertexs[j];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                priority += delta[u];
            }
        }
        priorityQueue.push(make_pair(i,priority));
    }


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
        vector<double> delta_temp(num,0);
        vector<double> value_temp(num,0);
        if(active_num > 4){//} && part < part_index){
            is_itr = true;
            int start_node = 0;
            for(int i = 0; i < vertexs.size(); i ++){
                int u = vertexs[i];
                if(i == 0)
                    start_node = u;
                u = node_inv_map[u];
                delta_temp[i] = delta[u];
            }
            while(!finished && active_num > 0){
                active_num = 0;
                finished = true;
                is_add = false;
                for(int i=0;i<vertexs.size();i++){
                    int u = vertexs[i];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        value_temp[i] += delta_temp[i];
                        vector<int> neigh = csr.neighbors(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(node_map[v],num);
                            if(part_v == part){
                                if(!is_add){
                                    cal_times += num * num;
                                    is_add = true;
                                }
                                delta_old[v] += delta_temp[i] * damping_factor / outDegree[u];
                                if(delta_old[v] > threshold){
                                    need_update.insert(v);
//                                    delta[v] += delta_old[v];
//                                    delta_old[v] = 0;
//                                    bm.set(node_map[v]);
                                    finished = false;
                                    active_num ++;
                                }
                            }

                        }
                        delta_temp[i] = 0;
                    }

                }
                for(auto nu : need_update){
                    if(delta_old[nu] > threshold){
                        delta_temp[node_map[nu]-start_node] += delta_old[nu];
                        delta_old[nu] = 0;
                        bm.set(node_map[nu]);
                    }
                }
            }
        }

        if(is_itr){
            for(int i = 0; i < vertexs.size(); i ++){
                int u = vertexs[i];
                u = node_inv_map[u];
                delta[u] = value_temp[i];
            }
        }

        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                u = node_inv_map[u];
                value[u] += delta[u];
                vector<int> neigh = csr.neighbors(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(node_map[v],num);
                    if(is_itr && part_v == part)
                        continue;
                    if(bm_cnt[part_v] == false){
                        cal_times += num * num;
                        bm_cnt[part_v] = true;
                    }
                    delta_old[v] += delta[u] * damping_factor / outDegree[u];
                    if((delta_old[v] + delta[v] > threshold && part != part_v) || (delta_old[v] > threshold && part == part_v)){
                        need_update.insert(v);
                    }
                }
                bm.reset(node_map[u]);
                delta[u] = 0;
            }

        }
        for(auto nu : need_update){
//            if(delta_old[nu] + delta[nu] > threshold){
            delta[nu] += delta_old[nu];
            delta_old[nu] = 0;
            int part_nu = get_partition(node_map[nu], num);
            if(bm_part.get(part_nu) == 0){
//                    q.push(part_nu);
                bm_part.set(part_nu);
            }
            if(bm.get(node_map[nu]) == 0){
                bm.set(node_map[nu]);

                vector<int> vertexs = vertex_n_partion(part_nu, n, num);
                double priority = 0;
                for(int j=0;j<vertexs.size();j++) {
                    int u = vertexs[j];
                    if(bm.get(u) == 1){
                        u = node_inv_map[u];
                        priority += delta[u];
                    }
                }
                priorityQueue.push(make_pair(part_nu,priority));

//                priorityQueue.push(make_pair(part_nu,delta[nu]));
            }
//            }
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_value[start_id+i] = value[i] + delta_old[i];
    }

}

#endif