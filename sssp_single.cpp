//
// Created by moyu on 2023/5/31.
//
#include "Graph.cpp"
#include "timer.h"
#include "timer.cpp"

void sssp_single(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node){
    int n = num_node;
    if(start.size() == 0)
        return;

    Timer timer;
    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);

    queue<int> q;
    Bitmap bm(n);
    for(auto i:start){
        bm.set(i-start_id);
        q.push(i-start_id);
    }

    while (!q.empty()) {
        int i = q.front();
        q.pop();
        set<int> need_update;
        if(bm.get(i) == 1){
            int outdegree = csr.row_ptr[i+1] - csr.row_ptr[i];
            int start_node = csr.row_ptr[i];
//            vector<int> neigh = csr.neighbors(i);
//            vector<int> neighbors_w = csr.neighbors_weights(i);
            for(int j=0;j<outdegree;j++){
                int v = csr.col_idx[start_node+j];
                int vweight = csr.weights[start_node+j];
//                timer.Start();
                dist_old[v] = min(dist_old[v], vweight + dist[i]);
//                process_time += timer.Finish();
                if(dist_old[v] < dist[v]){
                    dist[v] = dist_old[v];
                    if(bm.get(v) == 0){
                        q.push(v);
                        bm.set(v);
                    }

                }
            }
            bm.reset(i);
        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

//void sssp_async(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=1){
//    int n = num_node;
//    int partition = (n+num-1) / num;
//    if(start.size() == 0)
//        return;
//
//    vector<int> dist(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
//    vector<int> dist_old(global_dist.begin()+start_id,global_dist.begin()+start_id+num_node);
//
//    Bitmap bm_part(partition);
//    queue<int> q;
//    for(auto i:start){
//        int part = get_partition(i-start_id, num);
//        if(bm_part.get(part) == 0){
//            bm_part.set(part);
//            q.push(part);
//        }
//    }
//    while (!q.empty()) {
//        int part = q.front();
//        q.pop();
//        set<int> need_update;
////        vector<bool> bm_cnt(partition,false); //模拟去0
//
//        int u = part;
//        if(bm_part.get(u) == 1){
//            vector<int> neigh = csr.neighbors(u);
//            vector<int> neighbors_w = csr.neighbors_weights(u);
//            for(int j=0;j<neigh.size();j++){
//                int v = neigh[j];
////                int part_v = get_partition(v,num);
////                if(bm_cnt[part_v] == false){
////                    cal_times += num * num;
////                    bm_cnt[part_v] = true;
////                }
//                int vweight = neighbors_w[j];
//                dist_old[v] = min(dist_old[v], vweight + dist[u]);
//                if(dist_old[v] < dist[v]){
//                    need_update.insert(v);
//                }
////                cal_times ++;
//            }
//            bm_part.reset(u);
//        }
//
//        for(auto nu : need_update){
//            if(dist[nu] > dist_old[nu]){
//                dist[nu] = min(dist[nu], dist_old[nu]);
//                dist_old[nu] = dist[nu];
//                int part = get_partition(nu, num);
//                if(bm_part.get(part) == 0){
//                    q.push(part);
//                    bm_part.set(part);
//                }
//            }
//        }
//    }
//    for(int i = 0; i < num_node; i ++){
//        global_dist[start_id+i] = dist[i];
//    }
//
//}

void sssp_async(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=1){
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
        int part = get_partition(i-start_id, num);
        if(bm_part.get(part) == 0){
            bm_part.set(part);
            q.push(part);
        }
    }
    while (!q.empty()) {
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        unordered_set<int> need_update;
//        vector<bool> need_update(num_node,false);
        vector<int> vertexs = vertex_n_partion(part, n, num);
//        vector<bool> bm_cnt(partition,false); //模拟去0
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                int outdegree = csr.row_ptr[u+1] - csr.row_ptr[u];
                int start_node = csr.row_ptr[u];
//                vector<int> neigh = csr.neighbors(u);
//                vector<int> neighbors_w = csr.neighbors_weights(u);
                for(int j=0;j<outdegree;j++){
                    int v = csr.col_idx[start_node+j];
//                    int part_v = get_partition(v,num);
//                    if(bm_cnt[part_v] == false){
//                        cal_times += num * num;
//                        bm_cnt[part_v] = true;
//                    }
                    int vweight = csr.weights[start_node+j];
                    dist_old[v] = min(dist_old[v], vweight + dist[u]);
                    if(dist_old[v] < dist[v]){
//                        cal_times++;
                        need_update.insert(v);
//                        need_update[v] = true;
                    }
                }
                bm.reset(u);
            }
        }
//        for(int nu = 0; nu < num_node; nu ++){
//            if(need_update[nu]){
//                if(dist[nu] > dist_old[nu]){
//                    dist[nu] = min(dist[nu], dist_old[nu]);
//                    dist_old[nu] = dist[nu];
//                    int part = get_partition(nu, num);
//                    if(bm_part.get(part) == 0){
//                        q.push(part);
//                        bm_part.set(part);
//                    }
//                    bm.set(nu);
//                }
//            }
//        }
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


void sssp_sync(long long int &cal_times, float &process_time, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=1){
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
        int part = get_partition(i-start_id, num);
        if(bm_part.get(part) == 0){
            bm_part.set(part);
            q.push(part);
        }
    }
    while (!q.empty()) {
        int s = q.size();
        set<int> need_update;
//        vector<bool> need_update(num_node,false);
        while(s--){
            int part = q.front();
            q.pop();
            bm_part.reset(part);

            vector<int> vertexs = vertex_n_partion(part, n, num);
//        vector<bool> bm_cnt(partition,false); //模拟去0
            for(int i=0;i<vertexs.size();i++){
                int u = vertexs[i];
                if(bm.get(u) == 1){
                    vector<int> neigh = csr.neighbors(u);
                    vector<int> neighbors_w = csr.neighbors_weights(u);
                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int part_v = get_partition(v,num);
//                    if(bm_cnt[part_v] == false){
//                        cal_times += num * num;
//                        bm_cnt[part_v] = true;
//                    }
                        int vweight = neighbors_w[j];
                        dist_old[v] = min(dist_old[v], vweight + dist[u]);
                        if(dist_old[v] < dist[v]){
                            cal_times++;
                            need_update.insert(v);
//                        need_update[v] = true;
                        }
                    }
                    bm.reset(u);
                }
            }
        }

//        for(int nu = 0; nu < num_node; nu ++){
//            if(need_update[nu]){
//                if(dist[nu] > dist_old[nu]){
//                    dist[nu] = min(dist[nu], dist_old[nu]);
//                    dist_old[nu] = dist[nu];
//                    int part = get_partition(nu, num);
//                    if(bm_part.get(part) == 0){
//                        q.push(part);
//                        bm_part.set(part);
//                    }
//                    bm.set(nu);
//                }
//            }
//        }
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

int main(int argc, char ** argv){
    string filename(argv[1]);
    Graph<OutEdgeWeighted> graph(filename, true);
    int source_node = 0;    //以重排序前的0顶点为源顶点
    if(argc == 3){
        source_node = stoi(argv[2]);
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
    float process_time = 0;
    set<int> start;
    start.insert(source_node);
    Timer timer;
    timer.Start();
    sssp_single(cal_times,process_time,csr,start,global_dist,0,graph.num_nodes);
//    sssp_async(cal_times,process_time,csr,start,global_dist,0,graph.num_nodes);
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
