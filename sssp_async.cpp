//
// Created by moyu on 2023/5/10.
//
#include "Graph.cpp"
#include "timer.h"
#include "timer.cpp"

void sssp_async(long long int &cal_times, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num){
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
        set<int> need_update;
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

void sssp_async_stack(long long int &cal_times, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num){
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
    stack<int> s;
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
        set<int> need_update;
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
        for(auto nu : need_update){
            if(dist[nu] > dist_old[nu]){
                dist[nu] = min(dist[nu], dist_old[nu]);
                dist_old[nu] = dist[nu];
                int part = get_partition(nu, num);
                if(bm_part.get(part) == 0){
                    s.push(part);
                    bm_part.set(part);
                }
                bm.set(nu);
            }
        }
        while(!s.empty()){
            int k = s.top();
            s.pop();
            q.push(k);
        }
//        if(num_node > 1){
//            cal_times += num * num_node;
//        }
    }
    for(int i = 0; i < num_node; i ++){
        global_dist[start_id+i] = dist[i];
    }

}

void sssp_async_diag(long long int &cal_times, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num){
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
//                        u = node_inv_map[u];
                        vector<int> neigh = csr.neighbors(u);
                        vector<int> neighbors_w = csr.neighbors_weights(u);
                        for(int j=0;j<neigh.size();j++){
                            int v = neigh[j];
                            int part_v = get_partition(v,num);
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
                        bm.set(nu);
                    }
                }
            }
        }
        need_update.clear();

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
    Timer timer;
    timer.Start();
//    sssp_async(cal_times,csr,start,global_dist,0,graph.num_nodes,num);
    sssp_async_stack(cal_times,csr,start,global_dist,0,graph.num_nodes,num);
//    sssp_async_diag(cal_times,csr,start,global_dist,0,graph.num_nodes,num);
    float runtime = timer.Finish();
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