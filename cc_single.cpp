//
// Created by moyu on 2023/6/5.
//
#include "Graph.cpp"
#include "timer.h"
#include "timer.cpp"


void cc_sync(long long int &cal_times, const CSR& csr, set<int> start, vector<int> &global_dist, int start_id, int num_node ,int num=1){
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
//        bm_part.set(get_partition(i-start_id, num));
//        q.push(get_partition(i-start_id, num));
    }
    while (!q.empty()) {
        int part = q.front();
        q.pop();
        bm_part.reset(part);
        set<int> need_update;
        vector<int> vertexs = vertex_n_partion(part, n, num);
//        vector<bool> bm_cnt(partition,false); //模拟去0
        for(int i=0;i<vertexs.size();i++){
            int u = vertexs[i];
            if(bm.get(u) == 1){
                vector<int> neigh = csr.neighbors(u);
                for(int j=0;j<neigh.size();j++){
                    int v = neigh[j];
                    int part_v = get_partition(v,num);
//                    if(bm_cnt[part_v] == false){
//                        cal_times += num * num;
//                        bm_cnt[part_v] = true;
//                    }
                    dist_old[v] = min(dist_old[v], dist[u]);
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
    Graph<OutEdge> graph(filename, false);
    int source_node = 0;    //以重排序前的0顶点为源顶点
    if(argc == 3){
        source_node = stoi(argv[2]);
    }
    CSR csr;
    csr.n = graph.num_nodes;
    csr.row_ptr.assign(csr.n+1, 0);
    csr.col_idx.assign(graph.num_edges, 0);

    vector<int> global_dist = vector<int>(graph.num_nodes,0);
    for(int i = 0; i < graph.num_nodes; i ++)
        global_dist[i] = i;

    for(int i = 0; i <= graph.num_nodes; i ++)
        csr.row_ptr[i] = int(graph.offset[i]);
    for(int i = 0; i < graph.num_edges; i ++)
        csr.col_idx[i] = int(graph.edgeList[i].end);

    long long int cal_times = 0;
    set<int> start;
    for(int i = 0; i < graph.num_nodes; i ++)
        start.insert(i);
    Timer timer;
    timer.Start();
    cc_sync(cal_times,csr,start,global_dist,0,graph.num_nodes);
    float runtime = timer.Finish();
    cout << "Processing finished in " << runtime/1000 << " (s).\n";
    cout<<"calculation times: "<<cal_times<<endl;
    int update_num = 0;
    for(int i = 0; i < graph.num_nodes; i ++){
        if(global_dist[i] != i)
            update_num ++;
    }
    cout<<"updated vertex num: "<<update_num<<endl;
    for(int i = 0; i < min(100,(int)graph.num_nodes); i ++){
        cout<<global_dist[i]<<" ";
    }
    cout<<endl;


    return 0;

}