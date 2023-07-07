//
// Created by moyu on 2023/5/16.
//
#include "Graph.cpp"
#include "timer.h"
#include "timer.cpp"


void pr_sync(long long int &cal_times, const CSR& csr, set<int> start, vector<double> &global_delta, vector<double> &global_value, vector<int> outDegree, int start_id, int num_node ,int num){
    int n = num_node;
    int partition = (n+num-1) / num;
    if(start.size() == 0)
        return;

    vector<double> delta(global_delta.begin()+start_id,global_delta.begin()+start_id+num_node);
    vector<double> delta_old = vector<double>(num_node,0);
    vector<double> value(global_value.begin()+start_id,global_value.begin()+start_id+num_node);

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
                    value[u] += delta[u];
                    vector<int> neigh = csr.neighbors(u);
                    for(int j=0;j<neigh.size();j++){
                        int v = neigh[j];
                        int part_v = get_partition(v,num);
                        if(bm_cnt[part_v] == false){
                            cal_times += num * num;
                            bm_cnt[part_v] = true;
                        }
                        delta_old[v] += delta[u] * damping_factor / outDegree[u];
                        if(delta_old[v] > threshold){
                            need_update.insert(v);
                        }
                    }
                    bm.reset(u);
                    delta[u] = 0;
                }

            }
        }

        for(auto nu : need_update){
            if(delta_old[nu] > threshold){
                delta[nu] = delta_old[nu];
                delta_old[nu] = 0;
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
        global_value[start_id+i] = value[i] + delta_old[i];
    }

}

int main(int argc, char ** argv){
    string filename(argv[1]);
    int num = stoi(argv[2]);
    Graph<OutEdge> graph(filename, false);
    int source_node = 0;    //以重排序前的0顶点为源顶点
    if(argc == 4){
        source_node = stoi(argv[3]);
    }
    CSR csr;
    csr.n = graph.num_nodes;
    csr.row_ptr.assign(csr.n+1, 0);
    csr.col_idx.assign(graph.num_edges, 0);

    vector<double> global_delta = vector<double>(graph.num_nodes,0);
    vector<double> global_value = vector<double>(graph.num_nodes,0);
    vector<int> outDegree = vector<int>(graph.num_nodes,0);
    for(int i = 0; i < graph.num_nodes; i ++){
        global_delta[i] = 1 - damping_factor;
        outDegree[i] = graph.outDegree[i];
    }


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
    pr_sync(cal_times,csr,start,global_delta,global_value,outDegree,0,graph.num_nodes,num);
    float runtime = timer.Finish();
    cout << "Processing finished in " << runtime/1000 << " (s).\n";
    cout<<"calculation times: "<<cal_times<<endl;
    int update_num = 0;
    for(int i = 0; i < graph.num_nodes; i ++){
        if(global_value[i] != 1-damping_factor)
            update_num ++;
    }
    cout<<"updated vertex num: "<<update_num<<endl;
    for(int i = 0; i < min(100,(int)graph.num_nodes); i ++){
        cout<<global_value[i]<<" ";
    }
    cout<<endl;


    return 0;

}