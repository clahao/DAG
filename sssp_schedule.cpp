//
// Created by moyu on 2023/5/31.
//

//without scc
#include "scc_reorder.cpp"
#include "sssp_partition.cpp"
#include "timer.h"
#include "timer.cpp"
#include "rabbit_order/demo/reorder.cc"

using namespace std;

int main(int argc, char ** argv) {
    string filename(argv[1]);
    int num = stoi(argv[2]);
    Timer timer;
    Timer pre_timer;
    Timer process_timer;
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
    int scc_num = 1;
    pre_timer.Start();
    float pre_runtime = pre_timer.Finish();
    float processing_time = 0;
    long long int cal_times = 0;
    vector<set<int>> start_node_list(scc_num,set<int>());       //活跃顶点集合
    global_dist[source_node] = 0;
    start_node_list[0].insert(source_node);

    cout << "Processing start"<<endl;

    int node_start_index = 0;
            //构造csr和csc
            CSR csr;
            CSR csc;
            csr.n = graph.num_nodes;
            csc.n = graph.num_nodes;
            int num_edges = graph.num_edges;
            csr.row_ptr.assign(csr.n+1, 0);
            csr.col_idx.assign(num_edges, 0);
            csr.weights.assign(num_edges, 0);
            csc.row_ptr.assign(csc.n+1, 0);
            csc.col_idx.assign(num_edges, 0);
            csc.weights.assign(num_edges, 0);
            vector<int> indegree(csc.n,0);
            for(int j = 0; j < csr.n; j ++){
                int node_id = j;
                int outdegree = offset[node_id+1] - offset[node_id];
                csr.row_ptr[j+1] = csr.row_ptr[j] + outdegree;
                for(int k = 0; k < outdegree; k ++){
                    csr.col_idx[csr.row_ptr[j]+k] = neighbor[offset[node_id]+k];
                    csr.weights[csr.row_ptr[j]+k] = weight[offset[node_id]+k];
                    indegree[neighbor[offset[node_id]+k]] ++;
                }
            }

            vector<int> csc_index(csc.n,0);
            for(int j = 0; j < csc.n; j ++){
                csc.row_ptr[j+1] = csc.row_ptr[j] + indegree[j];
            }
            for(int j = 0; j < csc.n; j ++){
                int node_id = j;
                int outdegree = offset[node_id+1] - offset[node_id];
                for(int k = 0; k < outdegree; k ++){
                    int index = neighbor_scc[offset[node_id]+k];     //计算邻居顶点在csc对应的位置
                    csc.col_idx[csc.row_ptr[index] + csc_index[index]] = j;
                    csc.weights[csc.row_ptr[index] + csc_index[index]] = weight[offset[node_id]+k];
                    csc_index[index] ++;
                }
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
            long long int temp = cal_times;
            sssp_diag_iter_priority(cal_times,node_index,csr,csc,start_node_list[0],node_degree_map,node_degree_inv_map,global_dist,node_start_index,graph.num_nodes,num);
            processing_time += process_timer.Finish();
    float runtime = timer.Finish();
    cout << "Total runtime: " << runtime/1000 << " (s).\n";
    cout << "Pre-processing finished in " << pre_runtime/1000 << " (s).\n";
    cout << "Processing finished in " << processing_time/1000 << " (s).\n";
    cout << "Total processing time: " << (pre_runtime+processing_time)/1000 << " (s).\n";
    cout<<"calculation times: "<<cal_times<<endl;
//    cout<<"Proportion of calculation outside the community: "<<double(cal_times_between_communities)/cal_times_community<<endl;
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
