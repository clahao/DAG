//
// Created by moyu on 2023/5/16.
//
#include "scc_reorder.cpp"
#include "pr_partition.cpp"
#include "timer.h"
#include "timer.cpp"

using namespace std;

int main(int argc, char ** argv) {
    string filename(argv[1]);
    int num = stoi(argv[2]);
    Timer timer;
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

    vector<double> global_delta = vector<double>(graph.num_nodes,0);
    vector<double> global_value = vector<double>(graph.num_nodes,0);

    for(int i = 0; i < graph.num_nodes; i ++)
        global_delta[i] = 1-damping_factor;

    for(int i = 0; i <= graph.num_nodes; i ++)
        offset[i] = int(graph.offset[i]);
    for(int i = 0; i < graph.num_edges; i ++)
        neighbor[i] = int(graph.edgeList[i].end);
    for(int i = 0; i < graph.num_edges; i ++)
        weight[i] = int(graph.edgeList[i].w8);

    int scc_num;
    int max_scc_id = 0;
    vector<int> map = top_scc(source_node,max_scc_id,graph.num_nodes,graph.num_edges,offset,neighbor,weight,scc_num,scc,scc_node_num,offset_scc,neighbor_scc,weight_scc,split_scc);

    int node_start_index = 0;     //当前连通分量第一个顶点的id
    long long int cal_times = 0;
    long long int cal_times_community = 0;
    long long int cal_times_between_communities = 0;
    vector<set<int>> start_node_list(scc_num,set<int>());       //活跃顶点集合
    vector<set<int>> neighbor_node_list(scc_num,set<int>());    //连通分量的邻居顶点集合
//    start_node_list[scc[map[source_node]]].insert(map[source_node]);

    timer.Start();
    cout << "Processing start"<<endl;
    for(int i = 0; i < scc_num; i ++){
        for(int j = 0; j < scc_node_num[i]; j ++){
            start_node_list[i].insert(node_start_index+j);
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
            vector<int> ouDegree(csr.n,0);
            vector<int> indegree(csc.n,0);
            for(int j = 0; j < csr.n; j ++){
                int node_id = node_start_index + j;
                int outdegree = split_scc[node_id] - offset_scc[node_id];
                ouDegree[j] = offset_scc[node_id+1] - offset_scc[node_id];
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
            vector<vector<int>> map_degree = reorder_degree(csr,csc,num);
            vector<int> node_degree_inv_map = map_degree[0];        //重排序后的顶点id -> 原顶点id
            vector<int> node_degree_map = map_degree[1];            //原顶点id -> 重排序后的顶点id
            int node_index = map_degree[2][0];
            //当前连通分量内部的状态传递
            long long int temp = cal_times;
//            pr_diag_iter(cal_times,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_delta,global_value,ouDegree,node_start_index,scc_node_num[i],num);
//            pr_diag_iter_priority_plus(cal_times,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_delta,global_value,ouDegree,indegree,node_start_index,scc_node_num[i],num);
            pr_diag_iter_priority(cal_times,node_index,csr,csc,start_node_list[i],node_degree_map,node_degree_inv_map,global_delta,global_value,ouDegree,node_start_index,scc_node_num[i],num);
        }
        else if(scc_node_num[i] == 1){
            global_value[node_start_index] += global_delta[node_start_index];
        }

//        if(i != max_scc_id)
//            cal_times = temp;
//        cout<<i<<" "<<start_node_list[i].size()<<endl;
        //连通分量之间的状态传递
//        long long int cal_size = 0;   //连通分量间计算的矩阵维度
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
            if(global_value[node_id] > threshold){
                int outdegree = offset_scc[node_id+1] - split_scc[node_id];
                for(int k = 0; k < outdegree; k ++){
                    int neighbor_id = neighbor_scc[split_scc[node_id]+k];
                    if(visit[neighbor_id/num] == false){
                        visit[neighbor_id/num] = true;
                        part_num ++;
                    }
                    global_delta[neighbor_id] += global_value[node_id] * damping_factor / (offset_scc[node_id+1] - offset_scc[node_id]);

                }
            }

        }
//        cal_times += scc_node_num[i] * cal_size;

        cal_times += part_num * num * num;
        node_start_index += scc_node_num[i];
    }
    float runtime = timer.Finish();
    cout << "Processing finished in " << runtime/1000 << " (s).\n";
    cout<<"calculation times: "<<cal_times<<endl;
    cout<<"calculation times col: "<<cal_times_community<<endl;
//    cout<<"Proportion of calculation outside the community: "<<double(cal_times_between_communities)/cal_times_community<<endl;
    int update_num = 0;
    for(int i = 0; i < graph.num_nodes; i ++){
        if(global_value[map[i]] != 1 - damping_factor)
            update_num ++;
    }
    cout<<"updated vertex num: "<<update_num<<endl;
    for(int i = 0; i < min(100,(int)graph.num_nodes); i ++){
        cout<<global_value[map[i]]<<" ";
    }
    cout<<endl;


    return 0;
}