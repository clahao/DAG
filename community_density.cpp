//
// Created by moyu on 2023/6/13.
//
#include "scc_reorder.cpp"
#include "timer.h"
#include "timer.cpp"

using namespace std;

int main(int argc, char ** argv) {
    string filename(argv[1]);
    Graph<OutEdgeWeighted> graph(filename, true);

    vector<int> offset = vector<int>(graph.num_nodes + 1, 0);
    vector<int> neighbor = vector<int>(graph.num_edges, 0);
    vector<int> weight = vector<int>(graph.num_edges, 0);

    vector<int> scc = vector<int>(graph.num_nodes, 0);
    vector<int> scc_node_num = vector<int>(graph.num_nodes, 0);

    vector<int> global_dist = vector<int>(graph.num_nodes, INF);

    for (int i = 0; i <= graph.num_nodes; i++)
        offset[i] = int(graph.offset[i]);
    for (int i = 0; i < graph.num_edges; i++)
        neighbor[i] = int(graph.edgeList[i].end);
    for (int i = 0; i < graph.num_edges; i++)
        weight[i] = int(graph.edgeList[i].w8);

    CSR csr;
    csr.n = graph.num_nodes;
    csr.row_ptr.assign(csr.n+1, 0);
    csr.col_idx.assign(graph.num_edges, 0);
    csr.weights.assign(graph.num_edges, 0);
    CSR csc;
    csc.n = graph.num_nodes;
    csc.row_ptr.assign(csc.n+1, 0);
    csc.col_idx.assign(graph.num_edges, 0);
    csc.weights.assign(graph.num_edges, 0);
    vector<int> indegree(csc.n,0);
    int node_start_index = 0;
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
            int index = neighbor[offset[node_id]+k];
            csc.col_idx[csc.row_ptr[index] + csc_index[index]] = j;
            csc.weights[csc.row_ptr[index] + csc_index[index]] = weight[offset[node_id]+k];
            csc_index[index] ++;
        }
    }
    int num = 0;
    for(int i = 0; i < csr.n; i += 16){
        for(int j = i; j < min(i+16,csr.n); j ++){
            for(int k = csr.row_ptr[j]; k < csr.row_ptr[j+1]; k ++){
                if(csr.col_idx[k] > i && csr.col_idx[k] < min(i+16,csr.n)){
                    num ++;
                }
            }
        }
    }
    cout<<"edge num in community before reorder: "<<num<<endl;

    vector<vector<int>> map_degree = reorder_degree(csr,csc,16);
    vector<int> node_degree_inv_map = map_degree[0];        //重排序后的顶点id -> 原顶点id
    vector<int> node_degree_map = map_degree[1];            //原顶点id -> 重排序后的顶点id
//    for(int i = 0; i < 100; i ++)
//        cout<<node_degree_map[i]<<" ";
//    cout<<endl;
    num = 0;
    for(int i = 0; i < csr.n; i += 16){
        for(int j = i; j < min(i+16,csr.n); j ++){
            int map_node = node_degree_inv_map[j];
            for(int k = csr.row_ptr[map_node]; k < csr.row_ptr[map_node+1]; k ++){
                int node = node_degree_map[csr.col_idx[k]];
                if(node > i && node < min(i+16,csr.n)){
                    num ++;
                }
            }
        }
    }
    cout<<"edge num in community after reorder: "<<num<<endl;
    return 0;
}