#ifndef SCC_REORDER
#define SCC_REORDER
#include <iostream>
#include <vector>
#include <stack>
#include "sssp_partition.cpp"
#include "Graph.h"
#include "Graph.cpp"

using namespace std;
const unsigned int INT_INFINITY = std::numeric_limits<int>::max() - 1;

void tarjan(int u, vector<int> offset, vector<int> neighbor, stack<int> &S, stack<int> &s, vector<int> &dfn, vector<int> &low, vector<int> &scc, vector<int> parent, vector<bool> &inStack, int &dfs_cnt, int &scc_cnt) {
    dfn[u] = low[u] = ++dfs_cnt;
    inStack[u] = true;
    s.push(u);
    S.push(u);

    while (!s.empty()) {
//        cout<<s.size()<<" "<<endl;
//        for(int i = 0; i < dfn.size(); i ++)
//            cout<<"dfn "<<dfn[i]<<" low "<<low[i]<<endl;
        int v = s.top();
        bool done = true;
        for (int i = offset[v]; i < offset[v+1]; i ++) {
            int w = neighbor[i];
            if (dfn[w] == 0) {
                s.push(w);
                S.push(w);
                done = false;
                dfn[w] = low[w] = ++dfs_cnt;
                inStack[w] = true;
                parent[w] = v;
                break;
            } else if (inStack[w]) {
                low[v] = min(low[v], dfn[w]);
            }
        }

        if (done) {
//            inStack[v] = false;
//            scc[v] = scc_cnt;
            s.pop();
            low[parent[v]] = min(low[parent[v]], low[v]);
//            for (int i = offset[v]; i < offset[v+1]; i ++) {
//                int w = neighbor[i];
//                if(inStack[w]){
//                    low[w] = min(low[w], low[v]);
//                }
//            }
            if (low[v] == dfn[v]) {
                int w;
                do {
                    w = S.top();//cout<<w<<" "<<scc_cnt<<endl;
                    S.pop();
                    inStack[w] = false;
                    scc[w] = scc_cnt;
                } while (v != w);
                scc_cnt++;
            }
        }
    }
}

int find_scc(int n, vector<int> offset, vector<int> neighbor, stack<int> &S, stack<int> &s, vector<int> &dfn, vector<int> &low, vector<int> &scc, vector<int> parent, vector<bool> &inStack, int &dfsIndex, int &sccIndex) {

    for (int i = 0; i < n; i++) {
        if (dfn[i] == 0) {
//            tarjan(i,offset,neighbor,S,s,dfn,low,scc,parent,inStack,dfsIndex,sccIndex);
            dfn[i] = low[i] = ++dfsIndex;
            inStack[i] = true;
            s.push(i);
            S.push(i);
            while (!s.empty()) {
                int v = s.top();
                bool done = true;
                for (int j = offset[v]; j < offset[v+1]; j ++) {
                    int w = neighbor[j];
                    if (dfn[w] == 0) {
                        s.push(w);
                        S.push(w);
                        done = false;
                        dfn[w] = low[w] = ++dfsIndex;
                        inStack[w] = true;
                        parent[w] = v;
                        break;
                    } else if (inStack[w]) {
                        low[v] = min(low[v], dfn[w]);
                    }
                }

                if (done) {
                    s.pop();
                    low[parent[v]] = min(low[parent[v]], low[v]);
                    if (low[v] == dfn[v]) {
                        int w;
                        do {
                            w = S.top();//cout<<w<<" "<<scc_cnt<<endl;
                            S.pop();
                            inStack[w] = false;
                            scc[w] = sccIndex;
                        } while (v != w);
                        sccIndex++;
                    }
                }
            }
        }
    }
    int scc_num = 0;
    for (int i = 0; i < n; i++) {
//        cout << "Node " << i << " belongs to SCC " << scc[i] << endl;
        scc_num = max(scc_num,scc[i]+1);
    }
    return scc_num;
}


vector<int> topo_sort(int scc_num, int num_nodes, vector<int> offset, vector<int> neighbor, vector<int> scc, vector<vector<int>> &scc_set) {
    // 计算每个强连通分量的入度
    vector<int> in_degree(scc_num, 0);
    for(int i = 0; i < num_nodes; i ++){
        uint nbegin = offset[i];
        uint nend = offset[i+1];
        for(int j = nbegin; j < nend; j ++){
            uint dest = neighbor[j];
            //当边的源顶点和目的顶点处于不同的连通分量时才计算入度
            if(scc[i] != scc[dest]){
                in_degree[scc[dest]]++;
            }
        }
    }
    // 将入度为0的强连通分量加入队列
    queue<int> q;
    for (int i = 0; i < scc_num; i++) {
        if (in_degree[i] == 0) {
            q.push(i);
        }
    }
    //构建scc_num个vector，每个vector存储对应强连通分量所有顶点的id
    for(int i = 0; i < num_nodes; i ++){
        scc_set[scc[i]].push_back(i);
    }
    // 拓扑排序
    cout<<"topology order start"<<endl;
    vector<int> order;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        order.push_back(u);
        for(int i = 0; i < scc_set[u].size(); i ++){
            int node_id = scc_set[u][i];
            int nbegin = offset[node_id];
            int nend = offset[node_id+1];
            for(int j = nbegin; j < nend; j ++){
                uint dest = neighbor[j];
                //当边的源顶点和目的顶点处于不同的连通分量
                if(scc[node_id] != scc[dest]){
                    in_degree[scc[dest]]--;
                    if (in_degree[scc[dest]] == 0) {
                        q.push(scc[dest]);
                    }
                }
            }
        }
//        for(int i = 0; i < num_nodes; i ++){
//            //当前顶点所在的强连通分量从队列中取出
//            if(scc[i] == u){
//                uint nbegin = offset[i];
//                uint nend = offset[i+1];
//                for(int j = nbegin; j < nend; j ++){
//                    uint dest = neighbor[j];
//                    //当边的源顶点和目的顶点处于不同的连通分量
//                    if(scc[i] != scc[dest]){
//                        in_degree[scc[dest]]--;
//                        if (in_degree[scc[dest]] == 0) {
//                            q.push(scc[dest]);
//                        }
//                    }
//                }
//            }
//
//        }
    }
    cout<<"topology order finfshed"<<endl;
    return order;
}


vector<int> scc_reorder(int &source_node, int &max_scc_id, int num_nodes, int num_edges, int scc_num, vector<int> offset, vector<int> neighbor, vector<int> weight,
                 vector<int> &scc, vector<int> &scc_node_num, vector<int> &offset_scc, vector<int> &neighbor_scc, vector<int> &weight_scc, vector<int> &split_scc){

    vector<vector<int>> scc_set(scc_num,vector<int>());

    //拓扑排序
    vector<int> topology_order = topo_sort(scc_num,num_nodes,offset,neighbor,scc,scc_set);
//    for(int i = 0; i < topology_order.size(); i++){
//        cout<<topology_order[i]<<" ";
//    }


    //顶点重排序后的映射
    vector<int> map = vector<int>(num_nodes,0);
    vector<int> inv_map = vector<int>(num_nodes,0);
    int id_index = 0;
    for(int i = 0; i < topology_order.size(); i ++){
        int scc_id = topology_order[i];
        for(int j = 0; j < scc_set[scc_id].size(); j ++){
            int node_id = scc_set[scc_id][j];
            map[id_index] = node_id;        //重排序后的顶点id -> 原顶点id
            inv_map[node_id] = id_index;    //原顶点id -> 重排序后的顶点id
            id_index ++;
        }
    }
    for(int i = 0; i < scc_num; i ++){
        scc_node_num[i] = scc_set[topology_order[i]].size();
    }
    cout<<"vertex map finfshed"<<endl;
//    cout<<"source node id: "<<map[source_node]<<endl;
    cout<<"source node id: "<<source_node<<endl;

    //重新构造csr
    id_index = 0;
    offset_scc[0] = 0;
    for(int i = 0; i < topology_order.size(); i ++){
        int scc_id = topology_order[i];
        for(int j = 0; j < scc_set[scc_id].size(); j ++){
            int node_id = scc_set[scc_id][j];
            int node_outdegree = offset[node_id+1] - offset[node_id];
            offset_scc[id_index+1] = offset_scc[id_index] + node_outdegree;
            int nbegin = offset[node_id];
            int node_index = 0;
            for(int k = 0; k < node_outdegree; k ++){
                if(scc[neighbor[nbegin+k]] == scc_id){
                    neighbor_scc[node_index+offset_scc[id_index]] = inv_map[neighbor[nbegin+k]];
                    weight_scc[node_index+offset_scc[id_index]] = weight[nbegin+k];
                    node_index ++;
                }
            }
            split_scc[id_index] = node_index + offset_scc[id_index];
            for(int k = 0; k < node_outdegree; k ++){
                if(scc[neighbor[nbegin+k]] != scc_id){
                    neighbor_scc[node_index+offset_scc[id_index]] = inv_map[neighbor[nbegin+k]];
                    weight_scc[node_index+offset_scc[id_index]] = weight[nbegin+k];
                    node_index ++;
                }
            }
            id_index ++;
        }
    }
    //重排连通分量的id
    int max_scc_size = 0;
    max_scc_id = 0;
    id_index = 0;
    for(int i = 0; i < topology_order.size(); i ++){
        int scc_id = topology_order[i];
        if(max_scc_size < scc_set[scc_id].size()){
            max_scc_size = scc_set[scc_id].size();
            max_scc_id = i;
        }
        for(int j = 0; j < scc_set[scc_id].size(); j ++){
            scc[id_index] = i;
            id_index ++;
        }
    }
    cout<<"vertex reorder finfshed"<<endl;
    cout<<"largest scc id: "<<max_scc_id<<endl;
    cout<<"vertex ratio of the largest scc: "<<double(max_scc_size)/num_nodes<<endl;
    return inv_map;
}


vector<int> top_scc(int &source_node, int &max_scc_id, int num_nodes, int num_edges, vector<int> offset, vector<int> neighbor, vector<int> weight,
             int &scc_num, vector<int> &scc ,vector<int> &scc_node_num, vector<int> &offset_scc, vector<int> &neighbor_scc, vector<int> &weight_scc, vector<int> &split_scc){
    int dfsIndex = 0;
    int sccIndex = 0;
    stack<int> s;
    stack<int> S;
    vector<int> dfn = vector<int>(num_nodes,0);
    vector<int> low = vector<int>(num_nodes,INT_INFINITY);
    vector<int> parent = vector<int>(num_nodes,0);
    vector<bool> inStack = vector<bool>(num_nodes,0);

    cout<<"findSCC start"<<endl;
    scc_num = find_scc(num_nodes,offset,neighbor,S,s,dfn,low,scc,parent,inStack,dfsIndex,sccIndex);
    cout<<"scc num: "<<scc_num<<endl;
    cout<<"findSCC finished"<<endl;

    return scc_reorder(source_node,max_scc_id,num_nodes,num_edges,scc_num,offset,neighbor,weight,scc,scc_node_num,offset_scc,neighbor_scc,weight_scc,split_scc);

}

bool cmp(pair<int,int> &a, pair<int,int> &b){
    return a.second < b.second;
}

vector<vector<int>> reorder_degree(CSR csr, CSR csc, int num){
    int num_node = csr.n;
    vector<int> node_map(num_node,0);       //重排序后的顶点id -> 原顶点id
    vector<int> node_map_inv(num_node,0);   //原顶点id -> 重排序后的顶点id
    vector<bool> is_visited(num_node,false);
    vector<int> degree(num_node,0);
    vector<pair<int,int>> index(num_node,pair<int,int>(0,0));//first:id second:degree
    for(int i = 0; i < num_node; i ++){
        index[0].first = i;
    }
    for(int i = 0; i < num_node; i ++){
        for(int j = csr.row_ptr[i]; j < csr.row_ptr[i+1]; j ++){
            index[csr.col_idx[j]].second ++;
            degree[csr.col_idx[j]] ++;
        }
    }
    for(int i = 0; i < num_node; i ++){
        for(int j = csc.row_ptr[i]; j < csc.row_ptr[i+1]; j ++){
            index[csc.col_idx[j]].second ++;
            degree[csc.col_idx[j]] ++;
        }
    }
    sort(index.begin(),index.end(),cmp);
    int node_index = 0;
    for(int i = 0; i < num_node; i ++){
        int node_id = i;//index[i].first;
        if(!is_visited[node_id]){   //从当前未访问顶点中度数最高的顶点开始
            is_visited[node_id] = true;
            int node_in_block = 1;  //已经找到的顶点数(总共16个)
            queue<int> node_q;      //已经找到的顶点
            node_q.push(node_id);
            unordered_map<int,int> map;
            bool finished = false;
            while(node_in_block < num && !finished){
                finished = true;
                int max_connected = 0;
                int prior_node = node_id;
                //出边邻居
                for(int j = csr.row_ptr[node_id]; j < csr.row_ptr[node_id+1]; j ++){
                    int neighbor = csr.col_idx[j];
                    if(is_visited[neighbor])
                        continue;
                    if(map.count(neighbor) == 0){
                        map[neighbor] = 1;
                    }
                    else{
                        map[neighbor] ++;
                    }
                }
                //入边邻居
                for(int j = csc.row_ptr[node_id]; j < csc.row_ptr[node_id+1]; j ++){
                    int neighbor = csc.col_idx[j];
                    if(is_visited[neighbor])
                        continue;
                    if(map.count(neighbor) == 0){
                        map[neighbor] = 1;
                    }
                    else{
                        map[neighbor] ++;
                    }
                }
                for(auto it = map.begin(); it != map.end(); it ++){
                    if((*it).second > 0 && ((*it).second > max_connected )){ //|| ((*it).second == max_connected && degree[(*it).first] > degree[prior_node]))){
                        max_connected = (*it).second;
                        prior_node = (*it).first;
                        finished = false;
                    }
                }
                if(!finished){
                    node_in_block ++;
                    is_visited[prior_node] = true;
                    node_q.push(prior_node);
                    map[prior_node] = 0;
                    node_id = prior_node;
                }
            }
            if(node_q.size() == num){
                int q_size = node_q.size();
                for(int j = 0; j < q_size; j ++){
                    int node = node_q.front();
                    node_q.pop();
                    node_map[node_index] = node;
                    node_map_inv[node] = node_index ++;
                }
            }
            else{
                int q_size = node_q.size();
                for(int j = 0; j < q_size; j ++){
                    int node = node_q.front();
                    node_q.pop();
                    is_visited[node] = false;
                }
            }
        }
    }
    if(csr.n > 1000)
    cout<<"num of node in 16-size-block:"<<csr.n<<" "<<node_index<<endl;
    vector<int> index_num;
    index_num.push_back(node_index);
    if(node_index < num_node){
        for(int i = 0; i < num_node; i ++){
//            int node_id = index[i].first;
            if(!is_visited[i]){
                node_map[node_index] = i;
                node_map_inv[i] = node_index ++;
            }
        }
    }
    vector<vector<int>> result;
    result.push_back(node_map);
    result.push_back(node_map_inv);
    result.push_back(index_num);
    return result;
}

//int main(int argc, char ** argv) {
//    string filename(argv[1]);
//    string outputname(argv[2]);
//    Graph<OutEdgeWeighted> graph(filename, true);
//
//    vector<int> offset = vector<int>(graph.num_nodes+1,0);
//    vector<int> neighbor = vector<int>(graph.num_edges,0);
//    vector<int> weight = vector<int>(graph.num_edges,0);
//
//    vector<int> scc = vector<int>(graph.num_nodes,0);
//    vector<int> scc_node_num = vector<int>(graph.num_nodes,0);
//    vector<int> offset_scc = vector<int>(graph.num_nodes+1,0);
//    vector<int> neighbor_scc = vector<int>(graph.num_edges,0);
//    vector<int> weight_scc = vector<int>(graph.num_edges,0);
//    vector<int> split_scc = vector<int>(graph.num_nodes,0);
//
//    for(int i = 0; i <= graph.num_nodes; i ++)
//        offset[i] = int(graph.offset[i]);
//    for(int i = 0; i < graph.num_edges; i ++)
//        neighbor[i] = int(graph.edgeList[i].end);
//    for(int i = 0; i < graph.num_edges; i ++)
//        weight[i] = int(graph.edgeList[i].w8);
//
////    cout<<"findSCC start"<<endl;
//////    int scc_num = findSCC(n,offset,neighbor,s,dfn,low,scc,inStack,dfsIndex,sccIndex);
////    int scc_num = find_scc(n,offset,neighbor,S,s,dfn,low,scc,parent,inStack,dfsIndex,sccIndex);
////    cout<<"scc num: "<<scc_num<<endl;
////    cout<<"findSCC finished"<<endl;
////    for (int i = 0; i < n; i++) {
////        cout << "Node " << i << " belongs to SCC " << scc[i] << endl;
////        graph.scc_num = max(graph.scc_num,scc[i]+1);
////    }
//
////    scc_reorder(graph.num_nodes,graph.num_edges,scc_num,offset,neighbor,weight,scc,scc_node_num,offset_scc,neighbor_scc,weight_scc,split_scc);
//    int scc_num;
//    top_scc(graph.num_nodes,graph.num_edges,offset,neighbor,weight,scc_num,scc,scc_node_num,offset_scc,neighbor_scc,weight_scc,split_scc);
//    ofstream ofs;
//    ofs.open(outputname);
//    if(ofs.is_open()){
//        ofs<<graph.num_nodes<<" "<<graph.num_edges<<endl;
//        ofs<<scc_num<<endl;
//        for(int i = 0; i < scc_num; i ++){
//            ofs<<scc_node_num[i]<<" ";
//        }
//        ofs<<endl;
//        for(int i = 0; i <= graph.num_nodes; i ++){
//            ofs<<offset_scc[i]<<" ";
//        }
//        ofs<<endl;
//        for(int i = 0; i < graph.num_nodes; i ++){
//            ofs<<split_scc[i]<<" ";
//        }
//        ofs<<endl;
//        for(int i = 0; i < graph.num_edges; i ++){
//            ofs<<neighbor_scc[i]<<" ";
//        }
//        ofs<<endl;
//        for(int i = 0; i < graph.num_edges; i ++){
//            ofs<<weight_scc[i]<<" ";
//        }
//        ofs<<endl;
//    }
//    ofs.close();
//    return 0;
//}
#endif
