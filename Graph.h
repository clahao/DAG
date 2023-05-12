//
// Created by moyu on 2023/4/4.
//

#ifndef DAG_GRAPH_H
#define DAG_GRAPH_H
#include "global.h"

template <class E>
class Graph {
private:

public:
    bool isWeighted;
    bool islarge;

    uint num_nodes;     // 图顶点数
    uint num_edges;     // 边数
    int scc_num;       //强连通分量数量
    uint *offset;       // CSR顶点偏移量数组
    E *edgeList;        // 边数组
    uint *_offset;      // CSC顶点偏移量数组
    E *_edgeList;       // 边数组

    uint *weight;
    uint *inDegree;
    uint *outDegree;
    bool *label1;
    bool *label2;
    Adouble *delta1;
    Adouble *delta2;
    Adouble *pr_value;
    uint *value;
    uint *queue[package_num];
    uint *index;
    bool *scc_sort;

    Graph(string filename, bool _isWeighted);

    void AssignW8(uint w8, uint index);
    void _AssignW8(uint w8, uint index);

};

#endif //DAG_GRAPH_H
