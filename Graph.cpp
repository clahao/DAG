//
// Created by moyu on 2023/4/4.
//
#ifndef DAG_GRAPH_CPP
#define DAG_GRAPH_CPP
#include "Graph.h"
#include <string>

template <>
inline void Graph<OutEdgeWeighted>::AssignW8(uint w8, uint index) {
    edgeList[index].w8 = w8;
}

template <>
inline void Graph<OutEdge>::AssignW8(uint w8, uint index) {
    edgeList[index].end = edgeList[index].end;
}

template <class E>
Graph<E>::Graph(string filename, bool _isWeighted) {
    isWeighted = _isWeighted;

    srand((int)time(0));
    ifstream ifs;
    ifs.open(filename);
    string line;
    getline(ifs, line);
    getline(ifs, line);
    getline(ifs, line);
    line = line.substr(9);
    int pos = line.find(' ');
    num_nodes = stoi(line.substr(0, pos));
    pos = line.find_last_of(' ');
    num_edges = stoi(line.substr(pos + 1));
    getline(ifs, line);

    char c = '\t';
//    char c = ' ';

    for (int i = 0; i < num_edges; i++) {
        getline(ifs, line);
        pos = line.find(c);
        uint source = stoi(line.substr(0, pos));
        uint dest = stoi(line.substr(pos + 1));
        num_nodes = max(num_nodes, source + 1);
        num_nodes = max(num_nodes, dest + 1);
    }
    printf("num nodes:%d num edges:%d\n", num_nodes, num_edges);
    ifs.seekg(0);
    for(int i = 0; i < 4; i ++)
        getline(ifs, line);

    offset = (uint *)malloc((num_nodes+1) * sizeof(uint));
    edgeList = (E *)malloc(num_edges * sizeof(E));
    _offset = (uint *)malloc((num_nodes+1) * sizeof(uint));
    _edgeList = (E *)malloc(num_edges * sizeof(E));
    outDegree = (uint *)malloc(num_nodes * sizeof(uint));
    label1 = (bool *)malloc(num_nodes * sizeof(bool));
    label2 = (bool *)malloc(num_nodes * sizeof(bool));
    value = (uint *)malloc(num_nodes * sizeof(uint));
    scc_sort = (bool *)malloc(num_nodes * sizeof(bool));


    scc_num = 0;
    for (int i = 0; i < num_nodes; i++){
        outDegree[i] = 0;
        scc_sort[i] = false;
    }

    for (int i = 0; i < num_edges; i++) {
        getline(ifs, line);
        pos = line.find(c);
        uint source = stoi(line.substr(0, pos));
        uint dest = stoi(line.substr(pos + 1));
        outDegree[source] ++;
        edgeList[i].end = dest;
    }
    if(isWeighted)
        for (int i = 0; i < num_edges; i++) {
            getline(ifs, line);
//            weight[i] = stoi(line);
            AssignW8(stoi(line), i);
        }
    uint cnt = 0;
    for (int i = 0; i <= num_nodes; i++) {
        offset[i] = cnt;
        if(i<num_nodes)
            cnt += outDegree[i];
    }

    ifs.close();
    printf("read finished\n");
}
#endif