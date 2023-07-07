//
// Created by moyu on 2023/6/26.
//
#include "Graph.h"
#include "Graph.cpp"

int main(int argc, char **argv) {
    string input_file(argv[1]);
    string output_file(argv[2]);
    string output_file_w(argv[3]);
    Graph<OutEdgeWeighted> graph(input_file, true);

    ofstream ofs;
    ofs.open(output_file);
    srand((int)time(0));
    string line;

    cout<< output_file << " " << output_file_w << endl;
    cout<< graph.num_nodes << " " << graph.num_edges <<endl;
    ofs << "AdjacencyGraph" << endl;
    ofs << graph.num_nodes <<endl;
    ofs << graph.num_edges << endl;

    for(int i = 0 ; i < graph.num_nodes; i ++){
        ofs << graph.offset[i] << endl;
    }
    for(int i = 0 ; i < graph.num_edges; i ++){
        ofs << graph.edgeList[i].end << endl;
    }
    ofs.close();

    ofs.open(output_file_w);
    ofs << "WeightedAdjacencyGraph" << endl;
    ofs << graph.num_nodes <<endl;
    ofs << graph.num_edges << endl;
    for(int i = 0 ; i < graph.num_nodes; i ++){
        ofs << graph.offset[i] << endl;
    }
    for(int i = 0 ; i < graph.num_edges; i ++){
        ofs << graph.edgeList[i].end << endl;
    }
    for(int i = 0 ; i < graph.num_edges; i ++){
        ofs << graph.edgeList[i].w8 << endl;
    }
    ofs.close();
}