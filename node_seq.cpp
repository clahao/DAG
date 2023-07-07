//
// Created by moyu on 2023/5/23.
//
#include "global.h"

int main(int argc, char **argv) {
    string input_file(argv[1]);
    string output_file(argv[2]);

    ifstream ifs;
    ifs.open(input_file);
    ofstream ofs;
    ofs.open(output_file);

    srand((int)time(0));
    string line;
    getline(ifs, line);
    ofs << line << endl;
    getline(ifs, line);
    ofs << line << endl;
    getline(ifs, line);
    ofs << line << endl;

    line = line.substr(9);
    int pos = line.find(' ');
    int num_nodes = stoi(line.substr(0, pos));
    pos = line.find_last_of(' ');
    int num_edges = stoi(line.substr(pos + 1));

    getline(ifs, line);
    ofs << line << endl;

    char c = '\t';
    for (int i = 0; i < num_edges; i++) {
        getline(ifs, line);
        pos = line.find(c);
        int source = stoi(line.substr(0, pos));
        int dest = stoi(line.substr(pos + 1));
        num_nodes = max(num_nodes, source + 1);
        num_nodes = max(num_nodes, dest + 1);
    }
    printf("num nodes:%d num edges:%d\n", num_nodes, num_edges);

    ifs.seekg(0);
    for(int i = 0; i < 4; i ++)
        getline(ifs, line);

    vector<int> outDegree(num_nodes,0);
    vector<int> offset(num_nodes+1,0);
    vector<int> neighbor(num_edges,0);
    vector<int> weight(num_edges,0);
    vector<int> map(num_nodes,-1);
    int node_index = 0;
    int pre_source = -1;
    for (int i = 0; i < num_edges; i++) {
        getline(ifs, line);
        pos = line.find(c);
        uint source = stoi(line.substr(0, pos));
        uint dest = stoi(line.substr(pos + 1));
        if(pre_source == -1){
            pre_source = source;
            map[source] = node_index;
        }
        else if(source != pre_source){
            pre_source = source;
            node_index ++;
            map[source] = node_index;
        }
        weight[i] = rand() % 100 + 1;
        outDegree[node_index] ++;
        neighbor[i] = dest;

    }

    uint cnt = 0;
    for (int i = 0; i <= num_nodes; i++) {
        offset[i] = cnt;
        if(i<num_nodes)
            cnt += outDegree[i];
    }

    for(int i = 0; i < num_nodes; i ++){
        if(map[i] == -1)
            continue;
        for(int j = offset[map[i]]; j < offset[map[i] + 1]; j ++){
            ofs << i << '\t' << neighbor[j] << endl;
        }
    }

    for(int i = 0; i < num_edges; i ++)
        ofs << weight[i] << endl;

    ifs.close();
    ofs.close();


}