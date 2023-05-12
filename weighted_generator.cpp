//
// Created by moyu on 2022/9/14.
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
    uint num_edges = stoi(line.substr(pos + 1));

    getline(ifs, line);
    ofs << line << endl;

    printf("%d %d\n", num_nodes, num_edges);
    ifs.seekg(0);
    for(int i = 0; i < 4; i ++)
        getline(ifs, line);

    for (int i = 0; i < num_edges; i++) {
        getline(ifs, line);
//        ofs << line << endl;
        ofs << line << " " << rand() % 100 + 1 << endl;
    }
//    for (int i = 0; i < num_edges; i++) {
//        ofs << rand() % 100 + 1 << endl;
//    }
    ifs.close();
    ofs.close();


}