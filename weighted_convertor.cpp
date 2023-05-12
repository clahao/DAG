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
    cout<<num_nodes<<" "<<num_edges<<endl;

    getline(ifs, line);
    ofs << line << endl;

    vector<int> weight(num_edges,0);

    for (int i = 0; i < num_edges; i++) {
        getline(ifs, line);
        pos = line.find_last_of(' ');
        int w = stoi(line.substr(pos + 1));
        weight[i] = w;
        line = line.substr(0,pos);
        ofs << line << endl;

//        pos = line.find(' ');
//        ofs << line.substr(0,pos) << '\t' << line.substr(pos+1) << endl;
    }

    for(int i = 0; i < num_edges; i ++)
        ofs << weight[i] << endl;

    ifs.close();
    ofs.close();


}