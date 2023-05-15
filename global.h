//
// Created by moyu on 2023/4/4.
//

#ifndef DAG_GLOBAL_H
#define DAG_GLOBAL_H
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include <random>
#include <cstdio>
#include <iomanip>
#include <locale>
#include <cstring>
#include <vector>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <sstream>
#include <atomic>
#include <queue>
#include <stack>
#include <set>
#include <unordered_map>

using namespace std;
using Adouble = atomic<double>;
const unsigned int DIST_INFINITY = std::numeric_limits<int>::max() - 1;

typedef unsigned int uint;
typedef unsigned long long ll;
#define MAX_BDFS_DEPTH 5
#define damping_factor 0.85
#define threshold 0.001
//#define heap_factor 0.005
#define package_num 10
#define package_interval 5

struct OutEdge {
    uint end;
};

struct OutEdgeWeighted {
    uint end;
    uint w8;
};

struct Edge {
    uint source;
    uint end;
};

struct EdgeWeighted {
    uint source;
    uint end;
    uint w8;
};

struct Node {
    uint ID;
    uint value;
};

class Bitmap {
public:
    Bitmap(size_t size) : data_((size + 31) / 32) {}

    inline bool get(size_t index) const {
        return (data_[index / 32] >> (index % 32)) & 1;
    }

    inline void set(size_t index) {
        data_[index / 32] |= (1u << (index % 32));
    }

    inline void reset(size_t index) {
        data_[index / 32] &= ~(1u << (index % 32));
    }

public:
    vector<unsigned int> data_;
};

// CSR structure to store the graph in compressed sparse row format
struct CSR {
    int n;                   // number of vertices
    vector<int> row_ptr;     // row pointers of the CSR matrix
    vector<int> col_idx;     // column indices of the CSR matrix
    vector<int> weights;     // edge weights of the CSR matrix

    inline vector<int> neighbors(int u) const {
        vector<int> result;
        for (int i = row_ptr[u]; i < row_ptr[u+1]; ++i) {
            result.push_back(col_idx[i]);
        }
        return result;
    }

    inline vector<int> neighbors_weights(int u) const {
        vector<int> result;
        for (int i = row_ptr[u]; i < row_ptr[u+1]; ++i) {
            result.push_back(weights[i]);
        }
        return result;
    }
};

inline vector<int> vertex_n_partion(int index, int n, int num){
    // int num = (n+partiton-1) / partiton;
    // int k = index / num; //从0开始
    vector<int> src;
    int vertex_index = index*num;
    for(int i=0;i<num&&vertex_index+i<n;i++){
        src.push_back(vertex_index + i);
    }
    return src;
}

inline int get_partition(int i, int num){
    return i/num;
}

const int INF = numeric_limits<int>::max()/2 - 1;

#endif //DAG_GLOBAL_H
