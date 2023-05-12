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

#endif //DAG_GLOBAL_H
