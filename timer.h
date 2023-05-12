//
// Created by moyu on 2022/9/14.
//

#ifndef GRAPH_ALGORITHM_TIMER_H
#define GRAPH_ALGORITHM_TIMER_H
#include "global.h"
#include <stdlib.h>
#include <sys/time.h>

class Timer {
private:
    timeval StartingTime;
public:
    void Start();
    float Finish();
};
#endif //GRAPH_ALGORITHM_TIMER_H
