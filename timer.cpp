//
// Created by moyu on 2022/9/14.
//
#include "timer.h"

void Timer::Start() {
    gettimeofday(&StartingTime, NULL);
}

float Timer::Finish() {
    timeval PausingTime, ElapsedTime;
    gettimeofday(&PausingTime, NULL);
    ElapsedTime.tv_sec = PausingTime.tv_sec - StartingTime.tv_sec;
    ElapsedTime.tv_usec = PausingTime.tv_usec - StartingTime.tv_usec;
    float d = ElapsedTime.tv_sec * 1000.0 + ElapsedTime.tv_usec / 1000.0;
    return d;
}