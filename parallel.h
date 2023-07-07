//
// Created by moyu on 2023/6/1.
//

#ifndef DAG_PARALLEL_H
#define DAG_PARALLEL_H
#include <omp.h>
#include <iostream>
//#define parallel_for _Pragma("omp parallel for") for

template <class ET> inline bool CAS(ET *ptr, ET oldv, ET newv){
    if(sizeof(ET) == 1){
        return __sync_bool_compare_and_swap((bool *)ptr, *((bool *)&oldv),
                                            *((bool *)&newv));
    }else if(sizeof(ET) == 4){
        return __sync_bool_compare_and_swap((int *)ptr, *((int *)&oldv),
                                            *((int *)&newv));
    }else if(sizeof(ET) == 8){
        return __sync_bool_compare_and_swap((long *)ptr, *((long *)&oldv),
                                            *((long *)&newv));
    }else{
        std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
        abort();
    }
}
#endif //DAG_PARALLEL_H
