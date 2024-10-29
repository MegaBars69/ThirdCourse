#ifndef HEADER
#define HEADER
#include <string>
void* count_max_elements(void* arg);
struct Thread{
    std::string filename;
    int thread_id;
    int total_threads;
    int count_max;
};
#endif