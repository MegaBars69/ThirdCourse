#ifndef HEADER
#define HEADER
#include <string>
struct Thread{
    std::string filename;
    int thread_id;
    int total_threads;
    int count_max;
};

enum class io_status
{
    undef,
    error_open,
    error_read,
    succes
};

class Args{
    public:
        int p = 0;
        int k = 0;
        pthread_t tid = 0;
        const char* filename = nullptr;
        double res = 0;
        int count = 0;
        io_status error_type = io_status::undef;
        double error_flag = 0;
};

void* count_max_elements(void* arg);
void* thread_func(void *arg);
void reduce_sum(int p, double* a = nullptr, int n = 0);
int procces_args(Args *a);
io_status amount_of_max(FILE* file, double* result);


#endif