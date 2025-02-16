#ifndef HEADER
#define HEADER
#include <string>
template <typename T>
void reduce_sum(int p,T* a = nullptr, int n = 0)
{
    static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T* r = nullptr;
    int i;
    if (p <= 1)
    {
        return;
    }

    pthread_mutex_lock(&my_mutex);

    if (r==nullptr)
    {
        r = a;
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            r[i] += a[i];
        }
    }

    t_in++;

    if (t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while (t_in < p)
        {
            pthread_cond_wait(&c_in, &my_mutex);
        }
    }

    if(r != a)
    {
        for (i = 0; i < n; i++)
        {
            a[i] = r[i];
        }
    }
    t_out++;
    if (t_out >= p)
    {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while (t_out < p)
        {
            pthread_cond_wait(&c_out, &my_mutex);
        }
        
    }
    pthread_mutex_unlock(&my_mutex);
}

enum class io_status
{
    undef,
    error_open,
    error_read,
    error_av_doesnt_exist,
    succes
};

class Args{
    public:
        int p = 0;
        int k = 0;
        pthread_t tid = 0;
        const char* filename = nullptr;
        double sum = 0;
        int count = 0;
        io_status error_type = io_status::undef;
        double error_flag = 0;
        double avarage = 0;

        double first = -std::numeric_limits<double>::infinity();
        double last = -std::numeric_limits<double>::infinity();
        bool last_can_be_a_minimum = false;
        bool first_can_be_a_minimum = false;
        int amount_of_min = 0;
};

void* count_max_elements(void* arg);
void* thread_func(void *arg);
//void reduce_sum(int p, double* a = nullptr, int n = 0);
int procces_args(Args *a);
int procces_results(Args* Result, int p);
io_status amount_of_elements(FILE* file, double* Sum, int* Count);
io_status amount_of_local_min(FILE* file, Args* Result, double avarage);
io_status avarage(double* result, double* sum, int* amount);

#endif