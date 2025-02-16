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

void PrintMatrix(double* matrix, int n1, int n2);

enum class io_status
{
    undef,
    error_open,
    error_read,
    error_few_el,
    error_many_el,
    succes
};

class Args{
    public:
        int p = 0;
        int k = 0;
        int n1 = 0;
        int n2 = 0;
        int m1 = 0;
        int m2 = 0;        
        pthread_t tid = 0;
        double cpu_time_of_thread = 0;
        double cpu_time_of_all_threads = 0;

        double* A = nullptr;
        double* next_line = nullptr;
        double* prev_line = nullptr;
        io_status reading_file = io_status::undef;
        int m = 0;
        int l = 0;

        double sum = 0;
        int amount_of_el = 0;
        double avarage = 0;

        int amount_of_changed = 0;
        
        int res = 0;

        void PrintAll(bool more_info = false) const 
        {
            printf("\n");
            std::cout << "Args Data:" << std::endl;
            std::cout << "p: " << p << std::endl;
            std::cout << "k: " << k << std::endl;
            std::cout << "tid: " << tid << std::endl;
            std::cout << "sum: " << sum<< std::endl;
            std::cout << "el in sum: " << amount_of_el<< std::endl;
            std::cout<<"new element: "<<avarage<<std::endl;
            
            if(more_info)
            {
                std::cout << "m: " << m << std::endl;
                std::cout << "n1: " << n1 << std::endl;
                std::cout << "n2: " << n2 << std::endl;
                PrintMatrix(prev_line, 1, n2);
                PrintMatrix(A, m1,m2);
                PrintMatrix(next_line, 1, n2);
            }
            printf("cpu time of thread: %8.2e \n", cpu_time_of_thread);
            printf("cpu time of all threads: %8.2e \n", cpu_time_of_all_threads);
            std::cout << "amount of changed: " << amount_of_changed << std::endl;

            printf("\n");
        }
        
};
void UpdateElements(Args* a);
void get_avarage(Args *a);
void* thread_func(void *arg);



#endif