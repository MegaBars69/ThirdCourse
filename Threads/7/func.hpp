#ifndef HEADER
#define HEADER
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

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


class Args{
    public:
        int p = 0;
        int k = 0;
        int n = 0;        
        int chunk = 0;

        pthread_t tid = 0;

        double cpu_time_of_thread = 0;
        double cpu_time_of_all_threads = 0;  

        unsigned long long int L = 5;
        unsigned long long int answer = 0; 

        int founded_pairs = 0;

        int founded_pairs_in_the_last_interval = 0;

        void  PrintAll(bool only_cpu_info)
        {
            std::cout<<"----------------------------------------"<<std::endl;
            std::cout<<"k = "<<k<<std::endl;
            if (only_cpu_info)
            {
                std::cout<<"p = "<<p<<std::endl;
                std::cout<<"n = "<<n<<std::endl;
                std::cout<<"chunk = "<<chunk<<std::endl;
                std::cout<<"found pairs = "<<founded_pairs<<std::endl;
                std::cout<<"found pairs in the last interval = "<<founded_pairs_in_the_last_interval<<std::endl;
                std::cout<<"L = "<<L<<std::endl;
            
            }
            std::cout<<"Sum of CPU time of all threads = "<<cpu_time_of_all_threads<<std::endl;
            std::cout<<"CPU time of thread = "<<cpu_time_of_thread<<std::endl;
            std::cout<<"----------------------------------------"<<std::endl;

        }

};

bool NumberAndNextArePrime(unsigned long long int num);

void FindInterval(Args* a);

void ProccesResult(Args* a);

void* thread_func(void *arg);



#endif