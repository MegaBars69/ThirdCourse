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
        int n = 0;
    
        pthread_t tid = 0;
        double* array = nullptr;
        int m = 0;
        int l = 0;

        int amount_of_changed = 0;
        double pan = 0;  //prev array's elements
        double na0 = 0, na1 = 0;  //next array's elemnets
        
        void PrintAll(void) const 
        {
            printf("\n");
            std::cout << "Args Data:" << std::endl;
            std::cout << "p: " << p << std::endl;
            std::cout << "k: " << k << std::endl;
            std::cout << "n: " << n << std::endl;
            std::cout << "tid: " << tid << std::endl;
            
            // Печатаем массив, если он не нулевой
            if (array != nullptr) {
                std::cout << "array: ";
                for (int i = 0; i < m; ++i) {
                    std::cout << array[i] << " "; // Установка точности
                }
                std::cout << std::endl;
            } else {
                std::cout << "array: nullptr" << std::endl;
            }

            std::cout << "m: " << m << std::endl;
            std::cout << "amount_of_changed: " << amount_of_changed << std::endl;

            std::cout << "pan: " << pan << std::endl;
            std::cout << "na0: " << na0 << std::endl;
            std::cout << "na1: " << na1 << std::endl;
            printf("\n");
        }
        
};
void* thread_func(void *arg);
void UpdateElements(Args* a);


#endif