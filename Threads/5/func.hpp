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
    error_few_el,
    error_many_el,
    succes
};

class Args{
    public:
        int p = 0;
        int k = 0;
        int n = 0;
        pthread_t tid = 0;
        double cpu_time_of_thread = 0;
        double cpu_time_of_all_threads = 0;

        double* array = nullptr;
        char* name = nullptr;
        io_status reading_file = io_status::undef;
        int m = 0;
        int l = 0;

        int amount_of_changed = 0;
        double prev = 0;  //prev array's elements
        double next = 0;  //next array's elemnets
        double q_left = 0;
        double q_right = 0;
        double q = 0;
        double left_sum = 0;
        double right_sum = 0;
        int el_in_left_sum = 1;
        int el_in_right_sum = 1;

        bool left_can_connect = false;
        bool all_can_connect = true;
        int res = 0;


        void PrintAll(bool with_array = false) const 
        {
            printf("\n");
            std::cout << "Args Data:" << std::endl;
            std::cout << "p: " << p << std::endl;
            std::cout << "k: " << k << std::endl;
            std::cout << "n: " << n << std::endl;
            printf("cpu time of thread: %8.2e \n", cpu_time_of_thread);
            printf("cpu time of all threads: %8.2e \n", cpu_time_of_all_threads);

            std::cout << "tid: " << tid << std::endl;
            
            // Печатаем массив, если он не нулевой
            if(with_array)
            {
                if (array != nullptr) {
                    std::cout << "array: ";
                    for (int i = 0; i < m; ++i) {
                        std::cout << array[i] << " "; // Установка точности
                    }
                    std::cout << std::endl;
                } else {
                    std::cout << "array: nullptr" << std::endl;
                }
            }
            std::cout << "m: " << m << std::endl;
            std::cout << "amount of changed: " << amount_of_changed << std::endl;

            std::cout << "prev: " << prev << std::endl;
            std::cout << "next: " << next << std::endl;
            std::cout << "rigth sum: " << right_sum << std::endl;
            std::cout << "left sum: " << left_sum << std::endl;

            std::cout << "el in rigth sum: " << el_in_right_sum << std::endl;
            std::cout << "el in left sum: " << el_in_left_sum << std::endl;
            std::cout << "left can connect: " << left_can_connect << std::endl;

            std::cout << "q right: " << q_right << std::endl;
            std::cout << "q left: " << q_left << std::endl;
           
            printf("\n");
        }
        
};
void* thread_func(void *arg);
void ProccesElements(Args* a);
int ProccesResults(Args *a);

#endif