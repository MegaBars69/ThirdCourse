#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <time.h>
#include <sys/time.h>
#include <unistd.h> 
#include <sys/resource.h>
#include <cmath>
#include "func.hpp"
#include <sys/resource.h>
#include <sys/sysinfo.h>

#define EPSILON pow(10, -15)

void ProccesElements(Args* a)
{
    double* pa = a->array;
    int amount_of_changed = 0;
    int m = a->m;
    int k = a->k;
    double sum = 0;
    int el_in_sum = 1;
    double cur_q = 0;
    double prev_q = 0;
    double cur_el, prev_el;
    int start_of_seq = 0;
    double new_el = 0;
    bool in_seq = false;
   
    if (k > 0 && fabs(a->prev) > EPSILON)
    {
        prev_q = a->q = pa[0]/a->prev;
        sum = pa[0];
    } 
    if(k == 0 && m >= 2)
    {
        prev_q = (fabs(pa[0]) < EPSILON ? 0 : pa[1]/pa[0]); 

        sum = pa[0] + pa[1];
        el_in_sum = 2;
    }
    
    for (int i =(k > 0 ? 1 : 2); i < m; i++)
    {
        prev_el = pa[i-1];
        cur_el = pa[i];

        //cur_q = (fabs(prev_el) < EPSILON ? 0 : cur_el/prev_el);
        if(fabs(prev_el) > EPSILON)
        {
            cur_q =  cur_el/prev_el;
        
            if (fabs(cur_q - prev_q) < EPSILON)
            {
                sum += cur_el;
                el_in_sum++;
                if (el_in_sum == 2 || (!in_seq && el_in_sum == 3))
                {
                    start_of_seq = i-2; // Smth may go wrong
                    in_seq = true;
                }              
            }
            else if (start_of_seq >= 0 && el_in_sum > 2 && fabs(cur_q - prev_q) > EPSILON)
            {
                //End of seq
                new_el = sum/el_in_sum;
                for (int j = start_of_seq; j < i ; j++)
                {
                    pa[j] = new_el;
                    amount_of_changed++;
                }
                sum = cur_el;
                el_in_sum = 1;
                in_seq = false;
                a->all_is_seq = false;
            }
            else if (start_of_seq < 0 && el_in_sum > 2 && fabs(cur_q - prev_q) > EPSILON)
            {
                a->left_can_connect = true;
                a->left_sum = sum;
                a->el_in_left_sum = el_in_sum;
                a->q_left = cur_q;
                sum = cur_el;
                el_in_sum = 1;
                in_seq = false;
                a->all_is_seq = false;
            }
            else
            {
                sum = cur_el + prev_el;
                el_in_sum = 2;
            }
            
            if (i == m-1 && fabs(cur_q - prev_q) < EPSILON)
            {
                a->right_can_connect = true;
                a->right_sum = sum;
                a->el_in_right_sum = el_in_sum;
                a->q_right = cur_q;
            }
        }
        
        prev_q = cur_q;
    }

    a->amount_of_changed = amount_of_changed; 
}
/*
void ProccesResults(Args *a)
{
    int k = a->k;
    int p = a->p;
    
}*/

double get_cpu_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);

    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec/1e6;
}

void* thread_func(void *arg)
{
    Args *a = (Args *)arg;
    double t = 0;
    //int k = a->k;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    //int n_cpus = get_nprocs();
    //int cpu_id = n_cpus - 1 -(k%n_cpus);
    pthread_t tid = a->tid;
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
    
    t = get_cpu_time();

    ProccesElements(a);

    reduce_sum(a->p, &a->amount_of_changed, 1);


    t = get_cpu_time() - t;

    a->cpu_time_of_thread = t;
    a->cpu_time_of_all_threads = t;

    reduce_sum(a->p, &a->cpu_time_of_all_threads, 1);
  

    return nullptr;    
}


