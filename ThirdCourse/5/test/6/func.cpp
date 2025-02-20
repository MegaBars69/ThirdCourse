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
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include "func.hpp"

#define EPSILON pow(10, -16)

void get_avarage(Args *a)
{
    int n1 = a->m1;
    int n2 = a->m2;
    int k = a->k;
    int p = a->p;
    double *A = a->A;
    double *prev_line = a->prev_line;
    double *next_line = a->next_line;

    int amount_of_el = 0;
    double sum = 0;
    double mean_cur_sum = 0;
    double cur_el;
    if ((k == 0 || k == p - 1) && n1 == 1)
    {
        return;
    
    }
    for (int i = (k > 0 ? 0 : 1); i < (k < p-1 ? n1 : n1-1); i++)
    {
        for (int j = 0; j < n2 - 1; j++)
        {
            cur_el = A[i*n2 + j];
            
            if (i == 0 && k > 0 && a->n1 > 2)
            {
                mean_cur_sum = (prev_line[j] + prev_line[j+1] + A[j + 1] + (n1 > 1 ? A[(i+1)*n2 + j] + A[(i+1)*n2 + j + 1] : next_line[j] + next_line[j + 1])) / 5.0;
                if (cur_el >= mean_cur_sum)
                {
                    sum += cur_el;
                    amount_of_el++;                    
                    //std::cout<<"("<<i<<", "<<j<<")"<<std::endl;

                }
            }
            else if (i > 0 && i + 1 < n1 && a->n1 > 2 && i < n1 - 1)
            {
                mean_cur_sum = (A[(i-1)*n2 + j] + A[(i-1)*n2 + j + 1] + A[i*n2 + j + 1] +  A[(i+1)*n2 + j] + A[(i+1)*n2 + j + 1]) / 5.0;
                if (cur_el >= mean_cur_sum)
                {
                    sum += cur_el;
                    amount_of_el++;                    
                    //std::cout<<"("<<i<<", "<<j<<")"<<std::endl;

                }
            }
            else if (i == n1 - 1 && k < p - 1)
            {
                mean_cur_sum = (A[(i-1)*n2 + j] + A[(i-1)*n2 + j + 1] + A[i*n2 + j + 1] + next_line[j] + next_line[j + 1]) / 5.0;
                if (cur_el >= mean_cur_sum)
                {
                    sum += cur_el;
                    amount_of_el++;
                                        //std::cout<<"("<<i<<", "<<j<<")"<<std::endl;

                }
            }
        }
        
    }
    
    a->sum = sum;
    a->amount_of_el = amount_of_el;
}

void UpdateElements(Args *a)
{
    int n1 = a->m1;
    int n2 = a->m2;
    int k = a->k;
    int p = a->p;
    double *A = a->A;
    double *prev_line = a->prev_line;
    double *next_line = a->next_line;

    double mean_cur_sum = 0;
    double cur_el;
    double new_el = a->avarage;
    int amount_of_changed = 0; 

    if ((k == 0 || k == p - 1) && n1 == 1)
    {
        return;
    }
    
    printf("\n");
    for (int i = (k > 0 ? 0 : 1); i < (k < p-1 ? n1 : n1-1); i++)
    {
        if (k == 0 && i == 1)
        {
            memcpy(prev_line, A, n2*sizeof(double));
        }
        
        for (int j = 0; j < n2; j++)
        {
            cur_el = A[i*n2 + j];

            if (j + 1 < n2 && i == 0 && k > 0 && a->n1 > 2)
            {
                mean_cur_sum = (prev_line[j] + prev_line[j+1] + A[j + 1] + (n1 > 1 ? A[(i+1)*n2 + j] + A[(i+1)*n2 + j + 1] : next_line[j] + next_line[j + 1])) / 5.0;
                if (cur_el >= mean_cur_sum)
                {
                    A[i*n2 + j] = new_el;
                    amount_of_changed++;
                    //std::cout<<"("<<i<<", "<<j<<")"<<std::endl;
                }
            }
            else if (j + 1 < n2 && i > 0 && i + 1 < n1 && a->n1 > 2 && i < n1 - 1)
            {
                mean_cur_sum = (prev_line[j] + prev_line[j+1] + A[i*n2 + j + 1] +  A[(i+1)*n2 + j] + A[(i+1)*n2 + j + 1]) / 5.0;
                if (cur_el >= mean_cur_sum)
                {
                    A[i*n2 + j] = new_el;
                    amount_of_changed++;
                    //std::cout<<"("<<i<<", "<<j<<")"<<std::endl;
                }
            }
            else if (j + 1 < n2 && i == n1 - 1 && k < p - 1)
            {
                mean_cur_sum = (prev_line[j] + prev_line[j+1] + A[i*n2 + j + 1] + next_line[j] + next_line[j + 1]) / 5.0;
                if (cur_el >= mean_cur_sum)
                {
                    A[i*n2 + j] = new_el;
                    amount_of_changed++;
                    //std::cout<<"("<<i<<", "<<j<<")"<<std::endl;
                }
            }
            prev_line[j] = cur_el;
        }
        
    }
    a->amount_of_changed = amount_of_changed;

} 

void PrintMatrix(double* matrix, int n1, int n2) 
{
    printf("\n");
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            printf("%10.3e ", *(matrix + i * n2 + j)); 
        }
        printf("\n");
    }
    printf("\n");
}

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
    int p = a->p;
    
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    //int n_cpus = get_nprocs();
    //int cpu_id = n_cpus - 1 -(k%n_cpus);
    pthread_t tid = a->tid;
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    t = get_cpu_time();

    get_avarage(a);
    reduce_sum(p, &a->sum, 1);
    reduce_sum(p, &a->amount_of_el, 1);
    if (a->amount_of_el <= 0)
    { 
        a->res = 1;
    }
    else
    {
        a->avarage = a->sum/a->amount_of_el;
        UpdateElements(a);
    }

    t = get_cpu_time() - t;

    a->cpu_time_of_thread = t;
    a->cpu_time_of_all_threads = t;

    reduce_sum(a->p, &a->cpu_time_of_all_threads, 1);
  
    reduce_sum(a->p, &a->amount_of_changed, 1);

    return nullptr;    
}


