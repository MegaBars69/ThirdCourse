#include <iostream>
#include <fstream>
#include <pthread.h>
#include <cmath>
#include <stdio.h>
#include "func.hpp"

double get_cpu_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);

    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec/1e6;
}

bool NumberAndNextArePrime(unsigned long long int num)
{
    unsigned long long int next = num + 6;
    unsigned long long int i;
    if (num%2 == 0)
    {
        return false;
    }
    else
    {
        for (i = 3; i*i <= num; i += 2)
        {
            if (num % i == 0 || next % i == 0)
            {
                return false;
            }
        }
        for ( ; i*i <= next; i += 2)
        {
            if (next % i == 0)
            {
                return false;
            }
        }
    }
    return true;
}

void FindInterval(Args* a)
{
    int n = a->n;
    //int k = a->k;
    int p = a->p;
    int chunk = a->chunk;
    int founded_pairs = a->founded_pairs;
    int pairs_in_interval = 0;
    
    unsigned long long int L = a->L;
    unsigned long long int i = 0;
    unsigned long long int up_bound = L + chunk;
    unsigned long long int answer = 0;

    while (founded_pairs < n)  
    {
        pairs_in_interval = 0;

        for (i = L; i < up_bound; i++)
        {
            if (NumberAndNextArePrime(i))
            {
                pairs_in_interval++;
            }

        } 

        reduce_sum(p, &pairs_in_interval, 1);
        founded_pairs += pairs_in_interval;
        L += p*chunk;
        up_bound = L + chunk; 
    }

    a->L = L-p*chunk;
    a->founded_pairs = founded_pairs - pairs_in_interval;
    return;
}

void ProccesResult(Args* a)
{
    int p = a->p;
    unsigned long long int i = a->L;
    unsigned long long int up_bound = i + p*a->chunk; 
    int pair_in_interval = a->founded_pairs;
    int n = a->n; 
    for (; i < up_bound; i += 2)
    {
        if (NumberAndNextArePrime(i))
        {
            pair_in_interval++;
            if (pair_in_interval == n)
            {
                a->answer = i + 6;
            }
        }
        
    }
    
} 


void* thread_func(void *arg)
{
    Args *a = (Args *)arg;
    double t;
    
    t = get_cpu_time();

    FindInterval(a);

    if (a->k == 0)
    {
        ProccesResult(a);
    }

    t = get_cpu_time() - t;
    a->cpu_time_of_thread = t;
    a->cpu_time_of_all_threads = t;
    reduce_sum(a->p, &a->cpu_time_of_all_threads,1);

    return nullptr;    
}