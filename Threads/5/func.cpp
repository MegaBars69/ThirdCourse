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
    int p = a->p;
    double sum = pa[0];
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
    } 
    if(k == 0 && m >= 2)
    {
        if(fabs(pa[0]) > EPSILON)
        { 
            prev_q = pa[1]/pa[0]; 
            sum += pa[1];
            el_in_sum = 2;
            if (m == 2)
            {
                if (p > 1 && fabs(pa[1]) >EPSILON && fabs(a[0].next/pa[1] - prev_q) < EPSILON)
                {
                    a->el_in_right_sum = el_in_sum;
                    a->right_sum = sum;
                    a->q_right = prev_q;
                }
                else if (p > 1)
                {
                    a->el_in_right_sum = 1;
                    a->right_sum = pa[1];
                    a->q_right = 0;
                }
                return;
            }
        }
        else if(m == 2)
        {
            a->el_in_right_sum = 1;
            a->right_sum = pa[1];
            return;
        }
        
    }
    else if (m == 2 && k > 0)
    {
        if (fabs(pa[0]) > EPSILON)
        { 
            if (fabs(prev_q - pa[1]/pa[0]) < EPSILON)
            {
                a->right_sum = a->left_sum = pa[0] + pa[1];
                a->el_in_right_sum = a->el_in_left_sum = 2;
                a->left_can_connect = true;
                a->q_right = a->q_left = pa[1]/pa[0];
            }
            else if (k < p-1)
            {
                if (fabs(pa[1]) > EPSILON)
                {
                    if (fabs(a->next/pa[1] - pa[1]/pa[0])< EPSILON)
                    {
                        a->right_sum = pa[0] + pa[1];
                        a->el_in_right_sum = 2;
                        a->q_left = a->q_right = pa[1]/pa[0];
                        a->left_sum = pa[0];
                        a->el_in_left_sum = 1;
                    }
                    else
                    {
                        a->right_sum = pa[1];
                        a->q_left = a->q_right = pa[1]/pa[0];
                        a->left_sum = pa[0];
                        a->el_in_right_sum= a->el_in_left_sum = 1;
                    }
                    
                }
                
            }
            
        
            else
            {
                a->right_sum = pa[1];
                a->left_sum = pa[0];
                a->el_in_right_sum = a->el_in_left_sum = 1;

                a->q_right = a->q_left = pa[1]/pa[0];
            }
        }
        return;
    }
    
    
    else if (m == 1)
    {
        a->left_sum = a->right_sum = pa[0];
        a->el_in_right_sum = a->el_in_left_sum = 1; 
        return;
        //a->q_right = (fabs(pa[0]) > EPSILON && k < p-1 ? a->next/pa[0] : 0);
    }
    
    for (int i = (k > 0 ? 1 : 2); i < m; i++)
    {
        prev_el = pa[i-1];
        cur_el = pa[i];

        //cur_q = (fabs(prev_el) < EPSILON ? 0 : cur_el/prev_el);
        if(fabs(prev_el) > EPSILON)
        {
            cur_q =  cur_el/prev_el;
            if (fabs(cur_q - 1) > EPSILON)
            {
            
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
                    a->all_can_connect = false;
                    in_seq = false;
                }
                else if (start_of_seq < 0 && el_in_sum >= 2 && fabs(cur_q - prev_q) > EPSILON)
                {
                    a->left_can_connect = true;
                    a->left_sum = sum;
                    a->el_in_left_sum = el_in_sum;
                    a->q_left = cur_q;
                    sum = cur_el;
                    el_in_sum = 1;
                    in_seq = false;
                    a->all_can_connect = false;
                }
                else
                {
                    sum = cur_el + prev_el;
                    el_in_sum = 2;
                }
                
                if (i == m-1 && (k < p-1) && fabs(cur_q - prev_q) < EPSILON)
                {
                    if (fabs(cur_el)> EPSILON && fabs(a->next/cur_el - cur_q) < EPSILON)
                    {
                        if (start_of_seq < 0)
                        {
                            a->left_can_connect = true;
                            a->left_sum = sum;
                            a->el_in_left_sum = el_in_sum;
                            a->q_left = cur_q;
                        }
                        a->right_sum = sum ;
                        a->el_in_right_sum = el_in_sum;
                        a->q_right = cur_q;
                    }
                    else if (start_of_seq >= 0)
                    {
                        new_el = sum/el_in_sum;
                        for (int j = start_of_seq; j <= i ; j++)
                        {
                            pa[j] = new_el;
                            amount_of_changed++;
                        }
                    }
                    else if (start_of_seq < 0)
                    {
                        a->left_can_connect = true;
                        a->right_sum = a->left_sum = sum;
                        a->el_in_right_sum = a->el_in_left_sum = el_in_sum;
                        a->q_right = a->q_left = cur_q;
                    }
                    
                }
                else if (i == m-1 && (k < p-1))
                {
                    if (fabs(pa[m-1]) > EPSILON)
                    {
                        if (fabs(a->next/pa[m-1] - cur_q) < EPSILON)
                        {
                            if (start_of_seq < 0 && in_seq)
                            {
                                a->left_can_connect = true;
                                a->left_sum = sum;
                                a->el_in_left_sum = el_in_sum;
                                a->q_left = cur_q;
                            }
                            
                            a->right_sum = sum;
                            a->el_in_right_sum = el_in_sum;
                            a->q_right = cur_q;
                        }
                        else
                        {
                            a->right_sum = cur_el;
                            a->el_in_right_sum = 1;
                            a->q_right = cur_q; 
                        }
                        
                    }

                }
                else if (i == m-1 && fabs(cur_q - prev_q) < EPSILON && (k == p-1) && start_of_seq >= 0)
                {
                    /*sum += cur_el;
                    el_in_sum++;*/
                    new_el = sum/el_in_sum;
                    for (int j = start_of_seq; j <= i ; j++)
                    {
                        pa[j] = new_el;
                        amount_of_changed++;
                    }
                    sum = cur_el;
                    el_in_sum = 1;
                    in_seq = false;
                    a->all_can_connect = false; 
                }
                else if (i == m-1 && fabs(cur_q - prev_q) < EPSILON && (k == p-1) && start_of_seq < 0)
                {
                    a->left_can_connect = true;
                    a->left_sum = sum;
                    a->right_sum = sum;
                    a->el_in_right_sum = el_in_sum;
                    a->el_in_left_sum = el_in_sum;
                    a->q_right = a->q_left = cur_q;
                }
            }    
            
        }
        
        prev_q = cur_q;
        
    }

    a->amount_of_changed = amount_of_changed; 
}

int ProccesResults(Args *a)
{
    //int k = a->k;
    int p = a->p;
    int result = a[0].amount_of_changed; 
    double sum = 0;
    int el_in_sum = 0;
    int start_of_seq = 0;
    //double prev_q = 0;
    //double cur_q = 0;
    double new_el = 0;
    bool in_seq = false;
    
    for (int i = 1; i < p; i++)
    {
        a[i].PrintAll(true);

        if (a[i].m == 1 && i >= 1 && fabs(a[i].prev) > EPSILON)
        {
            a[i].q_right = a[i].array[0]/a[i].prev;
        }
        
        if (a[i].left_can_connect)
        {
            if (!in_seq)
            {
                in_seq = true;
                start_of_seq = i-1 - (i >= 2 && fabs(a[i-1].prev) > EPSILON && a[i-1].m == 1 && fabs(a[i].prev/a[i-1].prev - a[i].q_left) < EPSILON ? 1 : 0);
                sum += a[i-1].right_sum + (i >= 2 && fabs(a[i-1].prev) > EPSILON && a[i-1].m == 1 && fabs(a[i].prev/a[i-1].prev - a[i].q_left) < EPSILON? a[i-1].prev : 0);
                el_in_sum += a[i-1].el_in_right_sum + (i >= 2 && fabs(a[i-1].prev) > EPSILON && a[i-1].m == 1 && fabs(a[i].prev/a[i-1].prev - a[i].q_left) < EPSILON? 1 : 0);
            }

            sum += a[i].left_sum;
            el_in_sum += a[i].el_in_left_sum;

            /*if (fabs(a[i].q_left - a[i-1].q_right) <= EPSILON)
            {
                sum += a[i-1].right_sum + (a[i-1].m == 1 ? a[i-1].prev : 0);
                el_in_sum += a[i-1].el_in_right_sum + (a[i-1].m == 1 ? a[i-1].prev : 0);
            }
            else
            {
                sum += a[i].prev;
                el_in_sum++;
                a[i-1].el_in_right_sum = 1;
            }*/

            if (a[i].el_in_left_sum < a[i].m || (fabs(a[i].array[a[i].m-1]) > EPSILON && fabs(a[i].next/a[i].array[a[i].m-1] - a[i].q_left) > EPSILON))
            {
                /*sum += a[i].left_sum;
                el_in_sum += a[i].el_in_left_sum;*/
                new_el = sum/el_in_sum;
                for (int j = start_of_seq;  j < i; j++)
                {
                    for (int s = a[j].m - a[j].el_in_right_sum; s < a[j].m; s++)
                    {
                        a[j].array[s] = new_el;
                        result++;
                    }
                }
                for (int s = 0; s < a[i].el_in_left_sum; s++)
                {
                    a[i].array[s] = new_el;
                    result++;
                }
                el_in_sum = 0;
                sum = 0;
                in_seq = false;
            } 
        } 
          
        else 
        {
            if (fabs(a[i].prev) > 0)
            {
                if (fabs(a[i].array[0]/a[i].prev - a[i-1].q_right) < EPSILON)
                {
                    if (!in_seq)
                    {
                        in_seq = true;
                        start_of_seq = i-1 - (a[i-1].m == 1 ? 1 : 0);
                        sum += a[i-1].right_sum + (a[i-1].m == 1 ? a[i-1].prev : 0);
                        el_in_sum += a[i-1].el_in_right_sum + (a[i-1].m == 1 ? 1 : 0);
                    }
                    sum += a[i].left_sum;
                    el_in_sum += a[i].el_in_left_sum;
                    
                }
                else if(in_seq && el_in_sum >= 1)
                {
                    new_el = sum/el_in_sum;
                    for (int j = start_of_seq;  j < i; j++)
                    {
                        for (int s = a[j].m - a[j].el_in_right_sum; s < a[j].m; s++)
                        {
                            a[j].array[s] = new_el;
                            result++;
                        }
                    }
                    el_in_sum = 0;
                    sum = 0;
                    in_seq = false;
                } 
                
                if (a[i].m > 1 && el_in_sum >= 1)
                {
                    new_el = sum/el_in_sum;
                    a[i].array[0] = new_el;
                    result++;
                    for (int j = start_of_seq;  j < i; j++)
                    {
                        for (int s = a[j].m - a[j].el_in_right_sum; s < a[j].m; s++)
                        {
                            a[j].array[s] = new_el;
                            result++;
                        }
                    }
                    sum = 0;
                    el_in_sum = 0;
                    in_seq = false;
                }
                
            }
          
        }
        if (i == p-1 && in_seq && el_in_sum >= 1)
        {
            new_el = sum/el_in_sum;
            for (int j = start_of_seq;  j <= i; j++)
            {
                for (int s = a[j].m - a[j].el_in_right_sum; s < a[j].m; s++)
                {
                    a[j].array[s] = new_el;
                    result++;
                }
            }
            el_in_sum = 0;
            sum = 0;
            in_seq = false;
        }
    }
    return result;
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


