#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"
#define EPSILON pow(10, -16)
/*
template <typename T>
void reduce_sum(int p,T* a, int n)
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
}*/

io_status sum_of_increasing(FILE* file, double* result, int* amount)
{
    double number;
    int count = 0;
    double sum = 0;
    bool begin_of_file = true;
    bool begin_of_seq = true;
    double prev_num = std::numeric_limits<double>::lowest();

    while (true) 
    {
        int res = fscanf(file, "%lf", &number);
        if(res == EOF)
        {
            *result = sum;
            *amount = count;
            break;
        }
        if (res== 0)
        {
            return io_status :: error_read;
        }
        if (!begin_of_file)
        {
            if (number > prev_num) 
            {
                count += (begin_of_seq ? 2 : 1);
                sum += (begin_of_seq ? prev_num + number: number);
                begin_of_seq = false;
            } 
            else if ((fabs(number - prev_num) <= EPSILON) || (number < prev_num)) {
                begin_of_seq = true;
            }
        }
        else
        {
            begin_of_file = false;
        }
        prev_num = number;

    }
    *result = sum;
    *amount = count;
    return io_status::succes;
}

io_status amount_of_elements(FILE* file, double* result, double el)
{
    double number;
    double sum = 0;

    while (true) 
    {
        int res = fscanf(file, "%lf", &number);
        if(res == EOF)
        {
            *result = sum;
            break;
        }
        if (res== 0)
        {
            return io_status :: error_read;
        }

        if (number > el)
        {
            sum++;
        }
    } 

    *result = sum;
    return io_status::succes;
}

io_status avarage(double* result, double* sum, int* amount)
{
    if (*amount == 0)
    {
        return io_status :: error_av_doesnt_exist;
    }

    *result = *sum/(*amount);
    return io_status::succes;
}

void* thread_func(void *arg)
{
    Args *a = (Args *)arg;
    FILE *fp = fopen(a->filename, "r");

    if(fp == nullptr)
    {
        a->error_type = io_status :: error_open;
        a->error_flag = 1;
    }
    else
    {
        a->error_flag = 0;
    }

    reduce_sum(a->p ,&a->error_flag, 1);

    if(a->error_flag > 0)
    {
        if (fp != nullptr)
        {
            fclose(fp);
        }
        return nullptr;
    }
    
    a->error_type = sum_of_increasing(fp, &a->res, &a->count);

    if (a->error_type != io_status::succes)
    {
        a->error_flag = 1; 
    }
    else
    {
        a->error_flag = 0;
    }

    reduce_sum(a->p, &a->error_flag, 1);

    if(a->error_flag > 0)
    {
        if (fp != nullptr)
        {
            fclose(fp);
        }
        return nullptr;
    }

    fclose(fp);

    a -> error_type = io_status::succes;
    reduce_sum(a->p, &a->count, 1);

    //reduce_sum(a->p, &a->res, 1);
    //Second part. Finding avarage

    a->error_type = avarage(&a->avarage, &a->res, &a->count);

    if (a->error_type != io_status::succes)
    {
        a->error_flag = 1; 
    }
    else
    {
        a->error_flag = 0;
    }

    reduce_sum(a->p, &a->error_flag, 1);
    if(a->error_flag > 0)
    {
        a -> error_type = io_status::error_av_doesnt_exist;
        return nullptr;
    }

    a -> error_type = io_status::succes;
    reduce_sum(a->p, &a->avarage, 1);

    //Fird Part
    fp = fopen(a->filename, "r");
    if(fp == nullptr)
    {
        a->error_type = io_status :: error_open;
        a->error_flag = 1;
    }
    else
    {
        a->error_flag = 0;
    }

    reduce_sum(a->p ,&a->error_flag, 1);

    if(a->error_flag > 0)
    {
        if (fp != nullptr)
        {
            fclose(fp);
        }
        return nullptr;
    }
    
    a->error_type = amount_of_elements(fp, &a->res, a->avarage);

    if (a->error_type != io_status::succes)
    {
        a->error_flag = 1; 
    }
    else
    {
        a->error_flag = 0;
    }

    reduce_sum(a->p, &a->error_flag, 1);

    if(a->error_flag > 0)
    {
        if (fp != nullptr)
        {
            fclose(fp);
        }
        return nullptr;
    }

    fclose(fp);

    a -> error_type = io_status::succes;
    /*if (a->k == 0)
    {
        std::cout<<a->avarage<<std::endl;
    }*/
    
    reduce_sum(a->p, &a->res, 1);

    return nullptr;    
}

int procces_args(Args *a)
{

    for (int i = 0; i < a->p; i++)
    {
        if (a[i].error_type == io_status::error_read )
        {
            printf("Error in file %s \n", a[i].filename);
            printf("Result = %d \n", -2);
            return -2;
        }
        else if (a[i].error_type == io_status::error_open)
        {
            printf("Error opening file %s \n", a[i].filename);
            printf("Result = %d \n", -1);
            return -1;
        }
        else if (a[i].error_type == io_status::error_av_doesnt_exist)
        {
            printf("Error avarage number doesn't exist \n");
            printf("Result = %d \n", -3);
            return -3;
        }
    }

    return 0;
}