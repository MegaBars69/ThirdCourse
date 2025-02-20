#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"
#define EPSILON pow(10, -16)

io_status amount_of_elements(FILE* file, double* Sum, int* Count)
{
    double number;
    double sum = 0;
    int amount = 0;
    while (true) 
    {
        int res = fscanf(file, "%lf", &number);
        if(res == EOF)
        {
            *Sum = sum;
            *Count = amount;
            break;
        }
        if (res== 0)
        {
            return io_status :: error_read;
        }
        sum += number;
        amount++;        
    } 

    *Sum = sum;
    *Count = amount;

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

io_status amount_of_local_min(FILE* file, Args* Result, double avarage)
{
    double current = 0;
    double prev = 0;
    double next = 0;
    bool first = true;
    bool second = true;
    double sum = 0;
    int amount_ = 0;

    while (true) 
    {
        int res = fscanf(file, "%lf", &next);
        if(res == EOF)
        {
            if (!second)
            {
                Result->amount_of_min = sum;
                Result->last = current;
                if (prev >= current && current < avarage)
                {
                    Result->last_can_be_a_minimum = true;
                }
                else
                {
                    Result->last_can_be_a_minimum = false;
                }
            }

            break;
        }
        if (res== 0)
        {
            return io_status :: error_read;
        }
        if (first)
        {
            first = false;
            Result->first = next;
            Result->last = next;
            current = next;
        }
        else if (second)
        {
            second = false;
            Result->last = next;

            if (next >= current && current < avarage)
            {
                Result->first_can_be_a_minimum = true; 
            }
            if (next <= current && next < avarage)
            {
                Result->last_can_be_a_minimum = true;
            }
            
            prev = current;
            current = next;
        }
        else
        {
            if (next >= current && prev >= current && current < avarage)
            {
                sum++;
            }
            prev = current;
            current = next;
        }
        amount_ ++;          
    } 
    Result->count = amount_;
    Result->amount_of_min = sum;

    return io_status::succes;
}

int procces_results(Args* r, int p)
{
    int amount_of_min = r[0].amount_of_min;
    double avarage = r[0].avarage;

    if (r[0].error_type == io_status::error_read )
        {
            printf("Error in file %s \n", r[0].filename);
            printf("Result = %d \n", -2);
            return -2;
        }
        else if (r[0].error_type == io_status::error_open)
        {
            printf("Error opening file %s \n", r[0].filename);
            printf("Result = %d \n", -1);
            return -1;
        }
        else if (r[0].error_type == io_status::error_av_doesnt_exist)
        {
            printf("Error avarage number doesn't exist \n");
            printf("Result = %d \n", -3);
            return -3;
        }
    for (int k = 1; k < p; k++)
    {
        if (r[k].error_type == io_status::error_read )
        {
            printf("Error in file %s \n", r[k].filename);
            printf("Result = %d \n", -2);
            return -2;
        }
        else if (r[k].error_type == io_status::error_open)
        {
            printf("Error opening file %s \n", r[k].filename);
            printf("Result = %d \n", -1);
            return -1;
        }
        else if (r[k].error_type == io_status::error_av_doesnt_exist)
        {
            printf("Error avarage number doesn't exist \n");
            printf("Result = %d \n", -3);
            return -3;
        }
        if (r[k].count == 1)
        {
            if (r[k-1].last_can_be_a_minimum && r[k-1].last <= r[k].first)
            {
                amount_of_min++;
            }
            if (k + 1 < p)
            {
                if (r[k+1].first >= r[k].first && r[k-1].last >= r[k].first && r[k].first < avarage)
                {          
                    amount_of_min++;
                }
                
            }
            
        }
        else if (r[k].count == 0)
        {
            r[k].last_can_be_a_minimum = r[k-1].last_can_be_a_minimum;
            r[k].last = r[k-1].last;
        }
        
        else
        {
            amount_of_min += r[k].amount_of_min;
            
            if (r[k-1].last_can_be_a_minimum && r[k-1].last <= r[k].first)
            {
                amount_of_min++;
            }
            if (r[k].first_can_be_a_minimum && r[k-1].last >= r[k].first)
            {
                amount_of_min++;
            }
        }
    }
    return amount_of_min;
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
    
    a->error_type = amount_of_elements(fp, &a->sum, &a->count);

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

    a->error_type = avarage(&a->avarage, &a->sum, &a->count);
    
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
    a->count = a->count;
    a->avarage = a->avarage;
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
    
    a->error_type = amount_of_local_min(fp, a, a->avarage);

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
    //std::cout<<a->amount_of_min<<std::endl;
    reduce_sum<int>(a->p, nullptr, 0);


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