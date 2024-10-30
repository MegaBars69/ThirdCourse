#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"
#define EPSILON pow(10, -16)

// Функция для работы потоков
void* count_max_elements(void* arg) 
{
    Thread* data = static_cast<Thread*>(arg);
    const std::string& filename = data->filename;
    int& count_max = data->count_max;

    std::ifstream file(filename);
    if (!file.is_open()) {
        return reinterpret_cast<void*>(-1); // Ошибка открытия файла
    }

    double number;
    double max_value = std::numeric_limits<double>::lowest();
    count_max = 0;

    while (true) {
        file>>number;

        if (file.fail()) {
            if (file.eof()) {
                // Достигнут конец файла
                break;
            } 
            else {
                // Ошибка ввода (например, неверный тип данных)
                file.close();
                return reinterpret_cast<void*>(-2); // Ошибка чтения элемента
            }
        }        
        if (number > max_value) {
            max_value = number;
            count_max = 1;
        } else if (fabs(number-max_value) <= EPSILON) {
            count_max++;
        }
    }

    file.close();
    return nullptr; // Успешное завершение
}

void reduce_sum(int p, double* a, int n)
{
    static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double* r = nullptr;
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

io_status amount_of_max(FILE* file, double* result)
{
    double number;
    double max_value = std::numeric_limits<double>::lowest();
    int count = 0;

    while (true) 
    {
        int res = fscanf(file, "%lf", &number);
        if(res == EOF)
        {
            *result = count;
            break;
        }
        if (res== 0)
        {
            return io_status :: error_read;
        }
        
        if (number > max_value) {
            max_value = number; // Обновляем максимальное значение
            count = 1; // Сбрасываем счетчик
        } else if (fabs(number - max_value) <= EPSILON) {
            count+=1; // Увеличиваем счетчик, если число равно максимальному
        }
    }
    *result = count;
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
    
    a->error_type = amount_of_max(fp, &a->res);

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

    reduce_sum(a->p, &a->res, 1);
    
    return nullptr;    
}

int procces_args(Args *a)
{
    /*if (a->error_type != io_status :: succes)
    {*/
    for (int i = 0; i < a->p; i++)
    {
        if (a[i].error_type == io_status::error_read )
        {
            printf("Error in file %s \n", a[i].filename);
            return 1;
        }
        else if (a[i].error_type == io_status::error_open)
        {
            printf("Error opening file %s \n", a[i].filename);
            return 1;
        }
        else if (a[i].error_type == io_status::undef)
        {
            printf("Somthing went wrong %s \n", a[i].filename);
            return 1;
        }
    }
    //}
    return 0;
}