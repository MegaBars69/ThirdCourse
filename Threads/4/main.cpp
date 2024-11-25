#include <iostream>
#include <fstream>
#include <pthread.h>
#include <cmath>
#include "func.hpp"

int main(int argc, char* argv[])
{
    Args* a;
    double *array;
    double *pa = nullptr;
    int k = 0;

    if(argc < 4)
    {
        printf("RESULT :");

        printf("Usage ./a.out p n filename\n");
        return 1;
    }
    int n = atoi(argv[2]);
    int p = atoi(argv[1]); 

    if (p > n)
    {
        p = n; 
    }
    if (p <= 0)
    {
        printf("RESULT %2d:",p);

        printf(" p > 0 \n");
        return 4;
    }
    
    int m = n/p;
    int l = n%p;

    array = new double[n];     

    a = new Args[p];
    pa = array;
    double el;
    int i = 0;

    FILE *fp = fopen(argv[3], "r");
    if(fp == nullptr)
    {
        printf("RESULT %2d:",p);

        printf("Can't open file \n");
        delete[] array;
        delete[] a;

        return 1;
    }
    while (true) 
    {
        int res = fscanf(fp, "%lf", &el);
        if(res == EOF)
        {  
            if(i != n)
            {
                printf("RESULT %2d:",p);

                printf("To few elements in file(<n)\n");
                delete[] array;
                delete[] a;
                return 3;
            }
            break;
        }
        else if (res== 0)
        {
            printf("RESULT %2d:",p);

            printf("Problem with reading an element \n");
            
            delete[] array;
            delete[] a;
            return 2;
        }

        if (i >= n)
        {
            printf("RESULT %2d:",p);

            printf("To many elements in file(>n)\n");
            delete[] array;
            delete[] a;
            return 4;
        }
        else
        {
            array[i] = el;
        }
        i++;
    }
    //Обработка массива. Распределение памяти.
    for (k = 0; k < p; k++,pa += m)
    {
        a[k].k = k;
        a[k].p = p;
        a[k].n = n;
        a[k].l = l;
        a[k].name = argv[3];
        a[k].m = (k < p-1 ? m : m+l);
        if (k > 0)
        {
            a[k].pan = *(pa-1);
        }
        if (k < p-1)
        {
            a[k].na0 = *(pa + m);
            if (m > 1 || (((((k < p-1) && (l > 0)) || (k < p-2)) && (m==1))))
            {
                a[k].na1 = *(pa + m + 1);
            }     
        }
     
        a[k].array = pa;
    }

    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            
            delete[] array;
            return 1;
        }
    }
    a[0].tid = pthread_self();

    thread_func(a+0);

    //Прибиваем потоки молотком.
    for (k = 1; k < p; k++)
    {
        pthread_join(a[k].tid, nullptr);
    }
    


    if (a[0].res == 0)
    {

        for (k = 0; k < p; k++)
        {
            a[k].PrintAll(true);
        }

        printf("RESULT %2d:",p);
        for (int i = 0; i < n; i++)
        {
            printf(" %8.2e", array[i]);
        }
        printf("\n");
    }

    else
    {
        for (k = 0; k < p; k++)
        {
            a[k].PrintAll(false);
        }

        if (a[0].reading_file == io_status::error_read)
        {
            return 2;
        }
        else if (a[0].reading_file == io_status::error_open)
        {
            return 1;
        }
        else  if (a[0].reading_file == io_status::error_few_el)
        {
            return 3;
        }
        else if (a[0].reading_file == io_status::error_few_el)
        {
            return 4;
        }    
    } 
    
    delete[] array;
    delete[] a;
    return 0;
}