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
    for (k = 0; k < p; k++,pa += m + (k <= (p-l) ? 0 : 1))
    {
        a[k].k = k;
        a[k].p = p;
        a[k].n = n;
        a[k].l = l;
        a[k].name = argv[3];
        a[k].m = (k < p-l ? m : m+1);
        if (k > 0)
        {
            a[k].prev = *(pa-1);
        }
        if (k < p-l)
        {
            a[k].next = *(pa + m);
        }
        else if (k < p-1)
        {
            a[k].next = *(pa + m + 1);
        }

        a[k].array = pa;
        a[k].left_sum = *pa;
        a[k].right_sum = pa[a[k].m-1];

    }

    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            
            delete[] array;
            delete[] a;
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

    
    double res = ProccesResults(a);
    for (int i = 0; i < p; i++)
    {
        a[i].amount_of_changed = res;
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
            delete[] array;
            delete[] a;
            return 2;
        }
        else if (a[0].reading_file == io_status::error_open)
        {
            delete[] array;
            delete[] a;
            return 1;
        }
        else  if (a[0].reading_file == io_status::error_few_el)
        {
            delete[] array;
            delete[] a;
            return 3;
        }
        else if (a[0].reading_file == io_status::error_few_el)
        {
            delete[] array;
            delete[] a;
            return 4;
        }    
    } 
    
    delete[] array;
    delete[] a;
    return 0;
}