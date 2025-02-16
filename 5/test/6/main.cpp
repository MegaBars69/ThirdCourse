#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <pthread.h>
#include <cmath>
#include "func.hpp"

int main(int argc, char* argv[])
{
    Args* a;
    double *A;
    double *pa = nullptr;
    int k = 0;

    if(argc < 5)
    {
        printf("RESULT :");

        printf("Usage ./a.out p n1 n2 filename\n");
        return 1;
    }
    int n1 = atoi(argv[2]);
    int n2 = atoi(argv[3]);
    int p = atoi(argv[1]); 

    if (p > n1)
    {
        p = n1; 
    }
    if (p <= 0)
    {
        printf("RESULT %2d:",p);

        printf(" p > 0 \n");
        return 4;
    }
    
    int m = n1/p;
    int l = n1%p;

    A = new double[n1*n2];     
    a = new Args[p];

    pa = A;
    double el;
    int i = 0;

    FILE *fp = fopen(argv[4], "r");
    if(fp == nullptr)
    {
        printf("Can't open file \n");
        delete[] A;
        delete[] a;

        return 1;
    }
    while (true) 
    {
        int res = fscanf(fp, "%lf", &el);
        if(res == EOF)
        {  
            if(i != n1*n2)
            {
                printf("To few elements in file(<n)\n");
                delete[] A;
                delete[] a;
                fclose(fp);
                return 3;
            }
            break;
        }
        else if (res== 0)
        {
            printf("Problem with reading an element \n");
            
            delete[] A;
            delete[] a;
            fclose(fp);
            return 2;
        }

        if (i >= n1*n2)
        {
            printf("To many elements in file(>n)\n");
            delete[] A;
            delete[] a;
            fclose(fp);
            return 4;
        }
        else
        {
            A[i] = el;
        }
        i++;
    }
        
    //Обработка массива. Распределение памяти.
    for (k = 0; k < p; k++)
    {
        a[k].k = k;
        a[k].p = p;
        a[k].n1 = n1;
        a[k].n2 = n2;
        a[k].l = l;
        a[k].m = (k < p-l ? m : m+1);
        a[k].m1 = a[k].m;
        a[k].m2 = n2;

        a[k].prev_line = new double[n2];
        a[k].next_line = new double[n2];
        if (k == 0 || k == p-1)
        {
            memset(a[k].prev_line, 0, n2*sizeof(double));
            memset(a[k].next_line, 0, n2*sizeof(double));  
        }
        
        if (k > 0)
        {
            memcpy(a[k].prev_line, pa - n2, n2*sizeof(double));  
        }
        
        a[k].A = pa;
        pa += n2*a[k].m;

        if (k < p-1)
        {
            memcpy(a[k].next_line, pa, n2*sizeof(double));  
        }
    }
    PrintMatrix(A, n1,n2);
    /*
    for (k = 0; k < p; k++) 
    {
        a[k].PrintAll(true);
    }
    */
    //Запуск потоков.
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            
            for (int i = 0; i < p; i++) 
            {
                //a[k].PrintAll(false);
                delete[] a[i].prev_line;
                delete[] a[i].next_line;
            }
            delete[] A;
            delete[] a;
            fclose(fp);

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
    
    for (k = 0; k < p; k++) 
    {
        //a[k].PrintAll(false);
        delete[] a[k].prev_line;
        delete[] a[k].next_line;
    }

    a[0].PrintAll(false);

    /*if (a[0].res > 0)
    {
        printf("Avarage doesn't exist\n");
        delete[] A;
        delete[] a;
        return 3;
    }*/
    
    PrintMatrix(A, n1, n2);
    
    
    
    delete[] A;
    delete[] a;
    fclose(fp);
    return 0;
}