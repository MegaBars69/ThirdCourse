#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"
/*
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file1> <file2> ... <fileN>" << std::endl;
        return 1;
    }

    int total_threads = argc - 1;
    pthread_t* threads = new pthread_t[total_threads];
    Thread* thread_data = new Thread[total_threads];
    int* results = new int[total_threads];

    // Создаем потоки
    for (int i = 0; i < total_threads; ++i) {
        thread_data[i].filename = argv[i + 1];
        thread_data[i].thread_id = i + 1;
        thread_data[i].total_threads = total_threads;
        thread_data[i].count_max = 0;

        int ret = pthread_create(&threads[i], nullptr, count_max_elements, &thread_data[i]);
        if (ret) {
            std::cerr << "Error creating thread " << i << ": " << ret << std::endl;
            return 1;
        }
    }

    // Ждем завершения потоков и собираем результаты
    int total_max_count = 0;
    for (int i = 0; i < total_threads; ++i) {
        void* retval;
        pthread_join(threads[i], &retval);

        if (retval == reinterpret_cast<void*>(-1)) {
            std::cerr << "Error opening file: " << thread_data[i].filename << std::endl;
            delete[] threads;
            delete[] results;
            delete[] thread_data;
            return 1;
        } else if (retval == reinterpret_cast<void*>(-2)) {
            std::cerr << "Error reading element in file: " << thread_data[i].filename << std::endl;
            delete[] threads;
            delete[] results;
            delete[] thread_data;
            return 1;
        }

        total_max_count += thread_data[i].count_max;
    }

    std::cout << "Total number of maximum elements: " << total_max_count << std::endl;

    delete[] threads;
    delete[] results;
    delete[] thread_data;
    return 0;
}
*/
int main(int argc, char* argv[])
{
    Args* a;
    if(argc == 1)
    {
        printf("Usage %s files ", argv[0]);
        return 1;
    }

    int p = argc - 1;
    a = new Args[p];

    int k;
    for (k = 0; k < p; k++)
    {
        a[k].k = k;
        a[k].p = p;
        a[k].filename = argv[k+1];
    }
    for (k = 1; k < p; k++)
    {
        if (pthread_create(&a[k].tid, nullptr, thread_func, a+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            return 1;
        }
    }
    a[0].tid = pthread_self();

    thread_func(a+0);
    //reduce_sum(p);

    for (k = 1; k < p; k++)
    {
        pthread_join(a[k].tid, nullptr);
    }
    
    if (procces_args(a) == 0)
    {
        printf("Result = %lf \n", a[0].res);
    }
    
    delete[] a;
    return 0;
}