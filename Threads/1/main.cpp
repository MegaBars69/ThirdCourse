#include <iostream>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include "func.hpp"
/*
struct Thread{
    std::string filename;
    int thread_id;
    int total_threads;
    int count_max;
};

// Функция для работы потоков
void* count_max_elements(void* arg) {
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
}*/

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