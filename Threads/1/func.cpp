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
