#ifndef HEADER1  
#define HEADER1
#include <iostream>
#include <sys/sysinfo.h>
#include "initialize_matrix.hpp"

template <typename T>
void reduce_sum(int p,T* a = nullptr, int n = 0)
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
}

void PrintMatrix(double* A,  int n, int m, int r, int p = 0, int K = 0, bool full = true, bool exp_format = false, bool okruglenie = false);

enum class io_status
{
    working,
    not_working,
    succes
};

class Args{
    public:
        int p = 0;
        int k = 0;
        pthread_t tid = 0;
        double cpu_time = 0;
        double cpu_time_of_all_threads = 0;
        double astr_time = 0;
        io_status status = io_status::not_working;
        double res = 0;

        double* A = nullptr;
        double* B = nullptr; 
        double* U = nullptr;
        double* ProductResult = nullptr;
        double* ZeroMatrix = nullptr;

        double norm = 0;

        std::string name = "";

        int n = 0; //Matrix dim
        int m = 0; //Block dim
        int M = 0; //Amount of block strings
        int r = 0; //Amount of printing elements
        int s = 0;

        int cur_str = 0;
        int nomer_v_okne = 0;
        

        void PrintAll() const
        {
            printf("Number of thread: %d \n", k);
            printf("n: %d \n", n);            
            printf("p: %d \n", p);
            printf("m: %d \n", m);
            printf("M: %d \n", M);
            printf("r: %d \n", r);
            printf("s: %d \n", s);
            printf("cur_str: %d \n", cur_str);
            printf("cpu time: %lf \n", cpu_time);
            
            printf("res: %lf \n", res);

          std::cout<<"norm: "<<norm<<std::endl;

            printf("MATRIX A \n");
            PrintMatrix(A, n, m, r, p, k, false);

            printf("MATRIX B \n");
            PrintMatrix(B, n, m, r, p, k, false);

            printf("MATRIX U \n");
            PrintMatrix(U, m, 1, m);
        }
        //Args();
};

void FirstStep(double* A, double* B, double* U, double* ProductResult, double* ZeroMatrix, double norm, int n, int m, int p, int K, int s, Args *a);
void SecondStep(double* A, double* B, double* U, double* ProductResult, double* ZeroMatrix, double norm, int n, int m, int p, int K, int s, Args* aA);

void InverseMatrixParallel(double* A, double* B, double* U, double* ProductResult, double* ZeroMatrix, double norm, int n, int m, int p, int K, int M);

double get_cpu_time();
double get_fun_time();
void* thread_func(void *arg);
void ApplyVector(double* X, double* A, int row_num, int col_num, int k, bool inside = true);
int InverseTriungleBlock(double* A, double* B, int n, double norm);
void ReplaceWith(double*A, double*B, int row_size, int col_size);
void MatrixMinusEqual(double* A, double* B, int row_size, int col_size);
void MatrixPlusEqual(double* A, double* B, int row_size, int col_size);
void DifferenceOfBlockMatrix(double* A, double* B, int n, int m);
void MinusEqualBlockMulEqualBlockMul(double* Main, double *a, double* b, int n1, int m12, int n2);

int Triungulize(double* A, double* U, int row_num, int col_num, double norm);
int InverseMatrix(double* A, double* B, double* U, double* ProductResult, double* ZeroMatrix, double norm, int n, int m, int S);
void ApplyMatrix(double* U, double* A, int row_num, int col_num, int amount_of_vectors);
void ZeroOut(double* Diag, double* Down, double* U, int m, int row_size, double norm);
void ApplyMatrixToPair(double* U, double* Up, double* Down, int col_size, int row_size, int amount_of_vectors, bool down_is_zero = false);

double Norm(double* A, double* results, int n, int m);
void AddToMatrix(double* A, double* B, int n, int m);
void BlockMul(double *a, double* b, double* c, int n1, int m12, int n2);
void BlockMulOptimized(double *a, double* b, double* c, int n1, int m12, int n2);
double CalcError(double * A, double * InversedA, int n, int m);
double Discrepancy(double * A, double * InversedA, double* Column, double* ProductResult, double* Sum, int n, int m);

#endif 