#ifndef HEADER1  
#define HEADER1

#include "initialize_matrix.hpp"


class Args
{
    public:
        int p = 0;
        int k = 0;
        MPI_Comm comm = MPI_COMM_WORLD;
        double cpu_time = 0;
        double cpu_time_of_all_threads = 0;
        double astr_time = 0.0;
        double res = 0;

        double* A = nullptr;
        double* B = nullptr; 
        double* buf = nullptr;
        bool* ZerosMatrix = nullptr;
        double* U = nullptr;
        double* ProductResult = nullptr;
        double* ZeroMatrix = nullptr;
        double* results = nullptr;

        double norm = 0;
        bool last_line_isnt_fool = false;

        std::string name = "";

        int n = 0; //Matrix dim
        int m = 0; //Block dim
        int M = 0; //Amount of block strings
        int r = 0; //Amount of printing elements
        int s = 0;
        int l = 0;
        int K = 0;
        int rows = 0;
        int max_rows = 0;
        int active_procceses = 0; 

        int cur_str = 0;
        int cur_global_str = 0;
        int nomer_v_okne = 0;
        int shag = 0;
        int send_size = 0;

        void PrintInfo() 
        {
            printf("Number of proces: %d \n", k);
            printf("SEND SIZE: %d \n", send_size);            
            printf("SHAG: %d \n", shag);            

            printf("n: %d \n", n);            
            printf("p: %d \n", p);
            printf("m: %d \n", m);
            //printf("r: %d \n", r);
            //printf("s: %d \n", s);
            printf("l: %d \n", l);
            //printf("K: %d \n", K);
            //printf("Norm: %lf \n", norm);
            printf("nomer_v_okne: %d \n", nomer_v_okne);
            //printf("cur_str: %d \n", cur_str);
            //printf("cpu time: %lf \n", cpu_time);
            //printf("cpu time of all_threads: %lf \n", cpu_time_of_all_threads);
           // printf("astr time: %lf \n",astr_time);

            if (!last_line_isnt_fool)
            {
                printf("Last row is fine\n");
            }
            else
            {
                printf("Last row is bad\n");
            }
            
            //PrintLocalMatrix(A, n, m, p ,k, r);
            /*
            for(int bi = 0; bi< max_rows;bi++)
            {
                for (int bj = 0; bj < K+1; bj++)
                {
                    std::cout<<ZerosMatrix[bi*(K+1) + bj]<<" ";
                }
                std::cout<<std::endl;
            }*/
            //printf("res: %lf \n", res);
            /*std::cout<<"BUFFER"<<": ";
            for(cur_global_str = 0; cur_global_str < (2*n + m)*m; cur_global_str++)
            {
                std::cout<<buf[cur_global_str]<<" ";
            }
            std::cout<<std::endl;
            */

        }
};


double get_cpu_time();
double get_full_time();
int Triungulize(double* A, double* U, int row_num, int col_num, double norm);
int InverseMatrix(double* A, double* B, double* U, double* ProductResult, double* ZeroMatrix, double norm, int n, int m, int S);
void ApplyMatrix(double* U, double* A, int row_num, int col_num, int amount_of_vectors);
void ZeroOut(double* Diag, double* Down, double* U, int m, int row_size, double norm, bool down_is_triungle = false);
void ApplyMatrixToPair(double* U, double* Up, double* Down, int col_size, int row_size, int amount_of_vectors, bool down_is_zero = false, bool down_is_triungle = false);
void ApplyMatrixToPairPrev(double* U, double* Up, double* Down, int col_size, int row_size, int amount_of_vectors, bool down_is_zero = false, bool down_is_triungle = false);
void BlockMul(double *a, double* b, double* c, int n1, int m12, int n2);
double Norm(Args* a);
void FirstStep(double* A, double* B, double* U, double norm, int n, int m, int p, int K, int shag, Args *a);
double Discrepancy(double* A, double* B, Args* a) ;
void PrintZeros(Args *a);

int InverseMatrixParallel(Args* a);
#endif