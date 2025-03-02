#include "mpi.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <string.h>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <unistd.h> 
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include "algorithm.hpp"
#define EPSILON pow(10,-16)

using namespace std;

/*
double get_fun_time()
{
    struct timeval buf;
    gettimeofday(&buf,0);
    return buf.tv_sec + buf.tv_sec/1e6;
}*/

double get_full_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return (double)(buf.tv_sec) + (double)(buf.tv_usec)/1000000.;
}

double get_cpu_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);

    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec/1e6;
}

/*
double Discrepancy(Args* a)
{
    double * A = a->A, InversedA = a->B, Column = a->results, ProductResult = a->ProductResult, Sum = a->ZeroMatrix;
    int n = a->n, m = a->m;
    int l = n%m;
    int k = (n-l)/m;
    int block_size_row, block_size_col, m12;
    double final_answer = 0, sum = 0;

    for (int j = 0; j < m; j++)
    {
        Column[j] = 0;
    }
    
    for (int bj = 0; bj < k+1; bj++)
    {
        block_size_col = (bj < k ? m : l);
        for (int bi = 0; bi < k+1; bi++)
        {
            block_size_row = (bi < k ? m : l);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    Sum[i*m + j] = 0;
                    ProductResult[i*m +j] = 0;
                }  
            }       
            for (int s = 0; s < k+1; s++)
            {
                m12 = (s < k ? m : l);
                BlockMul((A + bi*m*n + s*m*block_size_row), (InversedA + s*m*n + bj*m*m12), ProductResult,block_size_row, m12, block_size_col);
                MatrixPlusEqual(Sum, ProductResult, block_size_row, block_size_col);
            }
            if (bi == bj)
            {
                for (int i = 0; i < block_size_row; i++)
                {
                    Sum[i*block_size_row + i] -= 1;
                }
                
            }
            
            for (int j = 0; j < block_size_col; j++)
            {
                sum = 0;
                for (int i = 0; i < block_size_row; i++)
                {
                    sum += fabs(Sum[i*m + j]);
                }
                Column[j] += sum;
            }
        }
        for (int j = 0; j < block_size_col; j++)
        { 
            double el = Column[j];
            final_answer = (final_answer > el ? final_answer : el);
        }
        for (int j = 0; j < block_size_col; j++)
        {
            Column[j] = 0;
        }  
    }

    return final_answer; 
}
*/
void block_mult_add(double* A, int av, int ag, double* B, int bg, double* C) {
    int av_reminder = av % 3;
    int bg_reminder = bg % 3;
    int _av = av - av_reminder;
    int _bg = bg - bg_reminder;
    int i = 0, j = 0, q = 0;
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    for(i = 0; i < _av; i += 3) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            c10 = 0; c11 = 0; c12 = 0;
            c20 = 0; c21 = 0; c22 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c11 += A[ag * (i + 1) + q] * B[bg * q + (j + 1)];
                c12 += A[ag * (i + 1) + q] * B[bg * q + (j + 2)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
                c21 += A[ag * (i + 2) + q] * B[bg * q + (j + 1)];
                c22 += A[ag * (i + 2) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i) + (j + 1)] += c01;
            C[bg * (i) + (j + 2)] +=c02;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 1) + (j + 1)] += c11;
            C[bg * (i + 1) + (j + 2)] += c12;
            C[bg * (i + 2) + (j)] += c20;
            C[bg * (i + 2) + (j + 1)] += c21;
            C[bg * (i + 2) + (j + 2)] += c22;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            c10 = 0;
            c20 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c10 += A[ag * (i + 1) + q] * B[bg * q + (j)];
                c20 += A[ag * (i + 2) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i + 1) + (j)] += c10;
            C[bg * (i + 2) + (j)] += c20;
        }
    }
    for(; i < av; ++i) {
        for(j = 0; j < _bg; j+= 3) {
            c00 = 0; c01 = 0; c02 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
                c01 += A[ag * (i) + q] * B[bg * q + (j + 1)];
                c02 += A[ag * (i) + q] * B[bg * q + (j + 2)];
            }
            C[bg * (i) + (j)] += c00;
            C[bg * (i) + (j + 1)] += c01;
            C[bg * (i) + (j + 2)] += c02;
        }
        for(;j < bg; ++j) {
            c00 = 0;
            for(q = 0; q < ag; ++q) {
                c00 += A[ag * (i) + q] * B[bg * q + (j)];
            }
            C[bg * (i) + (j)] += c00;
        }
    }
}
double calculate_discrepancy(double* matrix, double* inverse, int n, int m, double* tmp_block_m
        , double* tmp_line_n, int proc_num, int p, MPI_Comm comm, double* inverse_buf) 
{
    int k = n / m;
    int l = n - m * k;
    int el_in_blcok_line = k * m * m + l * m;
    int i_glob = 0, j = 0, q = 0;
    bool small_row = (l != 0 && k%p == proc_num) ? true : false;
    double sum = 0;
    int max_rows = get_max_rows(n, m, p);
    int rows = get_rows(n, m, p, proc_num);
    rows = (small_row == true) ? rows - 1 : rows;
    memset(tmp_line_n, 0, n * sizeof(double));
    for(i_glob = 0; i_glob < k; ++i_glob) {
        for(j = 0; j < rows; ++j) {
            memcpy(inverse_buf + j * m * m, inverse + j * el_in_blcok_line + i_glob * m * m, m * m * sizeof(double));
        }
        if (small_row == true) {
            memcpy(inverse_buf  + j * m * m, inverse + j * el_in_blcok_line + i_glob * m * l, m * l * sizeof(double));
        }
        MPI_Allgather(inverse_buf, max_rows * m * m, MPI_DOUBLE, matrix + max_rows * el_in_blcok_line
                , max_rows * m * m, MPI_DOUBLE, comm);
        for(q = 0; q < rows; ++q) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                        , m, tmp_block_m);
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                    , m, tmp_block_m); 
            if ((q * p + proc_num) == i_glob) {
                for(int s = 0; s < m; ++s) {
                    tmp_block_m[s * m + s] -= 1;
                }
            }
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
        if (small_row == true) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                        , m, tmp_block_m);
        
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * m * m + j_loc * m * m
                    , m, tmp_block_m); 
            for(int w = 0; w < m; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * m + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
    }

    if (l != 0) {
        for(j = 0; j < rows; ++j) {
            memcpy(inverse_buf + j * l * m, inverse + j * el_in_blcok_line + i_glob * m * m, l * m * sizeof(double));
        }
        if (small_row == true) {
            memcpy(inverse_buf  + j * l * m, inverse + j * el_in_blcok_line + i_glob * m * l, l * l * sizeof(double));
        }

        MPI_Allgather(inverse_buf, max_rows * m * l, MPI_DOUBLE, matrix + max_rows * el_in_blcok_line
                , max_rows * m * l, MPI_DOUBLE, comm);
        for(q = 0; q < rows; ++q) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * m * l + j_loc * l * m
                        , l, tmp_block_m);
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * m, m, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * l * m + j_loc * l * m
                    , l, tmp_block_m); 
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < m; ++v) {
                    sum += fabs(tmp_block_m[v * l + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
        if (small_row == true) {
            memset(tmp_block_m, 0, m * m * sizeof(double));
            for(j = 0; j < k; ++j) {
                int proc_j = j % p;
                int j_loc = j / p; 
                block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, m
                        , matrix  + max_rows * el_in_blcok_line + proc_j * max_rows * l * m + j_loc * l * m
                        , l, tmp_block_m);
        
            }
            int proc_j = j % p;
            int j_loc = j / p; 
            block_mult_add(matrix + q * el_in_blcok_line + j * m * l, l, l
                    , matrix + max_rows * el_in_blcok_line + proc_j * max_rows * l * m + j_loc * l * m
                    , l, tmp_block_m); 
            if ((q * p + proc_num) == i_glob) {
                for(int s = 0; s < l; ++s) {
                    tmp_block_m[s * l + s] -= 1;
                }
            }
            for(int w = 0; w < l; ++w) {
                sum = 0;
                for(int v = 0; v < l; ++v) {
                    sum += fabs(tmp_block_m[v * l + w]);
                }
                tmp_line_n[i_glob * m + w] += sum;
            }
        }
    }
    MPI_Allreduce(tmp_line_n, inverse_buf, n, MPI_DOUBLE, MPI_SUM, comm);

    double max = inverse_buf[0];
    for(int i = 0; i < n; ++i) {
        if(inverse_buf[i] > max) max = inverse_buf[i];
    }
    return max;
}

double Norm(Args *a)
{
    int n = a->n;
    int l = a->l;
    int m = a->m;
    int k = a->K;
    int p = a->p;
    int send_proc_num = (a->k + 1 + p)%p;
    int recv_proc_num = (a->k - 1 + p)%p;
    bool last_line_isnt_fool = a->last_line_isnt_fool;
    int block_size_col, block_size_row = m;
    double sum_in_col = 0, final_norm = 0, el;
    int rows = a->rows;
    double *A = a->A;
    double* results = a->results;
    double* buf = a->buf;
    double* result_copy = (buf + n);
    MPI_Comm comm = a->comm;

    for (int bj = 0; bj < k+1; bj++)
    {
        block_size_col = (bj < k ? m : l);
        for (int bi = 0; bi < rows ; bi++)
        {
            if(bi == rows-1)
            {
                if(last_line_isnt_fool)
                {
                    block_size_row = l;
                } 
                else
                {
                    block_size_row = m;
                }
            }

            double* pa = (A + m*n*bi + m*block_size_row*bj);
            for (int j = 0; j < block_size_col; j++)
            {
                sum_in_col = 0;
                for (int i = 0; i < block_size_row; i++)
                {
                    sum_in_col += fabs(pa[i*block_size_col + j]);
                }
                results[bj*m + j] += sum_in_col;
                result_copy[bj*m + j] += sum_in_col;
            }
        }
    }   
    for (int i = 1; i < p && rows > 0; i++)
    {   
        if(a->k != recv_proc_num)
        {
            MPI_Status st;
            //cout<<"Send: "<<send_proc_num<<"  Recv:"<<recv_proc_num<<endl;
            MPI_Sendrecv( result_copy, n, MPI_DOUBLE, send_proc_num,0, buf, n, MPI_DOUBLE, recv_proc_num, 0, comm, &st);

            //MPI_Recv(buf, n, MPI_DOUBLE, recv_proc_num, 0, comm, &st);
            
            for (int j = 0; j < n; j++)
            {
                results[j] += buf[j];
            }
            
            send_proc_num = (send_proc_num + 1 + p)%p;  
            recv_proc_num = (recv_proc_num - 1 + p)%p;          
        }
    }

    for (int j = 0; j < n; j++)
    {
        el = results[j];
        results[j] = 0;
        final_norm = (el < final_norm ? final_norm : el);
    }   
    return final_norm;
}

void MinusEqualBlockMul(double* Main, double *a, double* b, int n1, int m12, int n2)
{
    int l_col = n2%3;
    int l_row= n1%3;
    int k_col = (n2-l_col)/3;
    int k_row = (n1-l_row)/3; 
    int bi, bj, r;
    
    double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;

    for ( bi = 0; bi < k_row; bi++)
    {
        for ( bj = 0; bj < k_col; bj++)
        {
            s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;
            for ( r = 0; r < m12; r++)
            {
                s00 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + bj*3];
                s01 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 1];
                s02 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 2];
                s10 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + bj*3];
                s11 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 1];
                s12 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 2]; 
                s20 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + bj*3];
                s21 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + bj*3 + 1];
                s22 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + bj*3 + 2];
            }
            Main[bi*3*n2 + bj*3] -= s00; Main[bi*3*n2 + bj*3 + 1] -= s01; Main[bi*3*n2 + bj*3 + 2] -= s02;
            Main[(bi*3 + 1)*n2 + bj*3] -= s10; Main[(bi*3 + 1)*n2 + bj*3 + 1] -= s11; Main[(bi*3 + 1)*n2 + bj*3 + 2] -= s12;
            Main[(bi*3 + 2)*n2 + bj*3] -= s20; Main[(bi*3 + 2)*n2 + bj*3 + 1] -= s21; Main[(bi*3 + 2)*n2 + bj*3 + 2] -= s22;
        }
        if (l_col != 0)
        {

            s00 = 0, s01 = 0, s10 = 0, s11 = 0, s20 = 0, s21 = 0;
            for ( r = 0; r < m12; r++)
            {
                if(l_col > 1)
                { 
                    s01 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + k_col*3 + 1];
                    s11 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + k_col*3 + 1];
                    s21 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + k_col*3 + 1];

                }
                s00 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + k_col*3];
                s10 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + k_col*3];
                s20 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + k_col*3];
            }

            Main[bi*3*n2 + k_col*3] -= s00; 
            Main[(bi*3 + 1)*n2 + k_col*3] -= s10;
            Main[(bi*3 + 2)*n2 + k_col*3] -= s20;
            if(l_col > 1)
            { 
                Main[bi*3*n2 + k_col*3 + 1] -= s01; 
                Main[(bi*3 + 1)*n2 + k_col*3 + 1] -= s11;
                Main[(bi*3 + 2)*n2 + k_col*3 + 1] -= s21;
            }
        }
            
    }

    if(l_row != 0)
    {
        for ( bj = 0; bj < k_col; bj++)
        {
            s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0;
            for ( r = 0; r < m12; r++)
            {
                
                s00 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + bj*3];
                s01 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 1];
                s02 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 2]; 
                if(l_row > 1)
                {
                    s10 += a[(k_row*3 + 1)*m12 + r] * b[r*n2 + bj*3];
                    s11 += a[(k_row*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 1];
                    s12 += a[(k_row*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 2]; 
                }
            }

            Main[(k_row*3)*n2 + bj*3] -= s00; Main[(k_row*3 )*n2 + bj*3 + 1] -= s01; Main[(k_row*3 )*n2 + bj*3 + 2] -= s02;
            if (l_row > 1)
            {
                Main[(k_row*3  + 1)*n2 + bj*3] -= s10; Main[(k_row*3 + 1)*n2 + bj*3 + 1] -= s11; Main[(k_row*3 + 1)*n2 + bj*3 + 2] -= s12;
            }
            
        }
        if(l_col != 0)
        {
            s00 = 0, s01 = 0, s10 = 0, s11 = 0;
            for ( r = 0; r < m12; r++)
            {
                
                s00 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + k_col*3 ];
                if (l_col > 1)
                {
                    s01 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + k_col*3 + 1];
                }
                if(l_row > 1)
                {
                    s10 += a[(k_row*3  + 1)*m12 + r] * b[r*n2 + k_col*3 ];
                }
                if (l_col > 1 && l_row > 1)
                {
                    s11 += a[(k_row*3  + 1)*m12 + r] * b[r*n2 + k_col*3 + 1];
                }
                
            }

            Main[(k_row*3)*n2 + k_col*3] -= s00;
            if (l_col>1)
            {
                Main[(k_row*3)*n2 + k_col*3 + 1] -= s01;
            }
            
            if (l_row > 1)
            {
                Main[(k_row*3 + 1)*n2 + k_col*3] -= s10;
            }
            if (l_row > 1 && l_col>1)
            {
                Main[(k_row*3 + 1)*n2 + k_col*3 + 1] -= s11;
            }
            
        }     
    }       
}

double fabs(double a)
{
    return (a>0 ? a:(-a));
}
void ApplyVector(double* X, double* A, int row_num, int col_num, int k, bool inside)
{
    double s;
    double* px;
    int i, j;
    for (j = (inside ? k+1: 0); j < col_num; j++)
    {
        s = 0;
        px = X;
        for (i = k; i < row_num; i++,px++)
        {
            s += (*px) * (A[i*col_num + j]);
        }

        px = X;
        s *= 2;
        for (i = k; i < row_num; i++,px++)
        {
            A[i*col_num + j]-= s*(*px);
        }        
    }   
}

void ApplyMatrix(double* U, double* A, int row_num, int col_num, int amount_of_vectors)
{
    //PrintMatrix(U, row_num, 1, row_num);
    int k;
    for (k = 0; k < amount_of_vectors; k++)
    {
        ApplyVector((U + k*row_num), A, row_num, col_num, k, false);
    }
}

int Triungulize(double* A, double* U, int row_num, int col_num, double norm)
{
    double sk, ajk,akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;
    int i,j,k;

    for (k = 0; k < col_num; k++)
    {
        //Finding vector xk
        pu = (U + k*row_num);
        sk = 0;
        for (j = k+1; j < row_num; j++)
        {
            ajk = A[j*col_num + k];
            pu++;
            *pu = ajk;
            sk+= ajk*ajk;
            A[j*col_num + k] = 0;
        }
        akk = A[k*col_num + k];
        new_diag_el = sqrt(sk + akk*akk);
        
        A[k*col_num + k] = new_diag_el;

        first_in_x = akk - new_diag_el;
        
        pu = (U + k*row_num);
        *pu = first_in_x;
        
        norm_xk = sqrt(sk + first_in_x*first_in_x);

        if(!(fabs(norm_xk) < norm))
        {
            for (i = k; i < row_num; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
            }
        }

        //Applying Reflection
        
        ApplyVector((U + k*row_num), A, row_num, col_num, k, true); 
    }

    return 0;
}


void ZeroOut(double* Diag, double* Down, double* U, int m, int row_size, double norm, bool down_is_triungle)
{
    double sk, ajk, akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;
    double s;
    double* px;
    int i, j,jj;
    int up_bound = (!down_is_triungle ? row_size : 1);
    for (j = 0; j < m; j++)
    {
        sk = 0;
        pu = (U + j*(row_size+1));
        
        for (i = 0; i < up_bound; i++)
        {
            ajk = Down[i*m + j];
            pu++;
            *pu = ajk;
            sk+= ajk*ajk;
        }

        akk = Diag[j*m + j];
        new_diag_el = sqrt(sk + akk*akk);
        first_in_x = akk - new_diag_el;
        
        pu = (U + j*(row_size+1));
        
        *pu = first_in_x;
        norm_xk = sqrt(sk + first_in_x*first_in_x);
        
        if(!(fabs(norm_xk) < norm))
        {
            *pu = (*pu)/norm_xk;
            pu++;
            for (i = 0; i < up_bound; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
                Down[i*m + j] = 0;
            }
        
            //Applying Reflection
            
            Diag[j*m + j] = new_diag_el;
                
            for (jj = j+1; jj < m; jj++)
            {
                s = 0;
                px = U + j*(row_size+1);
                s += (*px) * (Diag[j*m + jj]);
                px++;
                for (i = 0; i < up_bound; i++, px++)
                {
                    s += (*px) * (Down[i*m + jj]);
                }
                px = U + j*(row_size+1);

                Diag[j*m + jj]-= 2*s*(*px);
                px++;
                s*=2;
                
                for (i = 0; i < up_bound; i++, px++)
                {
                    Down[i*m + jj]-= s*(*px);
                }
                
            }
        }
        if (down_is_triungle)
        {
            up_bound += (row_size > up_bound ? 1 : 0);
        }
    }
    
}

void ApplyMatrixToPair(double* U, double* Up, double* Down, int col_size, int row_size, int amount_of_vectors, bool down_is_zero, bool down_is_triungle)
{
    double s;
    double* pu;
    int vec_num = 0, j, i;
    int up_bound = (!down_is_triungle ? row_size : 1);

    for (j = 0; j < col_size; j++)
    {
        pu = U + vec_num*(row_size + 1);
        
        s = 0;

        s += (*pu) * (Up[vec_num*col_size + j]);
        if (!down_is_zero)
        {   
            for (i = 0; i < up_bound; i++)
            {
                pu++;
                s += (*pu) * (Down[i*col_size + j]);
            }
        }

        s*=2;

        pu = U + vec_num*(row_size + 1);
        Up[vec_num*col_size + j] -= s*(*pu);
        for (i = 0; i < up_bound; i++)
        {
            pu++;
            Down[i*col_size + j] -= s*(*pu);
        }
    }           
    for (vec_num = 1; vec_num < amount_of_vectors; vec_num++)
    {   
        if (down_is_triungle)
        {
            up_bound += (row_size > up_bound ? 1 : 0);
        }

        for (j = 0; j < col_size; j++)
        {
            pu = U + vec_num*(row_size + 1);
            
            s = 0;

            s += (*pu) * (Up[vec_num*col_size + j]);
        
            for (i = 0; i < up_bound; i++)
            {
                pu++;
                s += (*pu) * (Down[i*col_size + j]);
            }
        
            s*=2;
            pu = U + vec_num*(row_size + 1);
            Up[vec_num*col_size + j] -= s*(*pu);

            for (i = 0; i < up_bound; i++)
            {
                pu++;
                Down[i*col_size + j] -= s*(*pu);
            }
        }           
    }
    
}

int InverseTriungleBlock(double* A, double* B, int n, double norm)
{
    double* pa = A, *pb = B;
    double diag_el, sum;
    int i,j, s;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i*n+j] = (i != j ? 0 : 1);
        }
        
    }
    
    for (i = n-1; i >= 0; i--)
    {
        pa = A + i*n + i;
        pb = B + i*n;
        diag_el = *(pa);
        if (fabs(diag_el) <= norm)
        {
            return -1;
        }
    
        for (j = 0; j < n; j++)
        {
            if(j > i)
            {
                pa++;
                *pa = (*pa)/diag_el;
                
            }
            *pb = (*pb)/diag_el;
            pb++;
        }
                
        for (j = 0; j < n; j++)
        {
            sum = 0;
            for (s = i+1; s < n; s++)
            {
                sum += A[i*n + s]*B[s*n + j];
            }
            B[i*n+j]-= sum;
        }   
    }
    return 0;
}

double log2(double x)
{
    return log(x)/log(2);
}

void BlockMul(double *a, double* b, double* c, int n1, int m12, int n2)
{
    int l_col = n2%3;
    int l_row= n1%3;
    int k_col = (n2-l_col)/3;
    int k_row = (n1-l_row)/3; 
    int bi, bj, r;
    
    double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;

    for (bi = 0; bi < k_row; bi++)
    {
        for (bj = 0; bj < k_col; bj++)
        {
            s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;
            for (r = 0; r < m12; r++)
            {
                s00 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + bj*3];
                s01 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 1];
                s02 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 2];
                s10 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + bj*3];
                s11 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 1];
                s12 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 2]; 
                s20 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + bj*3];
                s21 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + bj*3 + 1];
                s22 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + bj*3 + 2];
            }
            c[bi*3*n2 + bj*3] = s00; c[bi*3*n2 + bj*3 + 1] = s01; c[bi*3*n2 + bj*3 + 2] = s02;
            c[(bi*3 + 1)*n2 + bj*3] = s10; c[(bi*3 + 1)*n2 + bj*3 + 1] = s11; c[(bi*3 + 1)*n2 + bj*3 + 2] = s12;
            c[(bi*3 + 2)*n2 + bj*3] = s20; c[(bi*3 + 2)*n2 + bj*3 + 1] = s21; c[(bi*3 + 2)*n2 + bj*3 + 2] = s22;
        }
        if (l_col != 0)
        {

            s00 = 0, s01 = 0, s10 = 0, s11 = 0, s20 = 0, s21 = 0;
            for (r = 0; r < m12; r++)
            {
                if(l_col > 1)
                { 
                    s01 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + k_col*3 + 1];
                    s11 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + k_col*3 + 1];
                    s21 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + k_col*3 + 1];

                }
                s00 += a[(bi*3 + 0)*m12 + r] * b[r*n2 + k_col*3];
                s10 += a[(bi*3 + 1)*m12 + r] * b[r*n2 + k_col*3];
                s20 += a[(bi*3 + 2)*m12 + r] * b[r*n2 + k_col*3];
            }

            c[bi*3*n2 + k_col*3] = s00; 
            c[(bi*3 + 1)*n2 + k_col*3] = s10;
            c[(bi*3 + 2)*n2 + k_col*3] = s20;
            if(l_col > 1)
            { 
                c[bi*3*n2 + k_col*3 + 1] = s01; 
                c[(bi*3 + 1)*n2 + k_col*3 + 1] = s11;
                c[(bi*3 + 2)*n2 + k_col*3 + 1] = s21;
            }
        }
            
    }

    if(l_row != 0)
    {
        for (bj = 0; bj < k_col; bj++)
        {
            s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0;
            for (r = 0; r < m12; r++)
            {
                
                s00 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + bj*3];
                s01 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 1];
                s02 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + bj*3 + 2]; 
                if(l_row > 1)
                {
                    s10 += a[(k_row*3 + 1)*m12 + r] * b[r*n2 + bj*3];
                    s11 += a[(k_row*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 1];
                    s12 += a[(k_row*3 + 1)*m12 + r] * b[r*n2 + bj*3 + 2]; 
                }
            }

            c[(k_row*3)*n2 + bj*3] = s00; c[(k_row*3 )*n2 + bj*3 + 1] = s01; c[(k_row*3 )*n2 + bj*3 + 2] = s02;
            if (l_row > 1)
            {
                c[(k_row*3  + 1)*n2 + bj*3] = s10; c[(k_row*3 + 1)*n2 + bj*3 + 1] = s11; c[(k_row*3 + 1)*n2 + bj*3 + 2] = s12;
            }
            
        }
        if(l_col != 0)
        {
            s00 = 0, s01 = 0, s10 = 0, s11 = 0;
            for (int r = 0; r < m12; r++)
            {
                
                s00 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + k_col*3 ];
                if (l_col > 1)
                {
                    s01 += a[(k_row*3 + 0)*m12 + r] * b[r*n2 + k_col*3 + 1];
                }
                if(l_row > 1)
                {
                    s10 += a[(k_row*3  + 1)*m12 + r] * b[r*n2 + k_col*3 ];
                }
                if (l_col > 1 && l_row > 1)
                {
                    s11 += a[(k_row*3  + 1)*m12 + r] * b[r*n2 + k_col*3 + 1];
                }
                
            }

            c[(k_row*3)*n2 + k_col*3] = s00;
            if (l_col>1)
            {
                c[(k_row*3)*n2 + k_col*3 + 1] = s01;
            }
            
            if (l_row > 1)
            {
                c[(k_row*3 + 1)*n2 + k_col*3] = s10;
                
            }
            if (l_row > 1 && l_col>1)
            {
                c[(k_row*3 + 1)*n2 + k_col*3 + 1] = s11;
            }
            
        }     
    }       
}

void ReplaceWith(double*A, double*B, int row_size, int col_size)
{
    for (int i = 0; i < row_size; i++)
    {
        for (int j = 0; j < col_size; j++, A++, B++)
        {
            *A = *B;
        }
        
    }
    
}

bool MatrixIsZero(double* A, int col_size, int row_size, double eps)
{
    int i,j;
    for (i = 0; i < row_size; i++)
    {
        for (j = 0; j < col_size; j++)
        {
            if (fabs(A[i*col_size + j]) > eps)
            {
                return false;
            }
            
        }
        
    }
    return true;
}
bool BothOfMatrixAreZero(double* A, int col_size_A, int row_size_A, double* B, int col_size_B, int row_size_B, double eps)
{
    int i1,j1,i2,j2;
    for (i1 = 0, i2 = 0; i1 < row_size_A || i2 < row_size_B; i1++, i2++)
    {
        for (j1 = 0, j2 = 0; j1 < col_size_A || j1 < col_size_B ; j1++, j2++)
        {
            if (i1 < row_size_A && j1 < col_size_A && fabs(A[i1*col_size_A + j1]) > eps)
            {
                return false;
            }
            if (i2 < row_size_B && j2 < col_size_B && fabs(B[i2*col_size_B + j2]) > eps)
            {
                return false;
            }           
        }
        
    }
    return true;
}

void FirstStep(Args *a)
{
    double* A = a->A, *B = a->B, *U = a->U;
    double norm = a->norm;
    int n = a->n, m = a->m, p = a->p, shag = a->shag;

    int l = a->l;
    int j, bj, bi;
    int k = a->K;
    int K = a->k;
    int rows = a->rows;
    int block_size_row = m, block_size_col, size = a->s, down_block_size_row = m, down_block_size_col;
    double* pa, *pa_side, *pa_down, *pa_down_side, *pb, *pb_down; 
    bool* ZerosMatrix = a->ZerosMatrix;
    int up_bound = rows;
    int rows_up_bound = rows;
    int col_up_bound = (l > 0 ? k+1: k);
    double eps = (a->norm > 10 && size != 4?  EPSILON : a->norm);
    bool last_line_isnt_fool = a->last_line_isnt_fool;

    int s = (a->cur_str);
    int prev_s = K + (a->cur_str)*p;;
    
    int dop_up_bound = (shag == 0 ? K + (a->cur_str)*p+1 : col_up_bound);

    if (s == rows - 1)
    {
        if(last_line_isnt_fool)
        {
            block_size_row = l; 
        }
    }
    
    if(prev_s < k || (prev_s == k && shag == k) )
    {
        pa = A + s*m*n + shag*block_size_row*m;
        
        // First part of algorithm

        Triungulize(pa, U,block_size_row, block_size_row, norm);

        for (j = shag + 1, pa_side = pa + block_size_row*m; j < col_up_bound; j++, pa_side += block_size_row*m)
        {
            block_size_col = (j < k ? m : l);

            ApplyMatrix(U, pa_side, block_size_row, block_size_col, block_size_row);
            
        }
        
        for (j = 0, pb = B + s*m*n; j < dop_up_bound; j++, pb += block_size_row*m)
        {
            block_size_col = (j < k ? m : l);
            
            /*if (ZerosMatrix[s*(k+1) + j])
            {*/
                ApplyMatrix(U, pb, block_size_row, block_size_col, block_size_row);
            //}
            
        
        }

        // Second part of algorithm
        for (bi = s+1; bi < rows_up_bound; bi++)
        {
            if(bi == rows - 1)
            {
                if(last_line_isnt_fool)
                {
                    down_block_size_row = l;
                }
            }
            pa_down = A + bi*m*n + shag*down_block_size_row*m;
            
            ZeroOut(pa, pa_down, U, m, down_block_size_row, norm);

            for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < col_up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
            {
                down_block_size_col = (bj < k ? m : l);
                          
                ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row);
            }
            
            for (bj = 0, pb = B + s*m*n, pb_down = B + bi*m*n; bj < col_up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
            {
                down_block_size_col = (bj < k ? m : l);

                /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                {*/
                    ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row);
                    /*ZerosMatrix[s*(k+1) + bj] = 1;
                    if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                    {
                        ZerosMatrix[bi*(k + 1) + bj] = 1; 
                    }
                    else
                    {
                        ZerosMatrix[bi*(k + 1) + bj] = 0; 
                    }
                    
                }*/
                
            }   
        } 
    }   
}

void SecondStepUpgr(Args* aA)
{
    double* A = aA->A, *B = aA->B, *U = aA->U, *ProductResult = aA->ProductResult;
    double norm = aA->norm;
    int n = aA->n, m = aA->m, p = aA->active_procceses, K = aA->k, shag = aA->shag;
    int a = int(log2(p));
    int two_in_the_power_of_a = pow(2,a);
    int b = p - two_in_the_power_of_a;
    int bi, bj, j;
    int x;
    MPI_Status st;
    MPI_Comm comm = aA->comm;

    int l = n%m;
    int k = (n-l)/m;
    int block_size_row, block_size_col, size = aA->s, down_block_size_row, down_block_size_col;
    double* pa, *pa_side, *pa_down, *pa_down_side, *pb, *pb_down; 
    double* buf_pa = aA->buf, *buf_pb = (aA->buf + n*m);
    double* ZeroMatrix = aA->ZeroMatrix;
    bool* ZerosMatrix = aA->ZerosMatrix;
    
    double eps = (aA->norm > 10 && size != 4 ? EPSILON : aA->norm);
    
    int up_bound = (l > 0 ? k+1: k);

    int s = aA->cur_str;
    int prev_s = K + p*s;

    int step;
    int send_size_a = aA->send_size;
    
    int send_size_b = n*m;

    int Nomer = aA->nomer_v_okne;
    int Proc_num_pair;

    if(p > 1)
    {
        bi = K + p*(aA->cur_str) + (Nomer >= two_in_the_power_of_a ? 0 : p - b);
        if((Nomer < b || Nomer >= two_in_the_power_of_a) && ((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || (bi == up_bound - 1 && l > 0)))
        {
            block_size_row = (prev_s < k ? m : l);
            pa = A + s*m*n + shag*block_size_row*m;
            pb = B + s*m*n;
            
            Proc_num_pair = (Nomer + shag + (Nomer < b? two_in_the_power_of_a : -two_in_the_power_of_a))%p;
            
            MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
            MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
            if(Nomer < b)
            {
                down_block_size_row = (bi < k ? m : l);
                pa_down = buf_pa;
                pb_down = buf_pb;  
            }
            else
            {
                down_block_size_row = (bi < k ? m : l);
                block_size_row = m;
                pa_down = pa;
                pb_down = pb;
                pa = buf_pa;
                pb = buf_pb;
            }
            ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, true) ;

            for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
            {
                down_block_size_col = (bj < k ? m : l);
                            
                ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, true);

            }
            
            for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
            {
                down_block_size_col = (bj < k ? m : l);

                /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                {*/
                    ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, true);
                /*ZerosMatrix[s*(k+1) + bj] = 1;
                    if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                    {
                        ZerosMatrix[bi*(k + 1) + bj] = 1; 
                    }
                    else
                    {
                        ZerosMatrix[bi*(k + 1) + bj] = 0; 
                    }
                    
                }*/
            }
            
            if(Nomer >= two_in_the_power_of_a)
            {
                return;
            }
        }

        Nomer++;
        bi = (K + p*(aA->cur_str) + (Nomer % 2 == 1 ? 1 : 0));

        if(((Nomer - 1) < (p-b)) && ((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || (bi == up_bound - 1 && l > 0)))
        {
            block_size_row = (prev_s < k ? m : l);
            pa = A + s*m*n + shag*block_size_row*m;
            pb = B + s*m*n;

            Proc_num_pair = (Nomer + shag - (Nomer % 2 == 1 ? 0 : 2))%p;
            
            MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
            MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);

            down_block_size_row = (bi < k ? m : l);

            if (Nomer % 2 == 1)
            {
                pa_down = buf_pa;
                pb_down = buf_pb; 
            }
            else
            {
                pa_down = pa;
                pb_down = pb;
                pa = buf_pa;
                pb = buf_pb;
                block_size_row = m;
            }

            ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, true);
                
            for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
            {
                down_block_size_col = (bj < k ? m : l);
                            
                ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, true);

            }
            
            for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
            {
                down_block_size_col = (bj < k ? m : l);
                
                /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                {*/
                    ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, true);
                    /*ZerosMatrix[s*(k+1) + bj] = 1;
                    if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                    {
                        ZerosMatrix[bi*(k + 1) + bj] = 1; 
                    }
                    else
                    {
                        ZerosMatrix[bi*(k + 1) + bj] = 0; 
                    }
                    
                }*/
            }
            
            if(Nomer % 2 == 0)
            {
                return;
            }
            
        }


        for (step = 1; step < a; step++)
        {
            x = pow(2,step);
            bi = (Nomer%(2*x) == 1 ? x : 0) + Nomer + shag - 1;
            block_size_row = (bi < k ? m : l);
            if(((Nomer - 1) < (p-b)) && ((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || (bi == up_bound - 1 && l > 0)))
            {
                pa = A + s*m*n + shag*block_size_row*m;
                pb = B + s*m*n;

                Proc_num_pair = (Nomer + shag - 1 + (Nomer%(2*x) == 1 ? x : -x))%p;
                
                MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                if(Nomer % (2*x) == 1)
                {
                    pa_down = buf_pa;
                    pb_down = buf_pb;
                    down_block_size_row = (bi < k ? m : l);
                }
                else
                {
                    pa_down = pa;
                    pb_down = pb;
                    pa = buf_pa;
                    pb = buf_pb;
                    down_block_size_row = block_size_row;
                    block_size_row = m;
                }
                
                ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, true);

                for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                                
                    ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, true);

                }
                
                for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                {
                    down_block_size_col = (bj < k ? m : l);

                    /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                    {*/
                        ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, true);
                        /*ZerosMatrix[s*(k+1) + bj] = 1;
                        if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 1; 
                        }
                        else
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 0; 
                        }   
                        
                    }*/
                } 
                
                if(Nomer %(2*x) != 1)
                {
                    return;
                }
            }

        }

    }

    //Fird part
    if(K == shag%p)
    {   
        block_size_row = (prev_s < k ? m : l);     
        pa = A + s*m*n + shag*block_size_row*m;

        if(InverseTriungleBlock(pa, U, block_size_row, norm) != 0)
        {
            cout<<"Matrix is singular"<<endl;
            aA->res = 1;
        }
        else
        {        
            for ( j = shag+1; j < up_bound; j++)
            {
                block_size_col = (j < k ? m : l);
                pa += m*m;
                
                BlockMul(U, pa, ProductResult, block_size_row, block_size_row, block_size_col); 
                ReplaceWith(pa, ProductResult, block_size_row, block_size_col);
            }


            pb = B + s*m*n;
            
            for( j = 0; j < up_bound; j++, pb += m*block_size_row)
            {
                block_size_col = (j < k ? m : l);
                if ((size == 0 || size == 4) || ((size == 3) && (j == 0)) || (j + 1 >= s))
                {
                    BlockMul(U, pb, ProductResult, block_size_row, block_size_row, block_size_col);
                    ReplaceWith(pb, ProductResult, block_size_row, block_size_col);
                }
                else
                {
                    ReplaceWith(pb, ZeroMatrix, block_size_row, block_size_col);
                }
                
            } 
        }
        
    }    
}

void SecondStep(Args* aA)
{
    double* A = aA->A, *B = aA->B, *U = aA->U, *ProductResult = aA->ProductResult;
    double norm = aA->norm;
    int n = aA->n, m = aA->m, p = aA->active_procceses, K = aA->k, shag = aA->shag;
    int a = int(log2(p));
    int two_in_the_power_of_a = pow(2,a);
    int b = p - two_in_the_power_of_a;
    int bi, bj, j;
    int x;
    MPI_Status st;
    MPI_Comm comm = aA->comm;

    int l = n%m;
    int k = (n-l)/m;
    int block_size_row, block_size_col, size = aA->s, down_block_size_row, down_block_size_col;
    double* pa, *pa_side, *pa_down, *pa_down_side, *pb, *pb_down; 
    double* buf_pa = aA->buf, *buf_pb = (aA->buf + n*m);
    double* ZeroMatrix = aA->ZeroMatrix;
    bool* ZerosMatrix = aA->ZerosMatrix;
    
    double eps = (aA->norm > 10 && size != 4 ? EPSILON : aA->norm);
    
    int up_bound = (l > 0 ? k+1: k);

    int s = aA->cur_str;
    int prev_s = K + p*s;

    int step;

    int send_size_a = /*n*m;*/((l == 0 ? k : k+1) - shag)*m*m;
    int send_size_b = n*m;

    int Nomer = aA->nomer_v_okne;
    int Proc_num_pair;
    bool down_is_triangle = false;

    if(p > 1)
    {
        if (Nomer < b)
        {
            bi = (K + p*(aA->cur_str) + p - b);
            block_size_row = (prev_s < k ? m : l);

            if((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || (bi == up_bound - 1 && l > 0))
            {   
                pa = A + s*m*n + shag*block_size_row*m;
                pb = B + s*m*n;

                Proc_num_pair = (Nomer + two_in_the_power_of_a + shag)%p;

                MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);

                down_block_size_row = (bi < k ? m : l);
                pa_down = buf_pa;
                pb_down = buf_pb;
                down_is_triangle = !(bi == up_bound - 1 && l > 0);         
                ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, down_is_triangle);

                for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                                
                    ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);

                }

                for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                {
                    down_block_size_col = (bj < k ? m : l);

                    /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                    {*/
                        ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);
                    /*ZerosMatrix[s*(k+1) + bj] = 1;
                        if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 1; 
                        }
                        else
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 0; 
                        }
                        
                    }*/
                }
            }
        }
        else if(Nomer >= two_in_the_power_of_a)
        {
            bi = (K + p*(aA->cur_str));
            block_size_row = (prev_s < k ? m : l);

            if((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || (bi == up_bound - 1 && l > 0))
            {  
                pa = A + s*m*n + shag*block_size_row*m;
                pb = B + s*m*n;

                Proc_num_pair = (Nomer - two_in_the_power_of_a + shag)%p;

                MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);  
                down_block_size_row = (bi < k ? m : l);
                block_size_row = m;
                pa_down = pa;
                pb_down = pb;
                pa = buf_pa;
                pb = buf_pb;    
                down_is_triangle = !(bi == up_bound - 1 && l > 0);         

                ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, down_is_triangle) ;

                for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                                
                    ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);

                }

                for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                {
                    down_block_size_col = (bj < k ? m : l);

                    /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                    {*/
                        ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);
                    /*ZerosMatrix[s*(k+1) + bj] = 1;
                        if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 1; 
                        }
                        else
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 0; 
                        }
                        
                    }*/
                }
            }
            
            return; 
        }

        Nomer++;

        if (Nomer % 2 == 1 && (Nomer - 1) < (p-b))
        {
            bi = (K + p*(aA->cur_str) + 1);
            block_size_row = m;

            if((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) ||  ((bi == up_bound - 1) && (l > 0)))
            { 
                pa = A + s*m*n + shag*block_size_row*m;
                pb = B + s*m*n;

                Proc_num_pair = (Nomer + shag)%p;
                
                MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                down_block_size_row = (bi < k ? m : l);
                pa_down = buf_pa;
                pb_down = buf_pb; 
                down_is_triangle = !(bi == up_bound - 1 && l > 0);         

                down_block_size_row = (bi < k ? m : l);
                //pa_down = A + bi*m*n + shag*down_block_size_row*m;
                
                ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, down_is_triangle);
                
                for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                                
                    ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);

                }

                for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                    
                    /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                    {*/
                        ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);
                        /*ZerosMatrix[s*(k+1) + bj] = 1;
                        if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 1; 
                        }
                        else
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 0; 
                        }
                        
                    }*/
                } 
            }
            
        }
        else if (Nomer % 2 == 0 && (Nomer - 1) < (p-b))
        {
            bi = (K + p*(aA->cur_str));
            block_size_row = (prev_s < k ? m : l);

            if((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || ((bi == up_bound - 1) && (l > 0)))
            {
                pa = A + s*m*n + shag*block_size_row*m;
                pb = B + s*m*n;

                Proc_num_pair = (Nomer - 2 + shag)%p;

                MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                
                pa_down = pa;
                pb_down = pb;
                pa = buf_pa;
                pb = buf_pb;
                block_size_row = m;
                down_block_size_row = (bi < k ? m : l);
                //pa_down = A + bi*m*n + shag*down_block_size_row*m;
                down_is_triangle = !(bi == up_bound - 1 && l > 0);         

                ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, down_is_triangle);
                
                for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                                
                    ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);

                }

                for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                {
                    down_block_size_col = (bj < k ? m : l);
                    
                    /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                    {*/
                        ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);
                        /*ZerosMatrix[s*(k+1) + bj] = 1;
                        if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 1; 
                        }
                        else
                        {
                            ZerosMatrix[bi*(k + 1) + bj] = 0; 
                        }
                        
                    }*/
                } 
            }
            return;
        }

        for (step = 1; step < a; step++)
        {
            x = pow(2,step);

            if (Nomer % (2*x) == 1 && (Nomer - 1) < (p - b))
            {
                bi = x + Nomer + shag - 1;
                block_size_row = m;

                if((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound) || (bi == up_bound - 1 && l > 0)) 
                {
                    pa = A + s*m*n + shag*block_size_row*m;
                    pb = B + s*m*n;

                    Proc_num_pair = (Nomer + x + shag - 1)%p;

                    MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                    MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                    
                    pa_down = buf_pa;
                    pb_down = buf_pb;
                    down_block_size_row = (bi < k ? m : l);
                    //pa_down = A + bi*m*n + shag*down_block_size_row*m;
                    down_is_triangle = !(bi == up_bound - 1 && l > 0);         

                    ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, down_is_triangle);

                    for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                    {
                        down_block_size_col = (bj < k ? m : l);
                                    
                        ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);

                    }

                    for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                    {
                        down_block_size_col = (bj < k ? m : l);

                        /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                        {*/
                            ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);
                            /*ZerosMatrix[s*(k+1) + bj] = 1;
                            if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                            {
                                ZerosMatrix[bi*(k + 1) + bj] = 1; 
                            }
                            else
                            {
                                ZerosMatrix[bi*(k + 1) + bj] = 0; 
                            }   
                            
                        }*/
                    } 
                }
            }
            else if (Nomer % (2*x) != 1 && (Nomer - 1) < (p - b))
            {
                bi = Nomer + shag - 1;
                block_size_row = (bi < k ? m : l);

                if((bi < up_bound - 1 && l > 0) || (l == 0 && bi < up_bound)||(bi == up_bound - 1 && l > 0) )
                {
                    pa = A + s*m*n + shag*block_size_row*m;
                    pb = B + s*m*n;

                    Proc_num_pair = (Nomer - x + shag - 1)%p;
                    
                    MPI_Sendrecv(pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, buf_pa, send_size_a, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                    MPI_Sendrecv(pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, buf_pb, send_size_b, MPI_DOUBLE, Proc_num_pair, 0, comm, &st);
                    
                    pa_down = pa;
                    pb_down = pb;
                    pa = buf_pa;
                    pb = buf_pb;
              
                    down_block_size_row = block_size_row;
                    block_size_row = m;
                    //pa_down = A + bi*m*n + shag*down_block_size_row*m;
                    down_is_triangle = !(bi == up_bound - 1 && l > 0);         

                    ZeroOut(pa, pa_down, U, m, down_block_size_row, norm, down_is_triangle);

                    for (bj = shag+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < up_bound; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
                    {
                        down_block_size_col = (bj < k ? m : l);
                                    
                        ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);

                    }

                    for (bj = 0; bj < up_bound; bj++, pb += m*m, pb_down += down_block_size_row*m)
                    {
                        down_block_size_col = (bj < k ? m : l);

                        /*if (ZerosMatrix[s*(k+1) + bj] || ZerosMatrix[bi*(k + 1) + bj])
                        {*/
                            ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row, false, down_is_triangle);
                            /*ZerosMatrix[s*(k+1) + bj] = 1;
                            if (!MatrixIsZero(pb_down, down_block_size_col, down_block_size_row, eps))
                            {
                                ZerosMatrix[bi*(k + 1) + bj] = 1; 
                            }
                            else
                            {
                                ZerosMatrix[bi*(k + 1) + bj] = 0; 
                            }   
                            
                        }*/
                    } 
                }
               
                return;
            }
            
        }
    }
    //Fird part
    if(K == shag%p)
    {   
        block_size_row = (prev_s < k ? m : l);     
        pa = A + s*m*n + shag*block_size_row*m;

        if(InverseTriungleBlock(pa, U, block_size_row, norm) != 0)
        {
            cout<<"Matrix is singular"<<endl;
            aA->res = 1;
        }
        else
        {        
            for ( j = shag+1; j < up_bound; j++)
            {
                block_size_col = (j < k ? m : l);
                pa += m*m;
                
                BlockMul(U, pa, ProductResult, block_size_row, block_size_row, block_size_col); 
                ReplaceWith(pa, ProductResult, block_size_row, block_size_col);
            }


            pb = B + s*m*n;
            
            for( j = 0; j < up_bound; j++, pb += m*block_size_row)
            {
                block_size_col = (j < k ? m : l);
                if ((size == 0 || size == 4) || ((size == 3) && (j == 0)) || (j + 1 >= s))
                {
                    BlockMul(U, pb, ProductResult, block_size_row, block_size_row, block_size_col);
                    ReplaceWith(pb, ProductResult, block_size_row, block_size_col);
                }
                else
                {
                    ReplaceWith(pb, ZeroMatrix, block_size_row, block_size_col);
                }
                
            } 
        }
        
    }    

}

void ThirdStep(Args *a)
{
    double* A = a->A, *B = a->B, *pb, *buf = a->buf;
    int n = a->n, m = a->m, p = a->p, K = a->k;
    int l = n%m;
    int bj, s;
    int k = (n-l)/m;
    int bi,  r, new_s;
    int rows = a->rows;
    int owner;
    int block_size_col, size = a->s, block_size_row;
    int up_bound = (l > 0 ? k+1 : k);
    MPI_Comm comm = a->comm;
    MPI_Status st;
   
    s = K + p*(a->cur_str - 1);
    new_s = a->cur_str-1;
    for (bi = up_bound - 1; bi > 0; bi -= 1)
    {
        owner = bi%p;
        block_size_row = (bi< k ? m : l);
        if(owner != K)
        {
            pb = buf;
        }
        else
        {
            pb = B + new_s*m*n;
        }

        MPI_Bcast(pb, n*block_size_row, MPI_DOUBLE, owner, comm);
        if(bi <= s)
        {
            s = s - p;
            new_s =  new_s - 1;
        }
        for (bj = 0; bj < up_bound; bj++)
        {
            block_size_col = (bj < k ? m : l);      

            for (r = new_s ; r >= 0; r --)
            {
                /*if ((size == 3 && (r <= bj || bj == 0)) || size != 3)
                {*/
                    MinusEqualBlockMul(B + r*m*n + bj*m*m, A + r*m*n + bi*m*m, pb + bj*block_size_row*m, m, block_size_row, block_size_col);
                //}
            }
            
        }
    }    
}

int InverseMatrixParallel(Args* a)
{
    int K = a->K, p = a->p;
    int k = a->k;
    int m = a->m;
    int l = a->l;
    int bi, rows = a->rows, up_bound = (a->l == 0 ? K : K+1);
    int res_l, res_g = 0; 
    int cur_global_str = k;
    int send_size;

    for (bi = 0; bi < up_bound; bi++)
    {
        a->shag = bi;
        a->send_size =(K - bi)*m*m + (l == 0 ? 0 : m*l);

        if(rows > 0)
        {
            FirstStep(a);
        
            SecondStep(a);
        }
        res_l = a->res;     
        
        MPI_Allreduce(&res_l, &res_g, 1, MPI_INT, MPI_MAX, a->comm);
        a->res = res_g;
        if (a->res > 0)
        {
            return 1;
        }

        if (k == bi%p)
        {
            a->cur_str++;
        }
        a->nomer_v_okne = (((a->nomer_v_okne-1)%p) + p)%p;

    }
    ThirdStep(a);
    
    return 0;
}
