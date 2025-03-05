#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include "initialize_matrix.hpp"
#include "algorithm.hpp"

using namespace std;

//Доступ на элемент a_{11}^(i,j) блоку A_{i,j} достигается обращением к указателю ∗(a+ i ∗ (k ∗m∗m+m∗ l) + j ∗ m∗block_size_row)

/*
Function to initialize matrix using formula
a_{i,j} = In Block_{i//m, j//m} in {i%col_size, j%row_size} position
i = bi*m + i'
j = bj*m + j'
*/

double f(int i, int j, int n, int s)
{
    switch (s) {
    case 1:
        return  n - ((i < j) ? j : i);
        break;
    case 2:
        return 1 + ((i >= j) ? i : j);
        break;
    case 3:
        return (i - j >= 0 ? (i - j) : (j - i));
        break;
    case 4:
        return 1/double(i+j+1);
        break;
    
    default:
        printf("Неправильный формат ввода s \n");
        break;
    }
    return -1;
}

int l2g(int m, int p, int k, int i_l)
{
    int i_l_m = i_l/m;
    int i_g_m = i_l_m*p + k;
    return i_g_m*m + i_l%m; 
}

int g2l(int m, int p, int i_g)
{
    int i_g_m = i_g/m;
    int i_l_m = i_g_m%p;
    return i_l_m*m+i_g%m;
}

int get_max_rows(int n, int m, int p)
{
    int b = (m+n-1)/m;
    return (b+p-1)/p;
}

int get_rows(int n, int m, int p, int k)
{
    int b = (m+n-1)/m;
    return (k >= b%p ? b/p : (b + p - 1)/p);
}

int get_k(int m, int p, int i_g)
{
    int i_g_m = i_g/m;
    return i_g_m % p;
}

void PrintMatrix(double* A, int n, int m, int p, int K, int r, double* buf, MPI_Comm comm)
{
    int l = n % m;
    int k = (n - l) / m;
    int bi, bj, i_, j_;
    int row_block_size, col_block_size;
    int line_in_block, elem_in_line;
    int row_bound = r / m;
    int col_bound = r - m * row_bound;
    int start_index, owner;

    if (K == 0) 
    {
        cout << endl;
    }

    for (bi = 0; bi < row_bound + 1; bi++) 
    {
        row_block_size = (bi < k) ? m : l;
        line_in_block = (bi < row_bound) ? m : col_bound;
        owner = bi % p;
        int line_of_blocks_loc = bi / p;

        if (K == 0) 
        {
            if (owner == 0) 
            {
                memcpy(buf, A + line_of_blocks_loc * (m * m * k + l * m), (k * row_block_size * m + row_block_size * l) * sizeof(double));
            } 
            else 
            {
                MPI_Status st;
                MPI_Recv(buf, k * row_block_size * m + row_block_size * l, MPI_DOUBLE, owner, 0, comm, &st);
            }
        } 
        else 
        {
            if (owner == K) 
            {
                MPI_Send(A + line_of_blocks_loc * (m * m * k + l * m), k * row_block_size * m + row_block_size * l, MPI_DOUBLE, 0, 0, comm);
            }
        }

        if (K == 0) 
        {
            for (i_ = 0; i_ < line_in_block; i_++) 
            {
                for (bj = 0; bj < row_bound + 1; bj++) 
                {
                    col_block_size = (bj < k) ? m : l;
                    elem_in_line = (bj < row_bound) ? m : col_bound;
                    
                    start_index = bj * m * row_block_size + i_ * col_block_size;

                    for (j_ = 0; j_ < elem_in_line; j_++) 
                    {
                        printf("%10.3e ", *(buf + start_index + j_));
                    }
                }
                cout << endl;
            }
        }
    }

    if (K == 0) 
    {
        cout << endl;
    }
}


void PrintLocalMatrix(double* A,  int n, int m, int p, int K, int r)
{
    if (K==0)
    {
        cout<<endl<<"Local Matrix"<<endl;
    }
    int l = n%m;
    int k = (n-l)/m;
    int max_p;
    int printed_strings = 0, printed_el = 0;

    for (int bi_loc = 0, bi = K; bi < k + 1; bi_loc++, bi += p)
    {
        max_p = (bi < k ? m : l);
        
        for (int p = 0; p < max_p; p++)
        {
            if(printed_strings == r)
            {
                break;
            }
            for (int bj = 0; bj < k+1; bj++)
            {
                int max_l = (bj<k ? m : l);
                double* pa = (A + bi_loc * (k*m*m + m*l) + bj*max_p*m + p*max_l); 
                for (int l = 0; l < max_l; l++)
                {
                    printf("%10.3e ", *(pa++));
                    printed_el++;

                    if(printed_el == r)
                    {
                        bj = k+1;
                        l = max_l;
                        printed_el = 0;
                    }
                }
            }

            printed_strings++;
            cout<<endl;
        }
    }    
    cout<<endl;
}


void PrintAllData(double* A, int n, int m, int p, int K)
{
    int rows = get_rows(n,m,p,K);
    for (int i = 0; i < rows*m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout<<A[i*n+j]<<" ";
        }
        cout<<endl;
    }
}

void FormulaMatrixInitialization(double* A, int n, int m, int p, int K, int s)
{
    int i_glob, j_glob, rows;
    int bi, bi_loc, bj, i_, j_;

    int l = n%m;
    int k = (n-l)/m;
    int col_block_size, row_block_size;

    double *pa = A;

    rows = get_rows(n, m, p, K);
    
    for (bi = K, bi_loc = 0; bi_loc < rows; bi+=p, bi_loc++)
    {
        row_block_size = (bi < k ? m : l);
        for (bj = 0; bj < k + 1; bj++)
        {
            col_block_size = (bj < k ? m : l);
            
            for (i_ = 0; i_ < row_block_size; i_++)
            {
                for (j_ = 0; j_ < col_block_size; j_++, pa++)
                {
                    i_glob = bi*m + i_;
                    j_glob = bj*m + j_;
                
                    *(pa) = f(i_glob, j_glob, n, s);
                    
                }
    
            }  
        }
        
    }
}

void BuildE(double* A, int n, int m, int p, int K)
{
    int i_glob, j_glob, rows;
    int bi, bi_loc, bj, i_, j_;

    int l = n%m;
    int k = (n-l)/m;
    int col_block_size, row_block_size;

    double *pa = A;

    rows = get_rows(n, m, p, K);
    
    for (bi = K, bi_loc = 0; bi_loc < rows; bi+=p, bi_loc++)
    {
        row_block_size = (bi < k ? m : l);
        for (bj = 0; bj < k + 1; bj++)
        {
            col_block_size = (bj < k ? m : l);
            
            for (i_ = 0; i_ < row_block_size; i_++)
            {
                for (j_ = 0; j_ < col_block_size; j_++, pa++)
                {
                    i_glob = bi*m + i_;
                    j_glob = bj*m + j_;
                
                    *(pa) = (i_glob != j_glob ? 0 : 1);
                    
                }
    
            }  
        }  
    }
}


int ReadMatrixFromFile(double* A, int n, int m, int p, int K, char* file_name, double* buf, MPI_Comm comm) 
{
    FILE* file = nullptr;
    int err = 0;
    int read_elements = 0;
    int read_strings = 0;
    int l = n % m;
    int k = (n - l) / m;
    int bi = 0, i_ = 0, j_ = 0, bj = 0, row_block_size = 0, col_block_size = 0, row_block_size_loc = 0;
    int el_in_line, owner;
    int res;

    if (K == 0) 
    {
        file = fopen(file_name, "r");
        if (!file) 
        {
            err = 1;
        }
    }

    MPI_Bcast(&err, 1, MPI_INT, 0, comm);

    if (err != 0) 
    {
        return err;
    }

    for (bi = 0; bi < k + 1; bi++) 
    {
        row_block_size = (bi < k) ? m : l;

        for (i_ = 0; i_ < row_block_size; i_++) 
        {
            for (bj = 0; bj <= k; bj++) 
            {
                col_block_size = (bj < k) ? m : l;

                el_in_line = bj * row_block_size * m + i_ * col_block_size;

                if (K == 0) 
                {
                    for (j_ = 0; j_ < col_block_size; j_++) 
                    {
                        res = fscanf(file, "%lf", buf + el_in_line + j_);
                        if (res == 1) 
                        {
                            read_elements++;
                        } 
                        else if(res == EOF) 
                        {
                            err = 2;
                            break; 
                        }
                        else if(res == 0)
                        {
                            err = 3;
                            break;
                        }
                        if (fgetc(file) == '\n') {
                            read_strings++;
                        }
                        ungetc('\n', file); 
                    }
                }

                if (err != 0) 
                {
                    break; 
                }
            }

            if (err != 0) 
            {
                break;
            }
        }

        MPI_Bcast(&err, 1, MPI_INT, 0, comm);
        MPI_Bcast(&read_elements, 1, MPI_INT, 0, comm);

        if (err != 0 && read_elements != n * n) 
        {
            if (K == 0) 
            {
                fclose(file);
            }
            return err;
        }

        owner = bi % p;
        row_block_size_loc = bi / p;

        if (K == 0) 
        {
            if (owner == 0) 
            {
                memcpy(A + row_block_size_loc * (m * m * k + l * m), buf, (k * row_block_size * m + l * row_block_size) * sizeof(double));
            } 
            else 
            {
                MPI_Send(buf, k * row_block_size * m + l * row_block_size, MPI_DOUBLE, owner, 0, comm);
            }
        } 
        else 
        {
            if (K == owner) 
            {
                MPI_Status st;
                MPI_Recv(A + row_block_size_loc * (m * m * k + l * m), k * row_block_size * m + l * row_block_size, MPI_DOUBLE, 0, 0, comm, &st);
            }
        }
    }
    
    double fin;

    if(K==0)
    {
        if(fscanf(file, "%lf",&fin) != EOF)
        {
            err = 4;
        }
        else if(read_strings != n)
        {
            err = 5;
        }
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, comm);
    

    if (K == 0) 
    {
        fclose(file);
    }

    return err;
}
