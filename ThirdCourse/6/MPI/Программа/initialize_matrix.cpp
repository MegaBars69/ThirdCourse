#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstring>

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

int l2g(int n, int m, int p, int k, int i_l)
{
    int i_l_m = i_l/m;
    int i_g_m = i_l_m*p + k;
    return i_g_m*m + i_l%m; 
}

int g2l(int n, int m, int p, int k, int i_g)
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

int get_k(int n, int m, int p, int i_g)
{
    int i_g_m = i_g/m;
    return i_g_m % p;
}

void PrintMatrix(double* matrix, int n, int m, int p, int K, int r, double* buf, MPI_Comm comm)
{
    int l = n%m;
    int k = (n-l)/m;
    int bi = 0, bj = 0, i_ = 0, j_ = 0;
    int row_block_size = 0, col_block_size = 0;
    int r_num_line_in_block = 0, r_num_elem_in_line = 0;

    int restrict_k = r / m;
    int restrict_l = r - m * restrict_k;
    int owner;
    
    if(K == 0)
    {
        cout<<endl;
    }
    
    if(r > n)
    {
        r = n;
    }

    for(bi =  0; bi < restrict_k + 1; bi++) 
    {
        row_block_size = (bi < k) ? m : l;
        r_num_line_in_block = (bi < restrict_k) ? m : restrict_l;
        owner = bi % p;
        int line_of_blocks_loc = bi / p;
        if (K == 0) 
        {
            if (owner == 0) 
            {
                memcpy(buf, matrix + line_of_blocks_loc * (m * m * k + l * m), (k * row_block_size * m + row_block_size * l) * sizeof(double));
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
                MPI_Send(matrix + line_of_blocks_loc * (m * m * k + l * m), k * row_block_size * m + row_block_size * l, MPI_DOUBLE, 0, 0, comm);
            }
        }

        if (K == 0) 
        {
            for(i_ = 0; i_ < r_num_line_in_block; i_++) 
            {
                for(bj = 0; bj < restrict_k + 1; bj++) 
                {
                    col_block_size = (bj < k) ? m : l;
                    r_num_elem_in_line = (bj < restrict_k) ? m : restrict_l;
                    int shift = bj * m * row_block_size + i_ * col_block_size;
                    for(j_ = 0; j_ < r_num_elem_in_line; j_++) 
                    {
                        printf(" %10.3e", *(buf + shift + j_));
                    }
                }
                cout<< endl;
            }
        }
    }
    if(K == 0)
    {
        cout<< endl;
    }
}


void PrintLocalMatrix(double* A,  int n, int m, int p, int K, int r, bool exp_format, bool okruglenie)
{
    if (K==0)
    {
        cout<<endl<<"Local Matrix"<<endl;
    }
    int l = n%m;
    int k = (n-l)/m;
    int max_p;
    int printed_strings = 0, printed_el = 0;
    int rows = get_rows(n,m,p,K);

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
                    if(exp_format)
                    {
                        printf("%10.3e ", *(pa++));
                    }
                    else
                    {
                        if(okruglenie)
                        {
                            cout<<fixed<<setprecision(2)<<*(pa++)<<" ";
                        }
                        else
                        {
                            cout<<*(pa++)<<" ";                         
                        }                
                    }
                    printed_el++;

                    if(printed_el == r)
                    {
                        bj = k+1;
                        l = max_l;
                        //p = max_p;
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
    int i_loc, j_loc, i_glob, j_glob, rows;
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
    int i_loc, j_loc, i_glob, j_glob, rows;
    int bi, bi_loc, bj, i_, j_;

    int l = n%m;
    int k = (n-l)/m;
    int col_block_size, row_block_size;

    double *pa = A;

    rows = get_rows(n, m, p, k);

    
    for (bi = K, bi_loc = 0; bi < k + 1; bi+=p, bi_loc++)
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
                    
                    *(pa) = (i_glob == j_glob ? 1 : 0);
                }
    
            }  
        }
        
    }
}


int ReadMatrixFromFile(double* A, int n, int m, int p, int K, char* file_name, double* buf, MPI_Comm comm) 
{
    FILE* file = nullptr;
    int error = 0;
    int read_strings = 0;
    int l = n%m;
    int k = (n-l)/m;
    int bi = 0, i_ = 0, j_ = 0,  bj = 0, row_block_size = 0, col_block_size = 0, row_block_size_loc = 0;
    int el_in_line, owner;
    
    if (K == 0) 
    {
        file = fopen(file_name, "r");
        if (!file) 
        {
            error = 1;
        }
    }

    MPI_Bcast(&error, 1, MPI_INT, 0, comm);

    if (error != 0) 
    {
        return error;
    }

    memset(buf, 0, n * m * sizeof(double));

    for(bi = 0; bi < k + 1; bi++) 
    {
        row_block_size = (bi < k) ? m : l;
        
        for(i_ = 0; i_ < row_block_size; i_++)
        {
            for(bj = 0; bj <= k; bj++) 
            {
                col_block_size = (bj < k) ? m : l;

                el_in_line = bj * row_block_size * m + i_ * col_block_size;
                if (K == 0) 
                {
                    for(j_ = 0; j_ < col_block_size; j_++)
                    {
                        if (fscanf(file, "%lf", buf + el_in_line + j_) == 1) 
                        {
                            read_strings++;
                        } 
                        else 
                        {
                            error = 2;
                            goto error_read_el;
                        }
                    }
                }
            }
        }

        error_read_el:
            MPI_Bcast(&error, 1, MPI_INT, 0, comm);
            MPI_Bcast(&read_strings, 1, MPI_INT, 0, comm);

            if (error != 0 && read_strings != n * n) 
            {
                if (K == 0) 
                {
                    fclose(file);
                }
                return error;
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
    if (K == 0) 
    {
        fclose(file);
        return 0;
    }
    
    return 0;
}

