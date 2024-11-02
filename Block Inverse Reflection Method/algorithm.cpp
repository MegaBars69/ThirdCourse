
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "algorithm.hpp"
#define EPSILON pow(10,-16)
using namespace std;

void PrintMatrix(double* A,  int n, int m, int r, bool exp_format, bool okruglenie)
{
    cout<<endl;
    int l = n%m;
    int k = (n-l)/m;
    int max_p;
    int printed_strings = 0, printed_el = 0;

    for (int bi = 0; bi < k+1; bi++)
    {
        max_p = (bi<k ? m : l);
        
        for (int p = 0; p < max_p; p++)
        {
            if(printed_strings == r)
            {
                break;
            }
            for (int bj = 0; bj < k+1; bj++)
            {
                int max_l = (bj<k ? m : l);
                double* pa = (A + bi * (k*m*m + m*l) + bj*max_p*m + p*max_l); 
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

/*double Norm(double* A, int n, int m)
{
    int l = n%m;
    int k = (n-l)/m;
    int i,j, bi, bj,i_,j_;
    int col_block_size, row_block_size;

    double sum = 0, max_sum = 0, a_ij;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            bi = i/m; bj = j/m;
            i_ = i % m; j_ = j % m;
            col_block_size = (bi<k ? m : l);
            row_block_size = (bj<k ? m : l);
            a_ij = A[bi*n*m + bj*m*col_block_size + i_*row_block_size + j_];
            sum += (a_ij >= 0 ? a_ij : -(a_ij)); 
        }  
        if (sum > max_sum)
        {
            max_sum = sum;
        }
        sum = 0;
        
    }
    return max_sum;
}
*/
double Norm(double* A, double* results, int n, int m)
{
    int l = n%m;
    int k = (n-l)/m;
    int block_size_col, block_size_row;
    double sum_in_col = 0, final_norm = 0, el;

    for (int bj = 0; bj < k+1; bj++)
    {
        block_size_col = (bj < k ? m : l);
        for (int bi = 0; bi < k+1 ; bi++)
        {
            block_size_row = (bi < k ? m : l);

            double* pa = (A + m*n*bi + m*block_size_row*bj);
            for (int j = 0; j < block_size_col; j++)
            {
                sum_in_col = 0;
                for (int i = 0; i < block_size_row; i++)
                {
                    sum_in_col += fabs(pa[i*block_size_col + j]);
                }
                results[j] += sum_in_col;
            }
        }
        for (int j = 0; j < block_size_col; j++)
        {
            el = results[j];
            results[j] = 0;
            final_norm = (el < final_norm ? final_norm : el);
        }   

    }   
    return final_norm;
}

void DifferenceOfBlockMatrix(double* A, double* B, int n, int m)
{
    int l = n%m;
    int k = (n-l)/m;
    int col_block_size, row_block_size;

    for (int bi = 0; bi < k+1; bi++)
    {
        col_block_size = (bi<k ? m : l);

        for (int bj = 0; bj < k+1; bj++)
        {
            row_block_size = (bj<k ? m : l);
            double* pb = (B + bi * m*n + bj*m*col_block_size);
            double* pa = (A + bi * m*n + bj*m*col_block_size);
            //double* pc = (Answer + bi * m*n + bj*m*col_block_size);

            for (int i_ = 0; i_ < col_block_size; i_++)
            {
                for (int j_ = 0; j_ < row_block_size; j_++, pa++, pb++)
                {
                    *(pa) -=  *pb;
                }
                
            }             
        }
        
    }
    
}

void MatrixMinusEqual(double* A, double* B, int row_size, int col_size)
{    
    double* pa = A;
    double* pb = B;

    for (int i = 0; i < row_size; i++)
    {
        for (int j= 0; j < col_size; j++)
        {
            *pa -= *pb;
            pa++;
            pb++;
        }  
    }    
}

void MatrixPlusEqual(double* A, double* B, int row_size, int col_size)
{
    double* pa = A;
    double* pb = B;

    for (int i = 0; i < row_size; i++)
    {
        for (int j= 0; j < col_size; j++)
        {
            *pa += *pb;
            pa++;
            pb++;

        }  
    } 
}

void AddToMatrix(double* A, double* B, int n, int m)
{
    int l = n%m;
    int k = (n-l)/m;
    int col_block_size, row_block_size;

    for (int bi = 0; bi < k+1; bi++)
    {
        col_block_size = (bi<k ? m : l);

        for (int bj = 0; bj < k+1; bj++)
        {
            row_block_size = (bj<k ? m : l);
            double* pb = (B + bi * m*n + bj*m*col_block_size);
            double* pa = (A + bi * m*n + bj*m*col_block_size);
            //double* pc = (Answer + bi * m*n + bj*m*col_block_size);

            for (int i_ = 0; i_ < col_block_size; i_++)
            {
                for (int j_ = 0; j_ < row_block_size; j_++, pa++, pb++)
                {
                    *(pa) = *pa + *pb;
                }
                
            }             
        }
        
    }
    
}

double fabs(double a)
{
    return (a>0 ? a:(-a));
}

double CalcError(double * A, double * InversedA, int n, int m)
{
    double max_sum = 0;
    int l = n%m;
    int kk = (n-l)/m;
    double sum = 0;
    int col_block_size1, col_block_size2, row_block_size1, row_block_size2;
    for (int j = 0; j < n; j++)
	{
        sum = 0;
        for (int i = 0; i < n; i++)
        {
            double cij = 0;
            for (int k = 0; k < n; k++)
            {
                int bi = i/m;
                int bk = k/m;
                int bj = j/m;
                col_block_size1 = (bi<kk ? m : l);
                col_block_size2 = (bk<kk ? m : l);
                row_block_size1 = col_block_size2;
                row_block_size2 = (bj<kk ? m : l);
                int i_ = i - bi*m, k_ = k - bk*m, j_ = j - bj*m;
                cij += A[bi * m*n + bk*m*col_block_size1 +i_*row_block_size1+k_]*InversedA[bk * m*n + bj*m*col_block_size2 +k_*row_block_size2+j_];
                //cout<<A[bi * m*n + bk*m*col_block_size1 +i_*row_block_size1+k_]<<"*"<<InversedA[bk * m*n + bj*m*col_block_size2 +k_*row_block_size2+j_]<<" ";
            } 
            sum+=fabs(cij);
            //cout<<endl<<sum<<"+="<<cij<<endl;
        }
        sum-=1;
        max_sum = (max_sum>sum ? max_sum : sum);
    }
    return max_sum;
}

double Discrepancy(double * A, double * InversedA, double* Column, double* ProductResult, double* Sum, int n, int m)
{
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


void ApplyVector(double* X, double* A, int row_num, int col_num, int k, bool inside)
{
    double s;
    double* px;
    for (int j = (inside ? k+1: 0); j < col_num; j++)
    {
        s = 0;
        px = X;
        for (int i = k; i < row_num; i++,px++)
        {
            s += (*px) * (A[i*col_num + j]);
        }

        px = X;
        s *= 2;
        for (int i = k; i < row_num; i++,px++)
        {
            A[i*col_num + j]-= s*(*px);
        }        
    }   
}

void ApplyMatrix(double* U, double* A, int row_num, int col_num, int amount_of_vectors)
{
    //PrintMatrix(U, row_num, 1, row_num);
    for (int k = 0; k < amount_of_vectors; k++)
    {
        ApplyVector((U + k*row_num), A, row_num, col_num, k, false);
    }
}

int Triungulize(double* A, double* U, int row_num, int col_num, double norm)
{
    double sk, ajk,akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;

    for (int k = 0; k < col_num; k++)
    {
        //Finding vector xk
        sk = 0;
        for (int j = k+1; j < row_num; j++)
        {
            ajk = A[j*col_num + k];
            sk+= ajk*ajk;
        }
        akk = A[k*col_num + k];
        new_diag_el = sqrt(sk + akk*akk);
        
        pu = (U + k*row_num);
        *pu = akk - new_diag_el;
        pu++;
        for (int i = k+1; i < row_num; i++, pu++)
        {
            *pu = A[i*col_num + k];
        }
        pu = (U + k*row_num);
        first_in_x = (*pu);
        norm_xk = sqrt(sk + first_in_x*first_in_x);

        if(!(fabs(norm_xk) < EPSILON*norm))
        {
            for (int i = k; i < row_num; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
            }
        }
        //Applying Reflection
        A[k*col_num + k] = new_diag_el;
        for (int i = k+1; i < row_num; i++)
        {
            A[i*col_num + k] = 0;
        }

        ApplyVector((U + k*row_num), A, row_num, col_num, k, true); 

    }

    return 0;
}


void ZeroOut(double* Diag, double* Down, double* U, int m, int row_size, double norm)
{
    double sk, ajk, akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;
    for (int j = 0; j < m; j++)
    {
        sk = 0;
        for (int i = 0; i < row_size; i++)
        {
            ajk = Down[i*m + j];
            sk+= ajk*ajk;
        }
        akk = Diag[j*m + j];
        new_diag_el = sqrt(sk + akk*akk);
        
        pu = (U + j*(row_size+1));
        *pu = akk - new_diag_el;
        pu++;
        for (int i = 0; i < row_size; i++, pu++)
        {
            *pu = Down[i*m + j];
        }
        pu = (U + j*(row_size+1));
        first_in_x = (*pu);
        norm_xk = sqrt(sk + first_in_x*first_in_x);
        
        if(!(fabs(norm_xk) < EPSILON*norm))
        {
            *pu = (*pu)/norm_xk;
            pu++;
            for (int i = 0; i < row_size; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
            }
        }
        
        //Applying Reflection
        
        Diag[j*m + j] = new_diag_el;
        for (int i = 0; i < row_size; i++)
        {
            Down[i*m + j] = 0;
        }
            
        double s;
        double* px;
        for (int jj = j+1; jj < m; jj++)
        {
            s = 0;
            px = U + j*(row_size+1);
            s += (*px) * (Diag[j*m + jj]);
            px++;
            for (int i = 0; i < row_size; i++, px++)
            {
                s += (*px) * (Down[i*m + jj]);
            }
            px = U + j*(row_size+1);

            Diag[j*m + jj]-= 2*s*(*px);
            px++;
            
            for (int i = 0; i < row_size; i++, px++)
            {
                Down[i*m + jj]-= 2*s*(*px);
            }
            
        }
    }
    
}

void ApplyMatrixToPair(double* U, double* Up, double* Down, int col_size, int row_size, int amount_of_vectors)
{
    double s;
    double* pu;
    for (int vec_num = 0; vec_num < amount_of_vectors; vec_num++)
    {   
        for (int j = 0; j < col_size; j++)
        {
            pu = U + vec_num*(row_size + 1);
            
            s = 0;

            s += (*pu) * (Up[vec_num*col_size + j]);
            for (int i = 0; i < row_size; i++)
            {
                pu++;
                s += (*pu) * (Down[i*col_size + j]);
            }
            s*=2;
            pu = U + vec_num*(row_size + 1);
            Up[vec_num*col_size + j] -= s*(*pu);
            for (int i = 0; i < row_size; i++)
            {
                pu++;
                Down[i*col_size + j] -= s*(*pu);
            }
        }           
    }
    
}

int InverseTriungleBlock(double* A, double* B, int n, double norm, bool right_part_to_E)
{
    double* pa = A, *pb = B;
    double diag_el, sum;
    if(right_part_to_E)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                B[i*n+j] = (i != j ? 0 : 1);
            }
            
        }
    }
    
    for (int i = n-1; i >= 0; i--)
    {
        pa = A + i*n + i;
        pb = B + i*n;
        diag_el = *(pa);
        if (fabs(diag_el) < EPSILON*norm)
        {
            return -1;
        }
        if(right_part_to_E){*(pa) = 1;}
    
        for (int j = 0; j < n; j++)
        {
            if(j > i)
            {
                pa++;
                *pa = (*pa)/diag_el;
                
            }
            *pb = (*pb)/diag_el;
            pb++;
        }
        
        
        for (int j = 0; j < n; j++)
        {
            sum = 0;
            for (int s = i+1; s < n; s++)
            {
                sum += A[i*n + s]*B[s*n + j];
            }
            B[i*n+j]-= sum;
        } 
        pa = A + i*n + i;

        if (right_part_to_E)
        {       
            for (int j = i+1; j < n; j++)
            {
                pa++;
                *pa = 0;
            }
        }
            
    }
    return 0;
}

/*
void BlockMul(double *a, double* b, double* c, int n1, int m12, int n2){
    int i, j, r;
    double s, air, brj;


    for (i = 0; i < n1; i++)
    {
        for (j = 0;  j < n2; j++)
        {
            s = 0;

            for (r = 0; r < m12; r++)
            {
                air = a[i*m12 + r];
                brj = b[r*n2 + j]; 
                s += air*brj; 
            }
            c[i*n2 + j] = s;
        }
        
    }
}
*/

void BlockMul(double *a, double* b, double* c, int n1, int m12, int n2)
{
    int l_col = n2%3;
    int l_row= n1%3;
    int k_col = (n2-l_col)/3;
    int k_row = (n1-l_row)/3; 
    
    double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;

    for (int bi = 0; bi < k_row; bi++)
    {
        for (int bj = 0; bj < k_col; bj++)
        {
            s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;
            for (int r = 0; r < m12; r++)
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
            for (int r = 0; r < m12; r++)
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
        for (int bj = 0; bj < k_col; bj++)
        {
            s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0;
            for (int r = 0; r < m12; r++)
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


int InverseMatrix(double* A, double* B, double* U, double* ProductResult, double norm, int n, int m)
{
    int l = n%m;
    int j, bj;
    int k = (n-l)/m;
    int m12;
    int block_size_row, block_size_col, down_block_size_row, down_block_size_col;
    double* pa, *pa_side, *pa_down, *pa_down_side, *pb, *pb_down; 
    for (int s = 0; s < k+1; s++)
    {
        block_size_row = (s < k ? m : l);
        
        pa = A + s*m*n + s*block_size_row*m;
        
        // First part of algorithm

        Triungulize(pa, U,block_size_row, block_size_row, norm);

        for (j = s+1, pa_side = pa + block_size_row*m; j < k+1; j++, pa_side += block_size_row*m)
        {
            block_size_col = (j < k ? m : l);
            ApplyMatrix(U, pa_side, block_size_row, block_size_col, block_size_row);
        }
        for (j = 0, pb = B + s*m*n; j < k+1; j++, pb += block_size_row*m)
        {
            block_size_col = (j < k ? m : l);
            ApplyMatrix(U, pb, block_size_row, block_size_col, block_size_row);
        }
    
        // Second part of algorithm

        for (int bi = s+1; bi < k+1; bi++)
        {
            down_block_size_row = (bi < k ? m : l);
            pa_down = (A + bi*m*n + s*down_block_size_row*m);

            ZeroOut(pa, pa_down, U, m, down_block_size_row, norm);

            for (bj = s+1, pa_down_side = pa_down + down_block_size_row*m, pa_side = pa + m*m; bj < k+1; bj++, pa_down_side += down_block_size_row*m, pa_side += m*m)
            {
                down_block_size_col = (bj < k ? m : l);
                ApplyMatrixToPair(U, pa_side, pa_down_side, down_block_size_col, down_block_size_row, block_size_row);
            }

            for (bj = 0, pb = B + s*m*n, pb_down = B + bi*m*n; bj < k+1; bj++, pb += m*m, pb_down += down_block_size_row*m)
            {
                down_block_size_col = (bj < k ? m : l);
                ApplyMatrixToPair(U, pb, pb_down, down_block_size_col, down_block_size_row, block_size_row);
            }   
        }  

        //Fird part
        
        pa = A + s*m*n + s*block_size_row*m;

        if(InverseTriungleBlock(pa, U, block_size_row, norm) != 0)
        {
            cout<<"Matrix is singular"<<endl;
            return -1;
        } 

        for (int j = s+1; j < k+1; j++)
        {
            block_size_col = (j < k ? m : l);
            pa += m*m;

            BlockMul(U, pa, ProductResult, block_size_row, block_size_row, block_size_col);
            ReplaceWith(pa, ProductResult, block_size_row, block_size_col);
        }

        pb = B + s*m*n;
        
        for(int j = 0; j < k+1; j++)
        {
            block_size_col = (j < k ? m : l);

            BlockMul(U, pb, ProductResult, block_size_row, block_size_row, block_size_col);
            ReplaceWith(pb, ProductResult, block_size_row, block_size_col);
            pb += m*block_size_row;
        } 
    }

    //Gauss Backward
    
    for (int bi = k-1; bi >= 0; bi--)
    {
        block_size_row = (bi< k ? m : l);

        for (int bj = 0; bj < k+1; bj++)
        {
            block_size_col = (bj < k ? m : l);      
            
            //Zero down U
            for (int i = 0; i < m+1; i++)
            {
                for (int j = 0; j < m+1; j++)
                {
                    U[i*(m+1)+j] = 0;
                } 
            }

            for (int r = bi + 1; r < k+1; r++)
            {
                m12 = (r < k ? m : l);

                BlockMul((A + bi*m*n + r*block_size_row*m), (B + r*m*n + bj*m12*m), ProductResult, block_size_row, m12, block_size_col);
                MatrixPlusEqual(U, ProductResult, block_size_row, block_size_col);
            }
            
            MatrixMinusEqual((B + bi*m*n + bj*block_size_row*m), U, block_size_row, block_size_col);
        } 
    }
    
    return 0;
}