#include <iostream>
#include <cmath>
#include "algorithm.hpp"
#include "initialize_matrix.hpp"

using namespace std;

double Trace(double* A, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        result += A[i*n + i];
    }
    return result;
}

double LengthOfMatrix(double* A, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            result += A[i*n + j]*A[j*n +i];
        }    
    }
    result = sqrt(result);
    return result;
}

double Norm(double* A, double* results, int n)
{
    int m = n;
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

void ApplyLeftVector(double* X, double* A, int n, int k, bool inside)
{
    double s;
    double* px;
    for (int j = (inside ? k+1: 0); j < n; j++)
    {
        s = 0;
        px = X;
        for (int i = k + 1 ; i < n; i++,px++)
        {
            s += (*px) * (A[i*n + j]);
        }

        px = X;
        s *= 2;
        for (int i = k + 1; i < n; i++,px++)
        {
            A[i*n + j]-= s*(*px);
        }        
    }   
}

void AppyRightVector(double* X, double* A, int n, int k, bool inside)
{
    double s;
    double* px;
    for (int j = (inside ? k+1: 0); j < n; j++)
    {
        s = 0;
        px = X;
        for (int i = k + 1 ; i < n; i++,px++)
        {
            s += (*px) * (A[j*n + i]);
        }

        px = X;
        s *= 2;
        for (int i = k + 1; i < n; i++,px++)
        {
            A[j*n + i]-= s*(*px);
        }        
    }   
}

void TriDiagonalize(double* A, double* U, int n, double mera)
{
    double sk, ajk,akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;

    for (int k = 0; k < n - 1; k++)
    {
        //Finding vector xk
        pu = U;
        sk = 0;

        for (int i = k+2; i < n; i++)
        {
            ajk = A[i*n + k];
            pu++;
            *pu = ajk;
            sk+= ajk*ajk;
            A[i*n + k] = 0;
            A[k*n + i] = 0;
        }

        akk = A[(k + 1)*n + k];
        new_diag_el = sqrt(sk + akk*akk);
        first_in_x = akk - new_diag_el;
        
        pu = U;
        *pu = first_in_x;
        
        norm_xk = sqrt(sk + first_in_x*first_in_x);

        if(!(fabs(norm_xk) < mera))
        {
            for (int i = 0; i < n - k - 1; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
            }
        }

        //Applying Reflection
        A[(k + 1)*n + k] = new_diag_el;
        A[k*n + k + 1] = new_diag_el;

        ApplyLeftVector(U, A, n, k, true); 
        AppyRightVector(U, A, n, k, true); 
    }
}

void UAU(double* A, int n, int k, double* X, double* Y)
{
    double yi = 0, sk = 0;
    int i, j;
    for ( i = 0; i < k; i++)
    {
        yi = 0;
        for ( j = 0; j < k; j++)
        {
            yi += A[i*n + j]*X[j];
        }
        Y[i] = yi;
    }
    for (i = 0; i < k; i++)
    {
        sk +=X[i]*Y[i];
    }
    
    sk = 2*sk; 
    
    for (i = 0; i < k; i++)
    {
        Y[i] = 2*Y[i] - sk*X[i];
    }

    for (i = 0; i < k; i++)
    {
        for ( j = 0; j < k; j++)
        {
            A[i*n +j] -= (Y[i]*X[j] + X[i]*Y[j]);
        }
    }
}

void TriDiagonalize(double* A, double* U, int n, double mera, double* Y)
{
    double *pa = A;
    double sk, ajk,akk;
    double new_diag_el, norm_xk, first_in_x;
    double* pu;

    for (int k = 0; k < n-1; k++)
    {
        pu = U;
        sk = 0;

        for (int i = k+2; i < n; i++)
        {
            ajk = A[i*n + k];
            pu++;
            *pu = ajk;
            sk+= ajk*ajk;
            A[i*n + k] = 0;
            A[k*n + i] = 0;
        }

        akk = A[(k + 1)*n + k];
        new_diag_el = sqrt(sk + akk*akk);
        first_in_x = akk - new_diag_el;
        
        pu = U;
        *pu = first_in_x;
        
        norm_xk = sqrt(sk + first_in_x*first_in_x);

        if(!(fabs(norm_xk) < mera))
        {
            for (int i = 0; i < n; i++, pu++)
            {
                *pu = (*pu)/norm_xk;
            }
        }
        A[(k + 1)*n + k] = new_diag_el;
        A[k*n + k + 1] = new_diag_el;
        pa += n + 1;
        UAU(pa, n, n - k - 1, U, Y);
    }  
}


int FindEigenValues(double* A, int n, double* X, double eps)
{
    int k;
    int iteration = 0;
    bool is_running = true;
    int up_bound = n;
    double s, ann, ann_1, half_ann_1;
    double x1 = 0, x2 = 0, y1 =  0, y2 = 0, new_diag_el =  0, ak1, ak2, norm_xk, scolar_sum;
    bool aply_to_prev = true, did_reflect = true;
    double a11, a12, a21, a22, D;

    while(is_running)
    {
        if (up_bound > 2)
        {   
            ann = A[(up_bound - 1)*n +  up_bound - 1];
            ann_1 = A[(up_bound - 1)*n +  up_bound  - 2];
            half_ann_1 = ann_1/2;
            s = ((ann > 0 && half_ann_1 > 0) || (ann < 0 && half_ann_1 < 0) ? ann + half_ann_1: ann - half_ann_1);            
            
            if (fabs(ann_1) < eps)
            {
                up_bound--;
                X[up_bound] = ann;
            }
            else
            {
                for (k = 0; k < up_bound - 1; k++)
                {
                    if(k == 0)
                    {
                        A[k*n + k] -= s;
                        A[(k + 1)*n + k + 1] -= s;

                        if (fabs(A[(k + 1)*n + k]) > eps)
                        {
                            //Finding vector xk

                            ak1 = A[k*n + k];
                            ak2 = A[(k + 1)*n + k];

                            new_diag_el = sqrt(ak1*ak1 + ak2*ak2);

                            x1 = ak1 - new_diag_el;
                            x2 = ak2;
                
                            A[k*n + k] = new_diag_el;
                            A[(k + 1)*n + k] = 0;
                            
                            norm_xk = sqrt(x1*x1 + x2*x2);

                            x1 = x1/norm_xk;
                            x2 = x2/norm_xk;

                            //Applying to second column

                            ak1 = A[k*n + k + 1];
                            ak2 = A[(k+1)*n + k + 1];

                            scolar_sum = 2*(ak1*x1 + ak2*x2);

                            A[k*n + k + 1] -= x1*scolar_sum;
                            A[(k+1)*n + k + 1] -= x2*scolar_sum;

                            //Applying to third column
                            
                            if (k + 2 < up_bound)
                            {
                                ak2 = A[(k+1)*n + k + 2];

                                scolar_sum = 2*ak2*x2;

                                A[(k+1)*n + k + 2] -= x2*scolar_sum;
                            }
                            y1 = x1;
                            y2 = x2;
                            aply_to_prev = true;
                        }
                        else
                        {
                            y1 = 0;
                            y2 = 0;
                            aply_to_prev = false;
                        }
                    }
                    else
                    {
                        A[(k + 1)*n + k + 1] -= s;
                    
                        if (fabs(A[(k + 1)*n + k]) > eps)
                        {
                            //Finding vector xk

                            ak1 = A[k*n + k];
                            ak2 = A[(k + 1)*n + k];

                            new_diag_el = sqrt(ak1*ak1 + ak2*ak2);

                            x1 = ak1 - new_diag_el;
                            x2 = ak2;
                
                            A[k*n + k] = new_diag_el;
                            A[(k + 1)*n + k] = 0;
                            
                            norm_xk = sqrt(x1*x1 + x2*x2);

                            x1 = x1/norm_xk;
                            x2 = x2/norm_xk;

                            //Applying to second column

                            ak1 = A[k*n + k + 1];
                            ak2 = A[(k+1)*n + k + 1];

                            scolar_sum = 2*(ak1*x1 + ak2*x2);

                            A[k*n + k + 1] -= x1*scolar_sum;
                            A[(k+1)*n + k + 1] -= x2*scolar_sum;

                            //Applying to third column
                            
                            if (k + 2 < up_bound)
                            {
                                ak2 = A[(k+1)*n + k + 2];

                                scolar_sum = 2*ak2*x2;

                                A[(k+1)*n + k + 2] -= x2*scolar_sum;
                            }
                            did_reflect = true;
                        }   
                        else
                        {
                            did_reflect = false;
                        }
                        //Aplying preview from the right side
                        if (aply_to_prev)
                        {
                            //Aplying to first row
                            ak1 = A[(k - 1)*n + k - 1];
                            ak2 = A[(k - 1)*n + k];

                            scolar_sum = 2*(ak1*y1 + ak2*y2);

                            A[(k - 1)*n + k - 1] -= y1*scolar_sum;
                            A[(k - 1)*n + k] -= y2*scolar_sum;

                            //Applying to second row
                            ak1 = A[k*n + k - 1];
                            ak2 = A[k*n + k];

                            scolar_sum = 2*(ak1*y1 + ak2*y2);

                            A[k*n + k - 1] -= y1*scolar_sum;
                            A[k*n + k] -= y2*scolar_sum;
                            A[(k - 1)*n + k] = A[k*n + k - 1]; 
                        }

                        if(did_reflect)
                        {
                            y1 = x1;
                            y2 = x2;
                            aply_to_prev = true;
                        }
                        else
                        {
                            y1 = 0;
                            y2 = 0;
                            aply_to_prev = false;
                        }
                        A[(k - 1)*n + k - 1] += s;
                    }   
                    if (k == up_bound - 2)
                    {
                        //Aplying preview from the right side

                        //Aplying to first row
                        ak1 = A[k*n + k];
                        ak2 = A[k*n + k + 1];

                        scolar_sum = 2*(ak1*y1 + ak2*y2);

                        A[k*n + k] -= y1*scolar_sum;
                        A[k*n + k + 1] -= y2*scolar_sum;

                        //Applying to second row
                        ak1 = A[(k + 1)*n + k];
                        ak2 = A[(k + 1)*n + k + 1];

                        scolar_sum = 2*(ak1*y1 + ak2*y2);

                        A[(k + 1)*n + k] -= y1*scolar_sum;
                        A[(k + 1)*n + k + 1] -= y2*scolar_sum;
                        //A[k*n + k + 1] = A[(k + 1)*n + k];

                        y1 = x1;
                        y2 = x2;
                    
                        A[k*n + k] += s;
                        A[(k + 1)*n + k + 1] += s;
                    }
                    
                    
                    //PrintMatrix(A, n, n);
                }
                iteration++;
                /*cout<<"Up_bound = "<<up_bound<<endl;
                double trace = Trace(A, n);
                double Length = LengthOfMatrix(A, n);
                cout<<"trA = "<<trace<<endl;
                cout<<"||A|| = "<<Length<<endl;*/
            }
        }
        else
        {
            a11 = *A;
            a12 = A[1];
            a21 = A[n];
            a22 = A[n + 1];

            if (fabs(a21) < eps)
            {
                X[0] = a11;
                X[1] = a22;
            }
            else
            {                
                D = (a11 + a22)*(a11 + a22) - 4*(a11*a22 - a12*a21);
                X[0] = (a11 + a22 + sqrt(D))/2;
                X[1] = (a11 + a22 - sqrt(D))/2;
            }
            is_running = false;
        }
    }

    return iteration;
}
