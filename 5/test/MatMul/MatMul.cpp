#include <iostream>
#include "MatMul.hpp"
#include <cmath>
using namespace std;

void MatMul1(double *a, double* b, double* c, int n){
    int i, j, m;
    double *pa, *pb, *pc, s;
    for (m = 0, pc = c; m < n; m++)
    {
        for (i = 0, pb = b; i < n; i++, pb++)
        {
            for (s = 0., j = 0, pa = a + m*n; j < n; j++)
                s += *(pa++) * pb[j*n];
                *(pc++) = s;                    
        }        
    }    
}

void MatMul2(double *a, double* b, double* c, int n){
    int i, j, m;
    double *pa, *pb, *pc, s;

    for (m = 0, pc = c; m < n; m++, pc++)
    {
        for (i = 0, pa = a, pb = b + m; i < n; i++)
        {
            for (s = 0., j = 0; j<n; j++)
                s+=*(pa++) * pb[j*n];
            pc[i*n] = s;
        }
        
    }
}

void MatMul3(double *a, double* b, double* c, int n){
    int i, bi, bm, nbi, nbm, j, m, N = 10;
    double *pa, *pb, *pc, s;

    for (bm = 0; bm < n; bm+=N)
    {
        nbm = (bm + N <= n ? bm + N: n);
        for (bi = 0; bi < n; bi += N)
        {
            nbi = (bi + N <= n ? bi + N: n);
            for (m = bm, pc =c+bm; m < nbm; m++, pc++)
            {
                for (i = bi, pa = a + bi*n, pb = b + m; i < nbi; i++)
                {
                    for (s = 0., j = 0; j < n; j++)                   
                        s += *(pa++) * pb[j*n];
                    pc[i*n] = s;
                }                
            }
            
        }
        
    }
}

void Print(double* a, int n){

    cout<<"=============="<<endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout<<a[i*n+j]<<" ";
        }
        cout<<endl;
    }   
    cout<<"=============="<<endl;
}

void Test(double*A, double*B, double* C, int n, void (*func)(double*, double*, double*, int))
{
    clock_t start = clock();
    func(A,B,C, n);
    clock_t end = clock();

    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Время выполнения: " << elapsed << " секунд" << std::endl;
}
;
void MatMul4(double *a, double* b, double* c, int n){
    int i,j,s,r,t,q, v, h, ah;
    int m = 3;
    double sum, s00, s01, s02, s10,s20, s11,s12,s21,s22;
    int k = n/m, l = n-k*m;
    int bl = (l!=0 ? k: k+1);//block amount;
    double *pa, *pb, *pc;
    for (i = 0; i < bl; i++)
    {
        for (j = 0; j < bl; j++)
        {
            v = (i<k ? m:l); h = (j<k ? m:l);//Size of c;
            pc = c + i*m*n + j*m;
            for (r = 0; r < v; r++)
            {
                for (t = 0; t < h; t++)
                {
                    pc[r*v+t] = 0;//C to zero;
                }
                
            }
            for (s = 0; s < bl; s++)
            {
                ah = (s < k ? m:l);
                pa = a + i*m*n + s*m;
                pb = b + s*m*n + j*m;
                int v3 = v%3;
                int h3 = h%3;
                for (r = 0; r < v3; r++)
                {
                    for (t = 0; t < h3; t++)
                    {
                        sum = 0;
                        for (q = 0; q < ah; q++)
                        {
                            sum +=pa[r*n + q]*pb[q*n + t];
                        }
                        pc[r*n+t]+=sum;
                    }
                    for (;t < h; t+=3)
                    {
                        s00 = 0; s01 = 0; s02 = 0;
                        for (q = 0; q < ah; q++)
                        {
                            s00 += pa[r*n + q]*pb[q*n + t];
                            s01 += pa[r*n + q]*pb[q*n + t + 1];
                            s02 += pa[r*n + q]*pb[q*n + t + 2];
                        }
                        pc[r*n + t] += s00; pc[r*n + t + 1] += s01; pc[r*n + t + 2] += s02;
                    }
                }
                for (; r < v; r+=3)
                {
                    for (t = 0; t < h3; t++)
                    {
                        s00 = 0; s10 = 0; s20 = 0;
                        for (q = 0; q < ah; q++)
                        {
                            s00 += pa[r*n + q]*pb[q*n + t];
                            s10 += pa[(r+1)*n + q]*pb[q*n + t];
                            s20 += pa[(r+2)*n + q]*pb[q*n + t];
                        }
                        pc[r*n + t] += s00; pc[(r+1)*n + t] += s10; pc[(r+2)*n + t] += s20;
                    }          
                }
                for (;t < h; t+=3)
                {
                    s00 = 0; s01 = 0; s02 = 0;
                    s10 = 0; s11 = 0; s12 = 0;
                    s20 = 0; s21 = 0; s22 = 0;
                    for (q = 0; q < ah; q++)
                    {
                        s00 += pa[r*n + q]*pb[q*n + t];
                        s01 += pa[r*n + q]*pb[q*n + t + 1];
                        s02 += pa[r*n + q]*pb[q*n + t + 2];
                        s10 += pa[(r+1)*n + q]*pb[q*n + t];
                        s11 += pa[(r+1)*n + q]*pb[q*n + t+1];
                        s12 += pa[(r+1)*n + q]*pb[q*n + t+2];
                        s20 += pa[(r+2)*n + q]*pb[q*n + t];
                        s21 += pa[(r+2)*n + q]*pb[q*n + t+1];
                        s22 += pa[(r+2)*n + q]*pb[q*n + t+2];
                    }
                    pc[r*n + t] += s00; pc[r*n + t + 1] += s01; pc[r*n + t + 2] += s02;
                    pc[(r+1)*n + t] += s10; pc[(r+1)*n + t+1] += s11; pc[(r+1)*n + t+2] += s12;
                    pc[(r+2)*n + t] += s20; pc[(r+2)*n + t + 1] += s21; pc[(r+2)*n + t + 2] += s22;
                }
                
            }
            
        }
        
    }
    
}







