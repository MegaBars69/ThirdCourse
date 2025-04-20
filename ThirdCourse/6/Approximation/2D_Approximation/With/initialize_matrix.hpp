#ifndef HEADER2  
#define HEADER2
#include <iostream>
#include "algorithm.hpp"

void ij2l (int nx, int /*ny*/, int i, int j, int &l);
void l2ij (int nx, int /*ny*/, int &i,	int &j,	int l); 
void thread_rows (int n, int p, int k, int &i1, int &i2);
int get_off_diag(int nx, int ny,int i, int j, int *I);
int get_len_msr_off_diag(int nx, int ny, double *A, int *I);
int get_len_msr (int nx, int ny);
void fill_I(int nx, int ny, int *I);
int fill_IA(int nx, int ny, double hx, double hy, int *I, double* A, int p, int k);
int IA_ij (int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int *I, double *A);
void fill_B(int n, int nx, int ny, double hx, double hy, double *b, double x0, double y0, int p, int k, double (*f)(double, double));


#endif 