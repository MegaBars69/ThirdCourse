#ifndef HEADER1  
#define HEADER1

double Trace(double* A, int n);
double LengthOfMatrix(double* A, int n);
double Norm(double* A, double* results, int n);
void ApplyLeftVector(double* X, double* A, int n, int k, bool inside);
void AppyRightVector(double* X, double* A, int n, int k, bool inside);
void TriDiagonalize(double* A, double* U, int n, double mera, double* Y);
void TriDiagonalize(double* A, double* U, int n, double mera);
int FindEigenValues(double* A, int n, double* X, double eps);
void UAU(double* A, int n, int k, double* X, double* Y);
//int FindEigenValues(double* A, int n, double* X, double* Y, double* Z , double eps);
#endif 