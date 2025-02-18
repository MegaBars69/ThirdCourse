#ifndef HEADER2  
#define HEADER2

int get_max_rows(int n, int m, int p);
int get_rows(int n, int m, int p, int k);
void FormulaMatrixInitialization(double* A, int n, int m, int p, int K, int s);
void PrintLocalMatrix(double* A,  int n, int m, int p, int K, int r, bool exp_format, bool okruglenie);
void PrintAllData(double* A, int n, int m, int p, int K);
int ReadMatrixFromFile(double* A, int n, int m, int p, int K, char* file_name, double* buf, MPI_Comm comm);
void BuildE(double* B, int n, int m, int p, int K);
void PrintMatrix(double* matrix, int n, int m, int p, int K, int r, double* buf, MPI_Comm comm);

#endif
