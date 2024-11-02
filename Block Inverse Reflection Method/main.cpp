#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <cmath>
#include "initialize_matrix.hpp"
#include "algorithm.hpp"

using namespace std;
/*
int main(int argc, char* argv[]) {
    
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " n m r s [filename]" << endl;
        return 1;
    }
    
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int r = atoi(argv[3]);
    int s = atoi(argv[4]);
    int result = 0;
    //int end_of_inverse = 0;
    double t1 = 0, t2 = 0, r1 = -1, r2 = -1;

    double* A = new double[n * n];
    double* B = new double[n * n];
    double* U = new double[(m+1)*(m+1)];
    double* ProductResult = new double[m*m];

    //Зануляем вспомогательную матрицу
    for (int i = 0; i < (m+1)*(m+1); i++)
    {
        U[i] = 0;
    }
    for (int i = 0; i < m*m; i++)
    {
        ProductResult[i] = 0;
    }
    
    
    if(m <= 0 || s < 0 || s > 4 || n <= 0)
    {
        cerr << "Usage: " << argv[0] << " n m r s [filename]" << endl;
        return 1;
    }

    //Инициализируем матрицы А и В
    if (s == 0) {
        string filename = argv[5];
        result = ReadMatrixFromFile(filename, A, n, m);
    } else {
        FormulaMatrixInitialization(A, n, m, s);
    }
    BuildE(B, n, m);

    //Проверка входного файла
    if (result == -1)
    {
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m);
        delete[] A;
        delete[] B;
        delete[] U;
        delete[] ProductResult;
        return 0;
    }
    PrintMatrix(A, n, m, r);
    PrintMatrix(B, n, m, r);
    double* X = new double[5*4];
    double* Z = new double[4*6];
    X[0] = 1; X[1] = 11; X[2] = 111 ; X[3] = 1111;
    X[4] = 2; X[5] = 22; X[6] = 222; X[7] = 2222;
    X[8] = 3; X[9] = 33; X[10] = 333; X[11] = 3333;
    X[12] = 4; X[13] = 44; X[14] = 444; X[15] = 4444;
    X[16] = 5; X[17] = 55; X[18] = 555; X[19] = 5555;

    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout<<X[i*4+j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl;
    Z[0] = 9; Z[1] = 8; Z[2] = 7; Z[3] = 6; Z[4] = 5; Z[5] = 4;
    Z[6] = 3; Z[7] = 2; Z[8] = 1; Z[9] = 0; Z[10] = -1; Z[11] = -2;
    Z[12] = -3; Z[13] = -4; Z[14] = -5; Z[15] = -6; Z[16] = -7; Z[17] = -8;
    Z[18] = -9; Z[19] = -10; Z[20] = -11; Z[21] = -12; Z[22] = -13; Z[23] = -14; 

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            cout<<Z[i*6+j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl;


    BlockMulOptimezed(X, Z, ProductResult, 5, 4, 6);
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            cout<<ProductResult[i*6+j]<<"  ";
        }
        cout<<endl;
    }

    
    cout<<endl;
    BlockMul(X, Z, ProductResult, 5,4,6);
    
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            cout<<ProductResult[i*6+j]<<"  ";
        }
        cout<<endl;
    }

    
    //Вычисление обратной матрицы
    clock_t start1 = clock();

    //end_of_inverse = InverseMatrix(A, B, U, ProductResult, n, m);

    clock_t end1 = clock();
    
    t1 = static_cast<double>(end1 - start1) / CLOCKS_PER_SEC;
    //Инициализируем матрицу для подсчета невязки
    
    if (s == 0) {
        string filename = argv[5];
        result = ReadMatrixFromFile(filename, A, n, m);
    } else {
        FormulaMatrixInitialization(A, n, m, s);
    }

    if (end_of_inverse != 0)
    {    
        PrintMatrix(A, n, m, r, true, false);
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m);
        delete[] A;
        delete[] B;
        delete[] U;
        delete[] ProductResult;
        return 0;
    }

    //Печать ответа
    PrintMatrix(A, n, m, r, true, false);
    PrintMatrix(B, n, m, r, true, false);
    
    //Подсчет невязки
    clock_t start2 = clock();

    //r1 = Discrepancy(A, B, n, m);
    //r2 = Discrepancy(B, A, n, m);
     
    clock_t end2 = clock();

    t2 = static_cast<double>(end2 - start2) / CLOCKS_PER_SEC;
    
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m);

    delete[] A;
    delete[] B;
    delete[] U;
    delete[] ProductResult;
    return 0;
}
*/

int main(int argc, char* argv[]) {
    
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " n m r s [filename]" << endl;
        return 1;
    }
    
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int r = atoi(argv[3]);
    int s = atoi(argv[4]);
    int result = 0;
    int end_of_inverse = 0;
    double t1 = 0, t2 = 0, r1 = -1, r2 = -1;

    double* A = new double[n * n];
    double* B = new double[n * n];
    double* U = new double[(m+1)*(m+1)];
    double* ProductResult = new double[m*m];
    double* results = new double[m];
    double* Sum = new double[m*m];
    double* Column = new double[m];

    //Зануляем вспомогательную матрицу
    for (int i = 0; i < (m+1)*(m+1); i++)
    {
        U[i] = 0;
    }
    for (int i = 0; i < m*m; i++)
    {
        ProductResult[i] = 0;
    }
    for (int i = 0; i < m; i++)
    {
        results[i] = 0;
    }
    
    if(m <= 0 || s < 0 || s > 4 || n <= 0)
    {
        cerr << "Usage: " << argv[0] << " n m r s [filename]" << endl;
        return 1;
    }

    //Инициализируем матрицы А и В
    if (s == 0) {
        string filename = argv[5];
        result = ReadMatrixFromFile(filename, A, n, m);
    } else {
        FormulaMatrixInitialization(A, n, m, s);
    }
    BuildE(B, n, m);

    //Проверка входного файла
    if (result == -1)
    {
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m);
        delete[] A;
        delete[] B;
        delete[] U;
        delete[] ProductResult;
        delete[] results;
        delete[] Sum;
        delete[] Column;
        return 0;
    }
        
    double norm = Norm(A, results, n, m);

    //Вычисление обратной матрицы
    clock_t start1 = clock();

    end_of_inverse = InverseMatrix(A, B, U, ProductResult, norm, n, m);

    clock_t end1 = clock();
    
    t1 = static_cast<double>(end1 - start1) / CLOCKS_PER_SEC;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    //Инициализируем матрицу для подсчета невязки
    
    if (s == 0) {
        string filename = argv[5];
        result = ReadMatrixFromFile(filename, A, n, m);
    } else {
        FormulaMatrixInitialization(A, n, m, s);
    }

    if (end_of_inverse != 0)
    {    
        PrintMatrix(A, n, m, r, true, false);
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m);
        delete[] A;
        delete[] B;
        delete[] U;
        delete[] ProductResult;
        delete[] results;
        delete[] Sum;
        delete[] Column;
        return 0;
    }

    //Печать ответа
    PrintMatrix(A, n, m, r, true, false);
    PrintMatrix(B, n, m, r, true, false);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Sum[i*m + j] = 0;
            ProductResult[i*m +j] = 0;
        }  
    }  
     for (int j = 0; j < m; j++)
    {
        Column[j] = 0;
    }     

    //Подсчет невязки
    clock_t start2 = clock();
    if (n <= 11000)
    {
        r1 = Discrepancy(A, B, Column, ProductResult, Sum, n, m);
        r2 = Discrepancy(B, A, Column, ProductResult, Sum, n, m);
    }
    else
    {
        r1 = 0;
        r2 = 0;
    }
    clock_t end2 = clock();

    t2 = static_cast<double>(end2 - start2) / CLOCKS_PER_SEC;
    
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 24, r1, r2, t1, t2, s, n, m);

    delete[] A;
    delete[] B;
    delete[] U;
    delete[] ProductResult;
    delete[] results;
    delete[] Sum;
    delete[] Column;
    return 0;
}

/*
int main()
{
    int n = 100;
    int m = 3;
    int r = n;
    int s = 2;
    double t1 = 0, t2 = 0, r1 = -1, r2 = -1;

    double* A = new double[n * n];
    double* B = new double[n * n];
    double* AA = new double[n * n];
    double* BB = new double[n * n];
    double* product_result = new double[m * m];
    double* C = new double[n * n];
    double* D = new double[n * n];

    double* U = new double[(m+1)*(m+1)];
    int result = 0;

    
    FormulaMatrixInitialization(A, n, m, s);
    FormulaMatrixInitialization(AA, n, 1, s);
    //FormulaMatrixInitialization(B, n, 1, 1);
    //double B[6] = {1, 2, 3, 4, 5, 6};
    //result = ReadMatrixFromFile("input.txt", A, n, m);
    //result = ReadMatrixFromFile("input.txt", AA, n, 1);
    if (result == -1)
    {
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", "a.out", 24, r1, r2, t1, t2, s, n, m);
        delete[] A;
        delete[] B;
        delete[] U;
        delete[] BB;
        delete[] AA;
        delete[] product_result;
        delete[] C;
        return 0;
    }
    
    //ReadMatrixFromFile("inverse.txt", B, n, m);

    //BuildE(A, n, m);

    //BuildE(D, n, m);
    BlockMul(AA, B,C, n,n,2);
    PrintMatrix(AA, n, 1, r, false);
    PrintMatrix(B, n, 1, r, false);
    PrintMatrix(C, n, 1, r, false);

    
    
    cout<<"Initialization"<<endl;
    BuildE(B, n, m);
    BuildE(BB, n, m);
    PrintMatrix(A,n, m, r);
    PrintMatrix(B, n, m, r);
    //PrintMatrix(AA,n, 1, r);
    PrintMatrix(BB,n, 1, r);
    
    
    PrintMatrix(C, n, m, r);
    PrintMatrix(D, n, m, r);

    
    cout<<"Triungulize"<<endl;
    if(Triungulize(A, U, n, n) == 0)
    {
        PrintMatrix(A, n, m, r);
        PrintMatrix(U, m, m, r);
        ApplyMatrix(U, B, n, n);
        PrintMatrix(B, n, m, r, false);
        
        cout<<"Second Step"<<endl;
        ZeroOut(A, C, U, m, m);
        PrintMatrix(A, n, m, r);
        PrintMatrix(C, n, m, r);
        PrintMatrix(U, m+1, m+1, r+1);
        ApplyMatrixToPair(U, B, D, n, n);
        PrintMatrix(B, n, m, r);
        PrintMatrix(D, n, m, r);
    }
    clock_t start1 = clock();
    
    if(InverseMatrix(A, B, U, product_result, n, m) == 0)
    {
        PrintMatrix(A,n, m, r);
        PrintMatrix(B, n, m, r);
    }
    
    U = new double[(n+1)*(n+1)];
    cout<<"ANSWER";
    Triungulize(AA, U, n, n);
    ApplyMatrix(U, BB, n, n, n);
    InverseTriungleBlock(AA, BB, n, false);

    PrintMatrix(AA, n, 1, r);
    PrintMatrix(A, n, m, r);
    //PrintMatrix(BB, n, 1, r);
    PrintMatrix(B, n, m, r);
    
    clock_t start1 = clock();
    if(InverseMatrix(A, B, U, product_result, n, m) == 0)
    {
        PrintMatrix(A,n, m, r);
        PrintMatrix(B, n, m, r);
    }
    
    
    Triungulize(AA, U, n, n);
    ApplyMatrix(U, BB, n, n, n);

    
    ReplaceWith(BB, AA, n, n);
    PrintMatrix(AA, n, 1, r);
    PrintMatrix(BB, n, 1, r);
    //PrintMatrix(U, n, 1, n);
   
    //result = ReadMatrixFromFile("input.txt", AA, n, 1);    
    clock_t end1 = clock();
    
    t1 = static_cast<double>(end1 - start1) / CLOCKS_PER_SEC;

    clock_t start2 = clock();

    //r1 = CalcError(A, B, n, m);
    FormulaMatrixInitialization(A, n, m, s);

    r1 = Discrepancy(A, B, n, m);
    //FormulaMatrixInitialization(A, n, m, s);

    //r2 = CalcError(B, A, n, m);
    r2 = Discrepancy(B, A, n, m);
     
    clock_t end2 = clock();
    
    t2 = static_cast<double>(end2 - start2) / CLOCKS_PER_SEC;
    
    cout<<endl;
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", "a.out", 24, r1, r2, t1, t2, s, n, m);

    delete[] A;
    delete[] B;
    delete[] U;
    delete[] BB;
    delete[] AA;
    delete[] product_result;
    delete[] C;


    return 0;
}*/