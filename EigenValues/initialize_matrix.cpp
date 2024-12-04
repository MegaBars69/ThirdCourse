#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <fstream>

#define EPSILON pow(10, -15)

int ReadMatrixFromFile(char* filename, double* A, int n)
{
    int i = 0;
    double el;
    char ch;
    int el_in_col = 0;
 
    FILE *fp = fopen(filename, "r");
    if(fp == nullptr)
    {
        printf("Can't open file \n");

        return 1;
    }
    while (true) 
    {
        int res = fscanf(fp, "%lf", &el);
        if(res == EOF)
        {  
            if(i != n*n)
            {
                printf("To few elements in file(<n)\n");
                fclose(fp);
                return 3;
            }
            break;
        }
        else if (res== 0)
        {
            printf("Problem with reading an element \n");
            fclose(fp);
            return 2;
        }
        

        if (i >= n*n)
        {
            printf("To many elements in file(>n)\n");
            fclose(fp);
            return 4;
        }
        else
        {
            A[i] = el;
            el_in_col++;
        }
        
        ch = fgetc(fp);
        if (ch == '\n' && (el_in_col < n) && (el_in_col > n)) 
        {
            printf("Error: Invalide file format \n");
            return 5;
        }
        else if (ch == '\n' && el_in_col == n)
        {
            el_in_col == 0;
        }
        else
        {
            ungetc(ch, fp);
        }
        i++;
    }
    fclose(fp);
    return 0;
}

double f(int i, int j, int n, int s)
{
    switch (s) {
    case 1:
        return  n - ((i < j) ? j : i);
        break;
    case 2:
        if (i == j)
        {
            return 2;
        }
        else if (abs(i - j) == 1)
        {
            return -1;
        }
        else
        {
            return 0;
        }  
        break;
    case 3:
        if (i == j && j < n-1)
        {
            return 1;
        }
        else if (j == n-1)
        {
            return i;
        }
        else if(i == n-1)
        {
            return j;
        }
        else
        {
            return 0;
        } 
        
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

void FormulaMatrixInitialization(double* A, int n, int s) 
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i*n + j] = f(i, j, n, s);
        }    
    }
    
}
/*
void FormulaMatrixInitialization(double* A, int n, int s) 
{
    int p = 1, m = 1, K = 0;
    int l = n%m;
    int k = (n-l)/m;
    int i,j;
    int col_block_size, row_block_size;


    for (int bi = K; bi < k+1; bi += p)
    {
        col_block_size = (bi<k ? m : l);
        for (int bj = 0; bj < k+1; bj++)
        {
            row_block_size = (bj<k ? m : l);
            double* pa = (A + bi * m*n + bj*m*col_block_size);

            for (int i_ = 0; i_ < col_block_size; i_++)
            {
                for (int j_ = 0; j_ < row_block_size; j_++, pa++)
                {
                    i = bi*m + i_;
                    j = bj*m + j_;
                    *(pa) = f(i, j, n, s);
                }
                
            }    
        }    
    }
}
*/

int min(int a, int b)
{
    return (a > b ? b : a);
}

void PrintMatrix(double* matrix, int n, int m) 
{
    printf("\n");
    for (int i = 0; i < min(m, n); i++) {
        for (int j = 0; j < min(m, n); j++) {
            printf("%10.3e ", *(matrix + i * n + j)); 
        }
        printf("\n");
    }
    printf("\n");
}

int CheckMatrix(double* A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                if (fabs(A[i*n + j] - A[j*n + i]) > EPSILON)
                {
                    return  1;
                }
                
            }
            
        }
        
    }
    
    return 0;
}
