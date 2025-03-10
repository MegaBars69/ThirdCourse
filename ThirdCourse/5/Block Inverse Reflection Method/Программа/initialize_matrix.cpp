#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include "initialize_matrix.hpp"

using namespace std;
// Function to build Matrix B (E)

// Function to read matrix from file and write it in a block form

//Доступ на элемент a_{11}^(i,j) блоку A_{i,j} достигается обращением к указателю ∗(a+ i ∗ (k ∗m∗m+m∗ l) + j ∗ m∗block_size_row)

/*
Function to initialize matrix using formula
a_{i,j} = In Block_{i//m, j//m} in {i%col_size, j%row_size} position
i = bi*m + i'
j = bj*m + j'
*/

int ReadMatrixFromFile(const string& filename, double* A,  int n, int m) 
{
    double el;
    int l = n%m;
    int k = (n-l)/m;
    int index, col = 0;
    string line;
    ifstream file(filename);
    if (!file) {
        cerr << "Ошибка: не удалось открыть файл" << std::endl;
        return -1;
    }
    for (int i = 0; i < k; i++)
    {
        for (int bi = 0; bi < m; bi++)
        {
            col = 0;
            if(!getline(file, line))
            {
                cerr<<"Error: invalid file format"<<endl;
                return -1;
            }
            istringstream iss(line);
            for (int bj = 0; bj < k; bj++)
            {
                for (int p = 0; p < m; p++)
                {
                    index = i*n*m + bi*m+bj*m*m+p;
                    if(!(iss >> A[index]))
                    {
                        cerr<<"Error: invalid file format"<<endl;
                        return -1;
                    }
                    col++;
                }
                
            }
            
            for (int lj = 0; lj < l; lj++)
            {
                index = i*n*m + bi*l+k*m*m+lj;
                if(!(iss >> A[index]))
                    {
                    cerr<<"Error: invalid file format"<<endl;
                    return -1;
                }
                col++;
            }
            if(col != n)
            {
                cerr<<"Error: invalid file format (length of line < n)"<<endl;
                return -1;
            }
        }
    }
    for (int li = 0; li < l; li++)
    {
        for (int lj = 0; lj < k; lj++)
        {
            for (int lp = 0; lp < m; lp++)
            {      
                index = k*n*m + li*m + lj*l*m + lp;
                if(!(file >> A[index]))
                    {
                    cerr<<"Error: invalid file format"<<endl;
                    return -1;
                }
            }

        }
        for (int lp = 0; lp < l; lp++)
        {
            index = k*n*m + k*m*l + li*l + lp;
            if(!(file >> A[index]))
            {
                cerr<<"Error: invalid file format"<<endl;
                return -1;
            }
        }
    }
    if(file >> el)
    {
        cerr<<"n*n < amount of elements"<<endl;
        return -1;
    }
    return 0;
}

void BuildE(double* B, int n, int m) {
    int l = n%m;
    int k = (n-l)/m;
    int max_col, max_row;

    for (int bi = 0; bi < k+1; bi++)
    {
        max_col = (bi<k ? m : l);
        for (int bj = 0; bj < k+1; bj++)
        {
            max_row = (bj<k ? m : l);
            double* pb = (B + bi * m*n + bj*m*max_col);
            for (int i = 0; i < max_col; i++)
            {
                for (int j = 0; j< max_row; j++, pb++)
                {
                    *(pb) = (!(bi == bj && i == j) ? 0 : 1); 
                }
                
            }
        }
        
    }   
}

void FormulaMatrixInitialization(double* A, int n, int m, int s) {
    int l = n%m;
    int k = (n-l)/m;
    int i,j;
    int col_block_size, row_block_size;

    switch (s) {
    case 1:  
        //a_{ij} = n - max{i,j}+1       
        for (int bi = 0; bi < k+1; bi++)
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
                        *(pa) = n - ((i < j) ? j : i);
                    }
                    
                }    
            }    
        }
        break;
    case 2:
        //max{i,j} 
        for (int bi = 0; bi < k+1; bi++)
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
                        *(pa) = 1 + ((i >= j) ? i : j);
                    }
                    
                }    
            }    
        }
        break;
    case 3:
        // |i-j|
        for (int bi = 0; bi < k+1; bi++)
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
                        *(pa) = (i - j >= 0 ? (i - j) : (j - i));
                    }
                    
                }    
            }    
        }
        break;
    case 4:
        // 1/(i+j-1)
        for (int bi = 0; bi < k+1; bi++)
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
                        *(pa) = 1/double(i+j+1);
                    }
                    
                }    
            }    
        }
        break;
    default:
        cout<<"Неправильный формат ввода s"<<endl;
        break;
    }
}

