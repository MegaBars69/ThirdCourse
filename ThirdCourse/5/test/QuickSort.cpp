#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void Print(vector<int> A)
{
    cout<<endl<<"[ ";
    for (int i = 1; i < A.size(); i++)
    {
        cout<<A[i]<<", ";
    }
    cout<<" ]"<<endl;
}


int Partion(vector<int>& A, int p, int r)
{

    int x = A[r];
    int i = p-1;
    for (int j = p; j < r; j++)
    {
        if(A[j] <= x)
        {
            i++;
            swap(A[i], A[j]);
        }
    }
    swap(A[i+1], A[r]);
    return i+1;
}

void QuickSort(vector<int>& A, int p, int r)
{
    if(p<r)
    {
        int q = Partion(A, p, r);
        QuickSort(A,p,q-1);
        QuickSort(A,q+1,r);
    }
}

int main(){
    vector<int> A = {-1,2,8,7,1,3,5,6,4,10,0,11,13,15};
    Print(A);
    QuickSort(A, 1, A.size()-1);
    Print(A);
    
    return 0;
}