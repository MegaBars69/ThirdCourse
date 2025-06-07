#include <mpi.h>
#include <cmath>
#include <iostream>

void upd_array(double *a, int n, int p, int k, MPI_Comm comm, MPI_Status *st);

void upd_array(double *a, int n, int p, int k, MPI_Comm comm, MPI_Status *st)
{
    int ln = n/p;
    int l = n%p;

    int count = 0;
    double sum = 0;
    
    int global_count = 0;
    double global_sum = 0;
    int s = ln*k, end = ln*(k+1);
    double har = 0;
    if(k == p-1)
    {
        end = n;
    }
    if(n <= 3)
    {
        return;
    }
    for (int i = s; i < end; i++)
    {
        if(i > 0 && i + 2 < n)
        {
            har = pow(a[i-1] * a[i]*  a[i+1] * a[i+2],0.25);
            if(a[i] < har)
            {
                sum += a[i];
                count++;
            }
        }
    }

    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&count, &global_count, 1, MPI_INT, MPI_SUM, comm);

    double av = global_sum/global_count;
    double anm1, an, anp1m, anp2;
    for (int i = s; i < end; i++)
    {
        if(i == s)
        {
            anm1 = a[i-1];
            an = a[i];
        }

        if(i > 0 && i + 2 < n)
        {
            
            an = a[i];
            har = pow(anm1*an*a[i+1]*a[i+2],0.25);
            if(a[i] < har)
            {
                a[i] = av;
            }
        }
        anm1 = an;
    }

    for (int pn = 1; pn < p; pn++)
    {
        int buf_size = pn*ln + (pn == p-1 ? l : 0);
        
        if(k == 0)
        {
            MPI_Recv(a + ln*pn, buf_size, MPI_DOUBLE, pn, 0, comm, st);
        }
        else if(pn == k)
        {
            MPI_Send(a + ln*pn, buf_size, MPI_DOUBLE, 0, 0, comm);
        }
    }
    
}

int main(int argc, char* argv[])
{
    int p, k;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status st;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &k);

    if(argc != 3)
    {
        MPI_Finalize(); 
        return 0; 
    }
    int n = atoi(argv[1]);
    double *array = new double [n];

    if(k == 0)
    {
        FILE *file = fopen(argv[2], "r");
        if(!file)
        {
            std::cout<<"Couldn't open file\n";
            delete[] array;
            MPI_Finalize(); 
            return 0; 
        }
        for(int i = 0; i<n;i++)
        {
            if(!fscanf(file, "%lf", &array[i]))
            {
                MPI_Finalize();
                return 0;
            }
        }
        std::cout<<"Array"<<" : ";
        for (int i = 0; i < n; i++)
        {
            std::cout<<array[i]<<' ';
        } 
        printf("\n");
    }
    MPI_Bcast(array, n, MPI_DOUBLE, 0, comm);
    double reduce_sum = 0;
    double time = MPI_Wtime();

    upd_array(array, n, p, k, comm, &st);

    time = MPI_Wtime() - time;

    if(k == 0)
    {
        std::cout<<"UPD Array: ";
        for (int i = 0; i < n; i++)
        {
            std::cout<<array[i]<<' ';
        }
        std::cout<<'\n';
    }
    for (int i = 0; i < p; i++)
    {
        MPI_Barrier(comm);
        if(i == k)
        {
            std::cout<<"Time from "<<k<<" proc: "<<time<<'\n';
        }
        MPI_Barrier(comm);
    }
    MPI_Reduce(&time, &reduce_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    
    if(k == 0)
    {
        std::cout<<"All time "<<reduce_sum<<'\n';
    }

    delete[] array;
    MPI_Finalize();
    return 0;
}