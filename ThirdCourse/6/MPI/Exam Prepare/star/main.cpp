#include <mpi.h>
#include <iostream>
#include <cmath>

#define EPSILON 1e-15

void upd_array(double* a, int n, int p, int k, MPI_Comm comm);


void upd_array(double* a, int n, int p, int k, MPI_Comm comm)
{
    int loc_len = n/p;
    int low = k*loc_len;
    int upb = (k+1)*loc_len;
    double s = 0;
    int count = 0;
    double prev_q  = 0, q = 0;
    bool in_seq = false;
    bool added = false;

    if(k == p-1)
    {
        upb = n; 
    }
    if(k > 0)
    {
        if(fabs(a[low-1]) > EPSILON && fabs(a[low])> EPSILON && fabs(a[low+1])> EPSILON)
        {
            prev_q = a[low]/a[low-1];
            q = a[low+1]/a[low];
            if(fabs(prev_q - q) < EPSILON)
            {
                s += a[low];
                in_seq = true;
                count++;
                
            }
            prev_q = q;
        }
    }
    if(k > 0 && loc_len > 1)
    {
        if(!in_seq && fabs(a[low-2]) > EPSILON && fabs(a[low-1])> EPSILON && fabs(a[low])> EPSILON)
        {
            prev_q = a[low-1]/a[low-2];
            q = a[low]/a[low-1];
            if(fabs(prev_q - q) < EPSILON)
            {
                s += a[low];
                in_seq = true;
                count++;
                
            }
            prev_q = q;
        }
    }
    if(k < p-1)
    {
        if(upb+1 < n)
        {
            if(fabs(a[upb-1]) > EPSILON && fabs(a[upb])> EPSILON && fabs(a[upb+1])> EPSILON)
            {
                prev_q = a[upb]/a[upb-1];
                q = a[upb+1]/a[upb];
                if(fabs(prev_q - q) < EPSILON)
                {
                    s += a[upb-1];
                    count++;
                    if(k == 0 && upb - low < 3)
                    {
                        s+=a[low];
                        count++;
                    }
                    added = true;
                }
            }
        }
    }

    for (int i = low+1; i < upb; i++)
    {
        if(i+1 < upb && fabs(a[i-1]) > EPSILON && fabs(a[i])> EPSILON && fabs(a[i+1])> EPSILON)
        {
            prev_q = a[i]/a[i-1];
            q = a[i+1]/a[i];
            if(fabs(prev_q - q) < EPSILON)
            {
                if(!in_seq)
                {
                    s += a[i-1];
                    count++; 
                }
                else if(k == 0 && i == low+1)
                {
                    s += a[i-1];
                    count++;
                }
                in_seq = true;
                if(upb - low>2)
                {
                    s += a[i];
                    count++;
                }
                
                
            }
            else
            {
                if(in_seq)
                {
                    s += a[i];
                    count++; 
                }
                in_seq = false;
            }
        }
        else if(i == upb-1 &&in_seq && !added)
        {
            s += a[i];
            count++;
        }
        else if(i == n-1 && in_seq)
        {
            s += a[i];
            count++;
        }
    }
    std::cout<<"In "<<k<<": "<<count<<" "<<s<<"\n";
    double reduce_s = 0;
    int reduce_count = 0;
    MPI_Allreduce(&s, &reduce_s, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&count, &reduce_count, 1, MPI_INT, MPI_SUM, comm);

    double av = reduce_s/reduce_count;
    if(count == 0)
    {
        return;
    }
    
    if(k > 0)
    {
        if(fabs(a[low-1]) > EPSILON && fabs(a[low])> EPSILON && fabs(a[low+1])> EPSILON)
        {
            prev_q = a[low]/a[low-1];
            q = a[low+1]/a[low];
            if(fabs(prev_q - q) < EPSILON)
            {
                a[low] = av;
                if(upb - low < 3)
                {
                    a[low+1] = av;
                }

            }
            prev_q = q;
        }
    }
    
    for (int i = low+1; i < upb; i++)
    {
        if(i< n-1 && fabs(a[i-1]) > EPSILON && fabs(a[i])> EPSILON && fabs(a[i+1])> EPSILON)
        {
            if(k == 0 && i == low+1)
            {
                prev_q = a[i]/a[i-1];
            }
            q = a[i+1]/a[i];
            if(fabs(prev_q - q) < EPSILON)
            {
                if(!in_seq)
                {
                    a[i-1] = av; 
                }
                in_seq = true;
                a[i] = av;
                if(k == 0 && i == low+1)
                {
                    a[i-1] = av;
                }
            }
            else
            {
                if(in_seq)
                {
                    a[i] = av;
                }
                in_seq = false;
            }
            prev_q = q;
        }
        else if(i == n-1 && in_seq)
        {
            a[i] = av;
        }
        
    }
    if(k < p-1)
    {
        if(upb+1 < n)
        {
            if(fabs(a[upb-1]) > EPSILON && fabs(a[upb])> EPSILON && fabs(a[upb+1])> EPSILON)
            {
                prev_q = a[upb]/a[upb-1];
                q = a[upb+1]/a[upb];
                if(fabs(prev_q - q) < EPSILON)
                {
                    a[upb-1] = av;
                }
            }
        }
    }
}


int main(int argc, char* argv[])
{
    int p, k, n;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status st;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &k);
    if(argc != 3)
    {
        MPI_Finalize();
        return 0;
    }

    n = atoi(argv[1]);
    double *array = new double[n];
    if(k == 0)
    {
        FILE *file = fopen(argv[2],"r");
        if(!file)
        {
            std::cout<<"Error opening file\n";
            delete[] array;
            MPI_Finalize();
            return 0;
        }
        for (int i = 0; i < n; i++)
        {
            if(!fscanf(file, "%lf", &array[i]))
            {
                std::cout<<"Error reading file\n";
                delete[] array;
                fclose(file);
                MPI_Finalize();
                return 0;
            } 
        }
        fclose(file);
        printf("Array: ");
        for (int i = 0; i < n; i++)
        {
            std::cout<<array[i]<<" ";
        }
        printf("\n");
    }

    MPI_Bcast(array, n, MPI_DOUBLE, 0, comm);

    double time = MPI_Wtime();
    upd_array(array, n, p, k, comm);
    time = MPI_Wtime() - time; 

    int loc_len = n/p;
    int l = n%p;
    for (int pn = 1; pn < p; pn++)
    {
        int index = pn*loc_len;
        int buf_size = loc_len + (pn < p-1 ? 0 : l);

        if(k == 0)
        {
            MPI_Recv(array + index, buf_size, MPI_DOUBLE, pn, 0, comm, &st);
        }
        else if(k == pn)
        {
            MPI_Send(array + index, buf_size, MPI_DOUBLE, 0, 0, comm);
        }
    }
    if(k == 0)
    {
        printf("UPD Array: ");
        for (int i = 0; i < n; i++)
        {
            std::cout<<array[i]<<" ";
        }
        printf("\n");
    }
    double reduce_time = 0;
    MPI_Reduce(&time, &reduce_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    /*
    for (int i = 0; i < p; i++)
    {
        MPI_Barrier(comm);
        if(i == k)
        {
            std::cout<<"Time from "<<k<<": "<<time<<"\n";
        }        
        MPI_Barrier(comm);
    }
    if(k == 0)
    {
        std::cout<<"Reduced time: "<<reduce_time<<"\n";
    }
    */

    delete[] array;
    MPI_Finalize();
    return 0;
}