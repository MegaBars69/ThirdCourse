#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void upd_array(double* a, int n, int p, int k)
{
    int loc_len = n/p;
    int upb = (k+1)*loc_len, lob = k*loc_len;
    if(k == p-1)
    {
        upb = n;
    }
    int count = 0;
    double s = 0, av;
    double beg_s = 0; int beg_count = 0;
    bool ifs = false, end = false;

    if(k > 0)
    {
        while (a[lob -1] < a[lob] && lob<n)
        {
            ifs = true;
            s+=a[lob];
            count++;
            lob++;
        }    
        beg_count = count;
        beg_s = s;
        count = 0;
        s = 0;
    }
    
    for (int i = lob; i < upb; i++)
    {
        if(i != upb-1 && a[i] < a[i+1])
        {
            count++;
            s+=a[i];
        }
        else if (i == upb-1 && k < p-1 && a[i] < a[i+1])
        {
            count++;
            s+=a[i];
            end = true;
        }
        else 
        {
            s+=a[i];
            count++;
            av = s/count;
            for(int i_loc = i; i_loc > i-count; i_loc--)
            {
                a[i_loc] = av;
            } 
            count = 0;
            s = 0;
        }
    }

    if(end)
    {
        while (a[upb] < a[upb+1] && upb < n-1)
        {
            count++;
            s+=a[upb];
            upb++;
        }
        s+=a[upb];
        count++;
        av = s/count;
        for (int i = upb; i > upb - count; i--)
        {
            a[i] = av;
        }
    }
    if(ifs && k>0)
    {
        lob = k*loc_len-1;

        while (lob>=0 && a[lob] <a[lob+1])
        {
            beg_s+=a[lob];
            beg_count++;
            lob--;
        }
        lob++;
        av = beg_s/beg_count;
        for (int i = lob; i < lob+beg_count; i++)
        {
            a[i] = av;
        }
    }
}

int main(int argc, char *argv[]) 
{
    int k, p;
    int n, loc_len, l;
    double time = 0;
    double *array = nullptr, *buf = nullptr;

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status st;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &k);
    MPI_Comm_size(comm, &p);

    // Проверка аргументов командной строки только на процессе 0
    if (k == 0) 
    {
        if (argc != 3) {
            fprintf(stderr, "Usage: %s <n> <filename>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        n = atoi(argv[1]);
        if (n <= 0) {
            fprintf(stderr, "Invalid array length n: %d\n", n);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    loc_len = n/p;
    l = n%p;
    array = new double[n];

    if(k==0)
    {
        FILE *file = fopen(argv[2], "r");
        if (!file) {
            perror("Failed to open file");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (!array) {
            fprintf(stderr, "Memory allocation failed\n");
            fclose(file);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < n; i++) {
            if (fscanf(file, "%lf", &array[i]) != 1) {
                fprintf(stderr, "Failed to read element %d from file\n", i);
                delete[] array;
                fclose(file);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
    
    MPI_Bcast(array, n, MPI_DOUBLE, 0, comm);

    if(k == 0)
    {
        // Каждый процесс выводит принятый массив
        printf("\nProcess %d rec array: ", k);
        for (int i = 0; i < n; i++) {
            printf("%.0f ", array[i]);
        }
        printf("\n");
    }

    time = MPI_Wtime();
    upd_array(array, n, p, k);
    time = MPI_Wtime() - time;

    for (int pn = p-1; pn > 0; pn--)
    {
        int ind = pn*loc_len;
        int amount = loc_len + (pn == p-1 ? l : 0);
        if(k == 0)
        {
            MPI_Recv(array + ind, amount, MPI_DOUBLE, pn, 0, comm, &st);
        }
        else if(k == pn)
        {
            MPI_Send(array + ind, amount, MPI_DOUBLE, 0, 0, comm);
        }
    }
    if(k == 0)
    {
        printf("Upd Array: ", k);
        for (int i = 0; i < n; i++) {
            printf("%.2f ", array[i]);
        }
        printf("\n");
    }
    delete[] array;

    for (int i = 0; i < p; i++)
    {
        MPI_Barrier(comm);
        if(i == k)
        {
            std::cout<<"Proc "<<k<<" time "<<time<<"\n";
        }
        MPI_Barrier(comm);
    }
    double all_time = 0;
    MPI_Reduce(&time, &all_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    if(k == 0)
    {
        std::cout<<"All time: "<<all_time<<'\n';
    }
    MPI_Finalize();
    return 0;
}
