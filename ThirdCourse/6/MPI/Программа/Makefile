CC=mpicxx
CFLAGS= mpicxx -isystem /opt/impi-5.1.3.223/intel64/include -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
a.out: main.o initialize_matrix.o algorithm.o
	$(CC)  main.o initialize_matrix.o algorithm.o -o a.out
main.o : main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o 
initialize_matrix.o: initialize_matrix.cpp
	$(CC) $(CFLAGS) -c initialize_matrix.cpp -o initialize_matrix.o 
algorithm.o : algorithm.cpp
	$(CC) $(CFLAGS) -c algorithm.cpp -o algorithm.o
clean : 
	rm -f *.o a.out