CC=mpicxx
CFLAGS= -isystem /opt/impi-5.1.3.223/intel64/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -O0 -g
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