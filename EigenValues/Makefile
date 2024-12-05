CC=g++
CFLAGS= -O0 -g
a.out: main.o initialize_matrix.o algorithm.o
	g++ -g main.o initialize_matrix.o algorithm.o -o a.out
main.o : main.cpp
	g++ $(CFLAGS) -c main.cpp -o main.o 
initialize_matrix.o: initialize_matrix.cpp
	g++ $(CFLAGS) -c initialize_matrix.cpp -o initialize_matrix.o 
algorithm.o : algorithm.cpp
	g++ $(CFLAGS) -c algorithm.cpp -o algorithm.o 