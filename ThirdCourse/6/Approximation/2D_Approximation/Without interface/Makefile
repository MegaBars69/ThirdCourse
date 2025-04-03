CC=g++
CFLAGS= -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -ffast-math -Wno-suggest-attribute=format
a.out: main.o initialize_matrix.o algorithm.o thread_function.o
	g++  main.o initialize_matrix.o algorithm.o thread_function -o a.out
main.o : main.cpp
	g++ $(CFLAGS) -c main.cpp -o main.o 
initialize_matrix.o: initialize_matrix.cpp
	g++ $(CFLAGS) -c initialize_matrix.cpp -o initialize_matrix.o 
algorithm.o : algorithm.cpp
	g++ $(CFLAGS) -c algorithm.cpp -o algorithm.o 
thread_function.o : thread_function.cpp
	g++ $(CFLAGS) -c thread_function.cpp -o thread_function.o 
