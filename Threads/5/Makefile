CFLAGS= -O0 -g
a.out: main.o func.o
	g++ main.o func.o -o a.out
main.o : main.cpp
	g++ $(CFLAGS) -c main.cpp -o main.o 
func.o: func.cpp
	g++ $(CFLAGS) -c func.cpp -o func.o 
