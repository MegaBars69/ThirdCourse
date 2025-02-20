prog: main.o MatMul.o
	g++ main.o MatMul.o -o prog
main.o : main.cpp
	g++ -g -c main.cpp -o main.o
MatMul.o: MatMul.cpp
	g++ -g -c MatMul.cpp -o MatMul.o   