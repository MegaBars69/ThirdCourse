CFLAGS= -O3 -O0 -mfpmath=sse -fstack-protector-all -lpthread -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
a.out: main.o func.o
	g++ main.o func.o -o a.out
main.o : main.cpp
	g++ $(CFLAGS) -c main.cpp -o main.o 
func.o: func.cpp
	g++ $(CFLAGS) -c func.cpp -o func.o 
