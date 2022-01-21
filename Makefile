
CC := g++

CFLAGS := -fopenmp -g -O3

all: main.o 
	@$(CC) $(CFLAGS) -o main main.o

main.o: main.cpp 
	@$(CC) $(CFLAGS) -c main.cpp

clean::
	@rm *.o *.bin *.dat
