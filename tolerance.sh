make clean
make

export OMP_NUM_THREADS=8

./main 0.01
./main 0.001
./main 0.0001
./main 0.00001
./main 0.000001
./main 0.0000001