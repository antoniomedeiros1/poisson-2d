make clean
make

export OMP_NUM_THREADS=1
./main 0.0000001

export OMP_NUM_THREADS=2
./main 0.0000001

export OMP_NUM_THREADS=4
./main 0.0000001

export OMP_NUM_THREADS=8
./main 0.0000001