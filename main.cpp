#include <stdio.h>
#include <chrono>
#include <iostream>
#include <cstring>
#include <fstream>
#include "omp.h"

#define OPENMP_SCHEDULE dynamic  

using namespace std;

// right side equation
float S(float x, float y){
  return -10*(x*x + y*y + 5);
}

/**
 * @brief Calculates next iteration of the linear system by the jacobi method
 * 
 * @param u_old current solution
 * @param u_next next iteration solution
 * @param N number of columns
 * @param M number of rows
 * @param h grid width
 * @return true When the tolerance is satisfied
 */
void jacobi(float* u_old, float* u_new, int N, int M, float h) {

  float h2 = h*h;

  // calculate next approximation for the solution, which is basically the mean of it's neighbours
  // note that it doesn't iterate through the boundary
  #pragma omp for collapse(1) // schedule (OPENMP_SCHEDULE)
  for (int j = 1; j < M - 1; j++) {
    for (int i = 1; i < N - 1; i++) {
      u_new[ j*N + i ] = ( 
        h2 * S (i*h, j*h) +
        u_old[ j*N + (i - 1) ] +
        u_old[ j*N + (i + 1) ] +
        u_old[ (j - 1)*N + i ] +
        u_old[ (j + 1)*N + i ] 
      )/4;
    }
  }
}

void saveDataBin(float* u, int* N, int* M, float* h) {

  FILE* file;

  file = fopen ("output.bin", "wb");

  fwrite(&h, sizeof(float), 4, file);
  fwrite(&N, sizeof(int), 4, file);
  fwrite(&M, sizeof(int), 4, file);
  fwrite(u, sizeof(float), sizeof(u), file);
  
  fclose(file);

}

void saveDataASCII(float* u, int N, int M, float h) {

  ofstream file;

  file.open("output.dat");

  file << h << "\n";
  file << N << "\n";
  file << M << "\n";

  for (int i = 0; i < N * M; i++) {
    file << u[i] << " ";
    if ((i + 1) % N == 0){
      file << "\n";
    }
  }

}

int main(){
  
  // Domain:
  //    [0,1] x [0,1]
  //
  // BC:
  //    u(x,0) = 0
  //    u(x,1) = 1
  //    u(0,y) = 0
  //    u(1,y) = 0

  int N = 400;           // number of columns
  int M = N;            // number of rows
  float h = 1.0/N;      // grid width

  int cont = 0;
  float* sum;
  float total = 0;
  float tol = 10;

  int size = N * M;     // size of vector

  float* u_old = new float[size];
  float* u_new = new float[size];

  // OpenMP config
  int t = 1;

#ifdef _OPENMP
  omp_set_num_threads(8);
  #pragma omp parallel
  t = omp_get_num_threads();
#endif

  cout << "Threads: " << t << endl;

  sum = new float[t];

  // setting grid values to 0 initially
  memset(u_old, 0.0, size * sizeof(float));
  memset(u_new, 0.0, size * sizeof(float));

  // boundary condition at u(x,1)
  for (int i = size - N; i < size; i++) {
    u_old[i] = 1.0;
    u_new[i] = 1.0;
  }

  // solving by the jacobi method (exercise 1)
  cout << "Executing jacobi...\n";
  auto inicio = chrono::high_resolution_clock::now();

  // approximates the result until a tolerance is satisfied
  for (cont = 0; tol > .000001; cont += 2) {

    // reseting values
    memset(sum, 0.0, t * sizeof(float));
    total = 0;

    // creating parallel region
    #pragma omp parallel 
    {
      int id = omp_get_thread_num();
    
      // calculates the two next iterations
      jacobi(u_old, u_new, N, M, h);
      jacobi(u_new, u_old, N, M, h);

      #pragma omp for collapse(1)
      for (int j = 1; j < M - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
          sum[id] += abs(u_new[j*N + i] - u_old[j*N + i]);
        }
      }
    }

    // calculates tolerance
    for (int i = 0; i < t; i++ ){
      total += sum[i];
    }
    tol = (total/((N - 2.0) * (M - 2.0)));

  }

  auto final = chrono::high_resolution_clock::now();
  chrono::duration<double> intervalo = final - inicio;

  cout << "\nTime spent (jacobi): " << intervalo.count() << "s\n";
  cout << "Iterations: " << cont << "\n";
  cout << "u(0.5, 0.5): " << u_old[ size/2 + N/2 ] << "\n\n";

  // saving data for plotting
  cout << "Saving data...\n";
  saveDataASCII(u_old, N, M, h);
  saveDataBin(u_old, &N, &M, &h);
  cout << "Data saved succesfully!\n";

  // solving by the SOR red-black method (exercise 2)

  return 0;
}
