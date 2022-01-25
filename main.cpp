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
 * @param t number of threads
 */
int jacobi(float* u_old, float* u_new, int N, int M, float h, int t) {

  float h2 = h*h;
  float * erro = new float[t];
  int cont;
  float tol = 1;

  // approximates the result until a tolerance is satisfied
  for (cont = 0; tol > .000001; cont += 2) {

    // cout << cont << endl;
    memset(erro, 0.0, t * sizeof(float));

    // creating parallel region
    #pragma omp parallel
    {
      int id = omp_get_thread_num();

      // calculates u_new from u_old
      #pragma omp for collapse(1) schedule (OPENMP_SCHEDULE)
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

      // calculates u_old from u_new
      #pragma omp for collapse(1) schedule (OPENMP_SCHEDULE)
      for (int j = 1; j < M - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
          u_old[ j*N + i ] = ( 
            h2 * S (i*h, j*h) +
            u_new[ j*N + (i - 1) ] +
            u_new[ j*N + (i + 1) ] +
            u_new[ (j - 1)*N + i ] +
            u_new[ (j + 1)*N + i ] 
          )/4;
          // erro[id] += abs(u_new[ j*N + i ] - u_old[ j*N + i ]);
        }
      }

      // calculating the error separetly providaded a better performance
      #pragma omp for collapse(1) schedule (OPENMP_SCHEDULE)
      for (int j = 1; j < M - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
          erro[id] += abs(u_new[ j*N + i ] - u_old[ j*N + i ]);
        }
      }
    }

    // calculating the tolerance based on last two iterations
    for (int i = 1; i < t; i++ ){
      erro[0] += erro[i];
    }
    
    tol = erro[0] / ( (N - 2) * (M - 2) );

  }

  return cont;

}

/**
 * @brief Calculates next iteration of the linear system by the SOR Back-Red method
 * 
 * @param u_old current solution
 * @param u_next next iteration solution
 * @param N number of columns
 * @param M number of rows
 * @param h grid width
 * @param w acceleration factor
 * @return the tolerance
 */
float SOR_BR(float* u_old, float* u_new, int N, int M, float h) {

  float h2 = h*h;
  float erro = 0.0;

  // calculate next approximation for the solution

  // iterate through all elements such that i + j is odd
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
      erro += abs(u_new[ j*N + i ] - u_old[ j*N + i ]);
    }
  }

  return erro / ( (N - 2) * (M - 2) );

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

  int N = 400;          // number of columns
  int M = N;            // number of rows
  float h = 1.0/N;      // grid width
  int cont;             // iterations counter

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

  cont = jacobi(u_old, u_new, N, M, h, t);

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
