#include <stdio.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include "omp.h"

#define OPENMP_SCHEDULE dynamic  

using namespace std;

// right side equation
float S(float x, float y){
  return 10*(x*x + y*y + 5);
}

/**
 * @brief Calculates an aproximation of the solution of a linear system by the jacobi method
 * 
 * @param tol tolerance for the aproximation
 * @param u_old current solution
 * @param u_next next iteration solution
 * @param N number of columns
 * @param M number of rows
 * @param h grid width
 * @param t number of threads
 * @return number of iterations
 */
int jacobi(float tol, float* u_old, float* u_new, int N, int M, float h, int t) {

  float h2 = h*h;
  float * erro = new float[t];
  int cont;
  float norm = 1.0f;

  // approximates the result until a tolerance is satisfied
  for (cont = 0; norm > tol; cont += 2) {

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
    
    norm = erro[0] / ( (N - 2) * (M - 2) );

  }

  return cont;

}

void SORaux(float* u_old, float* u_new, int N, int M, float h, float w) {

  float h2 = h*h;

  // first pass for i + j odd
  #pragma omp for collapse(1) schedule (OPENMP_SCHEDULE)
  for (int j = 1; j < M - 2; j += 2) {

    for (int i = 2; i < N - 1; i += 2) {
      
      u_new[ j*N + i ] = 
        (1 - w) * u_old[ j*N + i ] + (w * .25 * ( 
        h2 * S (i*h, j*h) +
        u_old[ j*N + (i - 1) ] +
        u_old[ j*N + (i + 1) ] +
        u_old[ (j - 1)*N + i ] +
        u_old[ (j + 1)*N + i ] 
      ));
      
      u_new[ (j + 1)*N + (i - 1) ] = 
        (1 - w) * u_old[ (j + 1)*N + (i - 1) ] + (w * .25 * ( 
        h2 * S (i*h, j*h) +
        u_old[ (j + 1)*N + (i - 2) ] +
        u_old[ (j + 1)*N + (i    ) ] +
        u_old[ (j    )*N + (i - 1) ] +
        u_old[ (j + 2)*N + (i - 1) ] 
      ));

    }

  }

  // second pass for i + j even
  #pragma omp for collapse(1) schedule (OPENMP_SCHEDULE)
  for (int j = 1; j < M - 2; j += 2) {

    for (int i = 1; i < N - 2; i += 2) {
      
      u_new[ j*N + i ] = 
        (1 - w) * u_old[ j*N + i ] + (w * .25 * ( 
        h2 * S (i*h, j*h) +
        u_new[ j*N + (i - 1) ] +
        u_new[ j*N + (i + 1) ] +
        u_new[ (j - 1)*N + i ] +
        u_new[ (j + 1)*N + i ] 
      ));
      
      u_new[ (j + 1)*N + (i + 1) ] = 
        (1 - w) * u_old[ (j + 1)*N + (i + 1) ] + (w * .25 * ( 
        h2 * S (i*h, j*h) +
        u_new[ (j + 1)*N + (i    ) ] +
        u_new[ (j + 1)*N + (i + 2) ] +
        u_new[ (j    )*N + (i + 1) ] +
        u_new[ (j + 2)*N + (i + 1) ] 
      ));

    }

  }
  
}

/**
 * @brief Calculates an aproximation of the solution of a linear system by the SOR Black-Red method
 * 
 * @param tol tolerance for the aproximation
 * @param u_old current solution
 * @param u_next next iteration solution
 * @param N number of columns
 * @param M number of rows
 * @param h grid width
 * @param t number of threads
 * @param w acceleration factor
 * @return number of iterations
 */
int SOR(float tol, float* u_old, float* u_new, int N, int M, float h, int t, float w) {

  float * erro = new float[t];
  int cont;
  float norm = 1.0f;

  // approximates the result until a tolerance is satisfied
  for (cont = 0; norm > tol; cont += 2) {

    // cout << cont << endl;
    memset(erro, 0.0, t * sizeof(float));

    // creating parallel region
    #pragma omp parallel
    {
      int id = omp_get_thread_num();

      SORaux(u_old, u_new, N, M, h, w);
      SORaux(u_new, u_old, N, M, h, w);

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
    
    norm = erro[0] / ( (N - 2) * (M - 2) );

  }

  return cont;

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

void saveDataASCII(string fileName, float* u, int N, int M, float h) {

  ofstream file;

  file.open(fileName + ".dat");

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

int main(int argc, char const *argv[]){
  
  // Domain:
  //    [0,1] x [0,1]
  //
  // BC:
  //    u(x,0) = 0
  //    u(x,1) = 1
  //    u(0,y) = 0
  //    u(1,y) = 0

  int N = 400;                // number of columns
  int M = N;                  // number of rows
  float h = 1.0/N;            // grid width
  float tol = atof(argv[1]);  // tolerance
  int cont;                   // iterations counter

  int size = N * M;           // size of vector

  // vector to store grid values
  float* u_old = new float[size];
  float* u_new = new float[size];

  // OpenMP config
  int t = 1;
#ifdef _OPENMP
  // omp_set_num_threads(1);
  #pragma omp parallel
  t = omp_get_num_threads();
#endif
  cout << "Threads: " << t << endl;

  // seting grid values to 0 initially
  memset(u_old, 0.0, size * sizeof(float));
  memset(u_new, 0.0, size * sizeof(float));

  // boundary condition at u(x,1)
  for (int i = size - N; i < size; i++) {
    u_old[i] = 1.0;
    u_new[i] = 1.0;
  }

  cout << "Tolerance: " << setprecision(10) << tol << endl;

  // solving by the jacobi method (exercise 1)
  cout << "Solving by jacobi...\n";
  auto inicio = chrono::high_resolution_clock::now();

  cont = jacobi(tol, u_old, u_new, N, M, h, t);

  auto final = chrono::high_resolution_clock::now();
  chrono::duration<double> intervalo = final - inicio;

  cout << "\nTime spent (jacobi): " << intervalo.count() << "s\n";
  cout << "Iterations: " << cont << "\n";
  cout << "u(0.5, 0.5): " << u_old[ size/2 + N/2 ] << "\n\n";

  // saving data for plotting
  cout << "Saving data...\n";
  saveDataASCII("jacobi", u_old, N, M, h);
  saveDataBin(u_old, &N, &M, &h);
  cout << "Data saved succesfully!\n";

  // reseting values
  memset(u_old, 0.0, size * sizeof(float));
  memset(u_new, 0.0, size * sizeof(float));
  for (int i = size - N; i < size; i++) {
    u_old[i] = 1.0;
    u_new[i] = 1.0;
  }

  // solving by the SOR black-red method (exercise 2)
  float w = 1.;

  cout << "\nSolving by SOR with w = " + to_string(w) + "...\n";
  inicio = chrono::high_resolution_clock::now();

  cont = SOR(tol, u_old, u_new, N, M, h, t, w);

  final = chrono::high_resolution_clock::now();
  intervalo = final - inicio;

  cout << "\nTime spent (SOR): " << intervalo.count() << "s\n";
  cout << "Iterations: " << cont << "\n";
  cout << "u(0.5, 0.5): " << u_old[ size/2 + N/2 ] << "\n\n";

  // saving data for plotting
  cout << "Saving data...\n";
  saveDataASCII("SOR", u_old, N, M, h);
  saveDataBin(u_old, &N, &M, &h);
  cout << "Data saved succesfully!\n";

  return 0;
}
