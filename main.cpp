#include <stdio.h>
#include <iostream>
#include <cstring>
#include <fstream>
#include "omp.h"

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
 * @return false When the tolerance of the error isn't enough
 */
bool jacobi(float* u_old, float* u_new, int N, int M, float h) {

  float tol;    // tolerance
  float sum;    // sum of the differences between each term of the vectors

  float h2 = h*h;

  // calculate next approximation for the solution, which is basically the mean of it's neighbours
  // note that it doesn't iterate through the boundary
  for (int j = 1; j < M - 1; j++) {
    for (int i = 1; i < N - 1; i++) {
      u_new[ j*N + i ] = ( 
        h2 * S (i*h, j*h) 
        + u_old[ j*N + (i - 1) ] 
        + u_old[ j*N + (i + 1) ] 
        + u_old[ (j - 1)*N + i ] 
        + u_old[ (j + 1)*N + i ] 
      )/4;
    }
  }

  tol = (1/((N - 2) * (M - 2))) * 1;  
  

  if (tol < .000001) {
    return true;
  }
  return false;

}

void saveResult(float* u, int N, int M, float h) {

  FILE* file;

  file = fopen ("output.bin", "wb");

  fwrite(&h, sizeof(float), 1, file);
  fwrite(&N, sizeof(int), 1, file);
  fwrite(&M, sizeof(int), 1, file);
  fwrite(u, sizeof(float), sizeof(u), file);

  fclose(file);

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

  int N = 10;           // number of columns
  int M = N;            // number of rows
  float h = 1/N;        // grid width

  int size = N * M;     // size of vector

  bool solved = false;

  float* u_old = new float[size];
  float* u_new = new float[size];

  // setting grid values to 0 initially
  memset(u_old, 0.0, size * sizeof(float));

  // boundary condition at u(x,1)
  for (int i = size - N; i < size; i++) {
    u_old[i] = 1.0;
  }

  // solving by the jacobi method
  while (!solved) {
    solved = jacobi(u_old, u_new, N, M, h);
  }

  saveResult(u_new, N, M, h);

  return 0;
}