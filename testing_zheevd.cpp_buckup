//LU factorization MAGMA public domain
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <cuda.h>
#include <cublas_v2.h>
#include <curand_kernel.h>
#include <cuda_runtime_api.h>
#include "magma_v2.h"
#include "magma_lapack.h"
//#include "testings.h"
//Matlab/Octave format
void printmat(int N, int M, double *A, int LDA) {
double mtmp;
 printf("[ ");
 for (int i = 0; i < N; i++) {
   printf("[ ");
   for (int j = 0; j < M; j++) {
     mtmp = A[i + j * LDA];
     printf("%5.2e", mtmp);
     if (j < M - 1) printf(", ");
   } if (i < N - 1) printf("]; ");
   else printf("] ");
 } printf("]");
}

int main() {
  magma_init();
  
  int M=3, N=3, lda;
  magma_int_t min_mn = 3;
  magma_int_t *ipiv, *ipiv_d, info, info_d;
  double *A, *A_d;
  magma_queue_t queue=NULL;
  magma_int_t dev = 0;
  magma_queue_create( dev, &queue );
    
  min_mn = (magma_int_t)fmin(M,N);
  lda = N;
  
  magma_dmalloc(&A_d, M*N);
  ipiv = (magma_int_t*)malloc( min_mn*sizeof(magma_int_t));
  A = (double*)malloc( M*N*sizeof(double) );

  A[0+0*lda]=1; A[0+1*lda]= 8; A[0+2*lda]= 3;
  A[1+0*lda]=2; A[1+1*lda]=10; A[1+2*lda]= 8;
  A[2+0*lda]=9; A[2+1*lda]=-5; A[2+2*lda]=-1;
  printf("A =");printmat(M,N,A,lda);printf("\n");
  printf("min_mn=%d\n", min_mn);

  magma_dsetmatrix(M, N, A, lda, A_d, lda, queue);
  magma_dgetrf_gpu(M, N, A_d, lda, ipiv, &info);
  printf("Matric Ad on GPU:");magma_dprint_gpu(M, N, A_d, lda, queue);

  // magma_dgetrf(M, N, A, lda, ipiv, &info);
  printf("Matrix A on CPU:");printmat(M,N,A,lda);printf("\n");

  magma_dgetmatrix(M, N, A_d, lda, A, lda, queue);
  printf("Matrix A on CPU:");printmat(M,N,A,lda);printf("\n");  

  magma_queue_destroy( queue );
  magma_free(A_d);
  free(ipiv);
  free(A);
  magma_finalize();
  return 0;
}
