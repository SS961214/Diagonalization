#include "mersenne_twister.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>

extern void zheevd_(char* JOBS, char* UPLO, int* N, double complex* A, int* LDA, double* W, double complex* WORK, int* LWORK, double* RWORK , int* LRWORK, int* IWORK, int* LIWORK, int* INFO);

extern void zhemm_(char* SIDE, char* UPLO, int* M, int* N, double complex* alpha, double complex* A, int* LDA, double complex* B, int* LDB, double complex* beta, double complex* C, int* LDC);

extern void zgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double complex* alpha, double complex* A, int* LDA, double complex* B, int* LDB, double complex* beta, double complex* C, int* LDC);

void print_matrix(int M, int N, double complex* mat) {
  printf("[\n");
  for(int i=0;i < M; ++i) {
    for(int j=0;j < N; ++j) printf("   %+.4f%+.4f*I", creal(mat[i+M*j]), cimag(mat[i+M*j]));
    printf("\n");
  }
  printf("];\n");
}

int main(int argc, char **argv) {
    if(argc != 2) {
      fprintf(stderr, "Usage: 1.This 2.Dmat\n");
      exit(1);
    }
    const int seed = 0;
    init_genrand(seed);

    int Dmat = atoi(argv[1]);
    double rand1, rand2, rand1_t, rand2_t;
    double complex *mat, *EigenVectors, *temp;
    
    double *w, *rwork;
    double complex *work;
    int lwork  = 2*Dmat+2*Dmat*Dmat;
    int liwork = 3+5*Dmat;
    int lrwork = 1 +5*Dmat +2*Dmat*Dmat;
    int *iwork;
    char JOBS='V';
    char UPLO='U';
    char SIDE='L';
    char nTrans='N';
    char Trans ='C';
    double complex alpha=1, beta=0;
    int info;

    struct timespec start_time, end_time;
    clock_t start, end;
    
    mat          = (double complex*)malloc( sizeof(double complex)*Dmat*Dmat );
    EigenVectors = (double complex*)malloc( sizeof(double complex)*Dmat*Dmat );
    temp         = (double complex*)malloc( sizeof(double complex)*Dmat*Dmat );
    w            = (double*)malloc( sizeof(double)*Dmat );
    work         = (double complex*)malloc( sizeof(double complex)*lwork );
    rwork        = (double*)malloc( sizeof(double)*lrwork );
    iwork        = (int*)malloc( sizeof(int)*liwork );
    
    for(int i=0;i < Dmat; ++i) {
      rand1 = genrand_real3();
      rand2 = genrand_real3();
      mat[i+Dmat*i] = sqrt(-2*log(rand1)) *cos(2*M_PI*rand2);
      for(int j=0;j < i; ++j) {
	rand1_t = genrand_real3();
	rand2_t = genrand_real3();
	rand1 = sqrt(-2*log(rand1_t)) *cos(2*M_PI*rand2_t);
	rand2 = sqrt(-2*log(rand1_t)) *sin(2*M_PI*rand2_t);
	mat[i+Dmat*j] = (rand1 +I*rand2) / sqrt(2.0);
	mat[j+Dmat*i] = conj(mat[i+Dmat*j]);
      }
    }
    for(int i=0;i < Dmat; ++i) for(int j=0;j < Dmat; ++j) EigenVectors[i+Dmat*j] = mat[i+Dmat*j];
    //    print_matrix(Dmat, Dmat, mat);

    start = clock();
    clock_gettime(CLOCK_REALTIME, &start_time);
    zheevd_(&JOBS, &UPLO, &Dmat, EigenVectors, &Dmat, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    end = clock();
    clock_gettime(CLOCK_REALTIME, &end_time);
    unsigned int sec;
    int nsec;
    double d_sec;
    sec = end_time.tv_sec - start_time.tv_sec;
    nsec = end_time.tv_nsec - start_time.tv_nsec;
    d_sec = (double)sec +(double)nsec / (1000 * 1000 * 1000);
    fprintf(stderr, "%d %f %f\n", Dmat, (double)(end-start)/CLOCKS_PER_SEC, d_sec);

    print_matrix(Dmat, Dmat, mat);
    
    zhemm_(&SIDE, &UPLO, &Dmat, &Dmat, &alpha, mat, &Dmat, EigenVectors, &Dmat, &beta, temp, &Dmat);
    zgemm_(&Trans, &nTrans, &Dmat, &Dmat, &Dmat, &alpha, EigenVectors, &Dmat, temp, &Dmat, &beta, mat, &Dmat);

    print_matrix(Dmat, Dmat, mat);
    
    return 0;
}
