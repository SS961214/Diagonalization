#include "magma_v2.h"
#include "magma_lapack.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sysexits.h>
#include <cublas_v2.h>
#include <curand_kernel.h>
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>


static inline double getETtime() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec +(double)tv.tv_usec *1e-6;
}

__global__ void setRandomMatrix(int Dmat, int Dmat32, magmaFloatComplex* mat, curandStateMtgp32_t* MTGPStates_d, int nBlock) {
  int index1 = blockIdx.x*blockDim.x +threadIdx.x;
  int index2 = blockIdx.y*blockDim.y +threadIdx.y;
  if( (index1>=Dmat32) || (index2>=Dmat32) ) return;
  if( (index1>=Dmat) || (index2>=Dmat) ) mat[index1 +Dmat32*index2] = MAGMA_C_ZERO;
  else if( index2<=index1 && index1<Dmat ){
    float rand1 = curand_normal(&MTGPStates_d[0]) /sqrt(2.0);
    float rand2 = curand_normal(&MTGPStates_d[0]) /sqrt(2.0);
    if(index1 == index2) mat[index1 +Dmat32*index2] = MAGMA_C_MAKE(sqrt(2.0)*rand1,0);
    else if(index1 < index2) {
      mat[index1 +Dmat32*index2] = MAGMA_C_MAKE(rand1,rand2);
      mat[index2 +Dmat32*index1] = MAGMA_C_CONJ(mat[index1 +Dmat32*index2]);
    }
  }
}

int main(int argc, char **argv) {
    if(argc != 2) {
      fprintf(stderr, "Usage: 1.This 2.Dmat\n");
      exit(EX_USAGE);
    }
    void (*funcPtr)(int, int, magmaFloatComplex*, curandStateMtgp32_t*, int);
    funcPtr = setRandomMatrix;
    struct cudaFuncAttributes attr;
    cudaError_t err;
    err = cudaFuncGetAttributes(&attr, funcPtr);
    if (err != cudaSuccess) {
      fprintf(stderr, "Error: cudaFuncGetAttributes failed with an exit code = %d.\n       %s", err, cudaGetErrorString(err));
      exit(err);
    }
    fprintf(stdout, "# constSizeBytes     = %zu\n", attr.constSizeBytes);
    fprintf(stdout, "# maxThreadsPerBlock = %d\n", attr.maxThreadsPerBlock);

    const int Dmat = atoi(argv[1]);
    int Dmat32;
    if(Dmat%32 == 0) Dmat32 = Dmat;
    else Dmat32 = (Dmat/32+1)*32;
    int nThread = (int)sqrt(attr.maxThreadsPerBlock);
    int nBlock = (int)Dmat32/nThread;
    if( Dmat32%nThread != 0 ) nBlock += 1;
    fprintf(stdout, "# (Dmat=%d, Dmat32=%d) nBlock=%d, nThread=%d *%d=%d\n", Dmat, Dmat32, nBlock, nThread, nThread, nThread*nThread);

    dim3 dimGrid(nBlock, nBlock, 1);
    dim3 dimBlock(nThread, nThread, 1);

    const unsigned long long seed = 0;
    curandStateMtgp32    *MTGPStates_d = NULL;
    mtgp32_kernel_params *KernelParams_d = NULL;

    magma_init();
    magmaFloatComplex *mat_d = NULL, *EigenVectors_d = NULL, *temp_d = NULL;
    magmaFloatComplex alpha = MAGMA_C_MAKE(1,0);
    magmaFloatComplex beta  = MAGMA_C_MAKE(0,0);
    float *EigenValues = NULL, *rwork = NULL, lrwork_t;
    magmaFloatComplex *wA = NULL, *work = NULL, lwork_t;
    magma_int_t *iwork = NULL, liwork_t;
    magma_int_t lwork  = -1;
    magma_int_t lrwork = -1;
    magma_int_t liwork = -1;
    magma_int_t info;
    magma_cheevd_gpu(MagmaVec, MagmaUpper, Dmat32, EigenVectors_d, Dmat32, EigenValues, wA, Dmat32, &lwork_t, lwork, &lrwork_t, lrwork, &liwork_t, liwork, &info);
    lwork  = (magma_int_t)MAGMA_C_REAL(lwork_t);
    lrwork = (magma_int_t)lrwork_t;
    liwork = liwork_t;

    magma_queue_t queue = NULL;
    magma_int_t   dev   = 0;

    double start, end, T_diag, T_prod;
    // ******************** Allocation ******************** //
    magma_queue_create( dev, &queue );
    magma_cmalloc( &mat_d, Dmat32*Dmat32 );
    magma_cmalloc( &EigenVectors_d, Dmat32*Dmat32 );
    magma_cmalloc( &temp_d, Dmat32*Dmat32 );

    magma_smalloc_cpu( &EigenValues, Dmat32 );
    magma_cmalloc_cpu( &wA, Dmat32*Dmat32 );
    magma_cmalloc_cpu( &work, lwork );
    magma_smalloc_cpu( &rwork, lrwork );
    magma_imalloc_cpu( &iwork, liwork );

    cudaMalloc( (void**)&MTGPStates_d, nBlock*nBlock *sizeof(curandStateMtgp32) );
    cudaMalloc( (void**)&KernelParams_d, sizeof(mtgp32_kernel_params) );
    // ******************** (END)Allocation ******************** //

    curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, KernelParams_d);
    curandMakeMTGP32KernelState(MTGPStates_d, mtgp32dc_params_fast_11213, KernelParams_d, 1, seed);
    setRandomMatrix<<<dimGrid,dimBlock>>>(Dmat, Dmat32, mat_d, MTGPStates_d, nBlock);
    // magma_cprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    magma_ccopymatrix(Dmat32, Dmat32, mat_d, Dmat32, EigenVectors_d, Dmat32, queue);

    start = getETtime();
      magma_cheevd_gpu(MagmaVec, MagmaUpper, Dmat32, EigenVectors_d, Dmat32, EigenValues, wA, Dmat32, work, lwork, rwork, lrwork, iwork, liwork, &info);
    end = getETtime();
      if(info != 0) {
        fprintf(stderr, "# Error: magma_cheevd_gpu failed. (info=%d)\n", info);
        magma_xerbla("magma_cheevd_gpu", info);
        magma_strerror(info);
        exit(EX_SOFTWARE);
      }
    T_diag = end-start;

    start = getETtime();
      magma_chemm(MagmaLeft, MagmaUpper, Dmat32, Dmat32, alpha, mat_d, Dmat32, EigenVectors_d, Dmat32, beta, temp_d, Dmat32, queue);
      magma_cgemm(MagmaConjTrans, MagmaNoTrans, Dmat32, Dmat32, Dmat32, alpha, EigenVectors_d, Dmat32, temp_d, Dmat32, beta, mat_d, Dmat32, queue);
    end = getETtime();
    T_prod = end-start;
    // magma_cprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    // magma_sprint(1, Dmat, EigenValues, 1);

    fprintf(stderr, "Real_time(sec): %d %d %lf %lf\n", Dmat, Dmat32, T_diag, T_prod);
    // ******************** Free ******************** //
    magma_queue_destroy(queue);
    magma_free(mat_d);
    magma_free(EigenVectors_d);
    magma_free(temp_d);

    magma_free_cpu(EigenValues);
    magma_free_cpu(wA);
    magma_free_cpu(work);
    magma_free_cpu(rwork);
    magma_free_cpu(iwork);

    cudaFree(MTGPStates_d);
    cudaFree(KernelParams_d);
    // ******************** (END)Free ******************** //
    magma_finalize();
    return 0;
}
