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

__global__ void setRandomMatrix(int Dmat, magmaFloatComplex* mat, curandStateMtgp32_t* MTGPStates_d, int nBlock) {
  int index1 = blockIdx.x*blockDim.x +threadIdx.x;
  int index2 = blockIdx.y*blockDim.y +threadIdx.y;
  //int blockId = blockIdx.x +nBlock*blockIdx.y;
  if( index1<Dmat && index2<Dmat ){
    float rand1 = curand_normal(&MTGPStates_d[0]) /sqrt(2.0);
    float rand2 = curand_normal(&MTGPStates_d[0]) /sqrt(2.0);
    if(index1 == index2) mat[index1 +Dmat*index2] = MAGMA_C_MAKE(sqrt(2.0)*rand1,0);
    else if(i < j) {
      mat[index1 +Dmat*index2] = MAGMA_C_MAKE(rand1,rand2);
      mat[index2 +Dmat*index1] = MAGMA_C_CONJ(mat[index1 +Dmat*index2]);
    }
  }
}

int main(int argc, char **argv) {
    if(argc != 2) {
      fprintf(stderr, "Usage: 1.This 2.Dmat\n");
      exit(EX_USAGE);
    }
    void (*funcPtr)(int, magmaFloatComplex*, curandStateMtgp32_t*, int);
    funcPtr = SetRandomMatrix;
    struct cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr, funcPtr);
    printf("constSizeBytes     = %d\n", attr.constSizeBytes);
    printf("maxThreadsPerBlock = %d\n", attr.maxThreadsPerBlock);

    const int Dmat = atoi(argv[1]);
    int nThread = attr.maxThreadsPerBlock;
    int nBlock = (int)Dmat/nThread;
    if( Dmat%nThread != 0 ) nBlock += 1;
    printf("nBlock=%d, nThread=%d *%d=%d\n", nBlock, nThread, nThread, nThread*nThread);

    dim3 dimGrid(nBlock, nBlock, 1);
    dim3 dimBlock(nThread, nThread, 1);

    const unsigned long long seed = 0;
    curandStateMtgp32    *MTGPStates_d;
    mtgp32_kernel_params *KernelParams_d;

    magma_init();
    magmaFloatComplex *mat_d, *EigenVectors_d, *temp_d;
    magmaFloatComplex alpha = MAGMA_C_MAKE(1,0);
    magmaFloatComplex beta  = MAGMA_C_MAKE(0,0);
    float *EigenValues, *rwork, lrwork_t;
    magmaFloatComplex *wA, *work, lwork_t;
    magma_int_t *iwork, liwork_t;
    magma_int_t lwork  = -1;
    magma_int_t lrwork = -1;
    magma_int_t liwork = -1;
    magma_int_t info;
    magma_cheevd_gpu(MagmaVec, MagmaUpper, Dmat, EigenVectors_d, Dmat, w, wA, Dmat, &lwork_t, lwork, &lrwork_t, lrwork, &liwork_t, liwork, &info);
    lwork  = lwork_t;
    lrwork = lrwork_t;
    liwork = liwork_t;

    magma_queue_t queue =  dNULL;
    magma_int_t   dev   = 0;

    double start, end, T_diag, T_prod;
    // ******************** Allocation ******************** //
    magma_queue_create( dev, &queue );
    magma_cmalloc( &mat_d, Dmat*Dmat );
    magma_cmalloc( &EigenVectors_d, Dmat*Dmat );
    magma_cmalloc( &temp_d, Dmat*Dmat );

    magma_smalloc_cpu( &EigenValues, Dmat );
    magma_cmalloc_cpu( &wA, Dmat*Dmat );
    magma_cmalloc_cpu( &work, lwork );
    magma_smalloc_cpu( &rwork, lrwork );
    magma_imalloc_cpu( &iwork, liwork );

    cudaMalloc( (void**)&MTGPStates_d, nBlock*nBlock *sizeof(curandStateMtgp32) );
    cudaMalloc( (void**)&KernelParams_d, sizeof(mtgp32_kernel_params) );
    // ******************** (END)Allocation ******************** //

    curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, KernelParams_d);
    curandMakeMTGP32KernelState(MTGPStates_d, mtgp32dc_params_fast_11213, KernelParams_d, 1, seed);
    SetRandomMatrix<<<dimGrid,dimBlock>>>(Dmat, mat_d, MTGPStates_d, nBlock);
    magma_cprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    magma_ccopymatrix(Dmat, Dmat, mat_d, Dmat, EigenVectors_d, Dmat, queue);

    start = getETtime();
      magma_cheevd_gpu(MagmaVec, MagmaUpper, Dmat, EigenVectors_d, Dmat, EigenValues, wA, Dmat, work, lwork, rwork, lrwork, iwork, liwork, &info);
    end = getETtime();
    T_diag = end-start;

    start = getETtime();
      magma_chemm(MagmaLeft, MagmaUpper, Dmat, Dmat, alpha, mat_d, Dmat, EigenVectors_d, Dmat, beta, temp_d, Dmat, queue);
      magma_cgemm(MagmaConjTrans, MagmaNoTrans, Dmat, Dmat, Dmat, alpha, EigenVectors_d, Dmat, temp_d, Dmat, beta, mat_d, Dmat, queue);
    end = getETtime();
    T_prod = end-start;
    magma_cprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    magma_sprint(1, Dmat, EigenValues, 1);

    fprintf(stderr, "Real_time(sec): %d %lf %lf\n", Dmat, T_diag, T_prod);
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
