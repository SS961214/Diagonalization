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

__global__ void setRandomMatrix(int Dmat, magmaDoubleComplex* mat_d, curandStateMtgp32_t* MTGPStates_d, int nBlock) {
  int index1 = blockIdx.x*blockDim.x +threadIdx.x;
  int index2 = blockIdx.y*blockDim.y +threadIdx.y;
  //int blockId = blockIdx.x +nBlock*blockIdx.y;
  if( index1<Dmat && index2<Dmat ){
    double rand1 = curand_normal_double(&MTGPStates_d[0]) /sqrt(2.0);
    double rand2 = curand_normal_double(&MTGPStates_d[0]) /sqrt(2.0);
    //printf("%+.4lf, %+.4lf\n", rand1, rand2);
    if(index1 == index2) mat_d[index1 +Dmat*index2] = MAGMA_Z_MAKE(sqrt(2.0)*rand1,0);
    else if(index1 < index2) {
      mat_d[index1 +Dmat*index1] = MAGMA_Z_MAKE(rand1,rand2);
      mat_d[index2 +Dmat*index1] = MAGMA_Z_CONJ(mat[index1 +Dmat*index2]);
    }
  }
}

int main(int argc, char **argv) {
    if(argc != 2) {
      fprintf(stderr, "Usage: 1.This 2.Dmat\n");
      exit(EX_USAGE);
    }
    void (*funcPtr)(int, magmaDoubleComplex*, curandStateMtgp32_t*, int);
    funcPtr = setRandomMatrix;
    struct cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr, funcPtr);
    fprintf(stdout, "# constSizeBytes     = %zu\n", attr.constSizeBytes);
    fprintf(stdout, "# maxThreadsPerBlock = %d\n", attr.maxThreadsPerBlock);

    const int Dmat = atoi(argv[1]);
    int nThread = attr.maxThreadsPerBlock;
    int nBlock = (int)Dmat/nThread;
    if( Dmat%nThread != 0 ) nBlock += 1;
    fprintf(stdout, "# (Dmat=%d) nBlock=%d, nThread=%d *%d=%d\n", Dmat, nBlock, nThread, nThread, nThread*nThread);

    dim3 dimGrid(nBlock, nBlock, 1);
    dim3 dimBlock(nThread, nThread, 1);

    const unsigned long long seed = 0;
    curandStateMtgp32    *MTGPStates_d = NULL;
    mtgp32_kernel_params *KernelParams_d = NULL;

    magma_init();
    magmaDoubleComplex *mat_d = NULL, *EigenVectors_d = NULL, *temp_d = NULL;
    magmaDoubleComplex alpha = MAGMA_Z_MAKE(1,0);
    magmaDoubleComplex beta  = MAGMA_Z_MAKE(0,0);
    double *EigenValues = NULL, *rwork = NULL, lrwork_t;
    magmaDoubleComplex *wA = NULL, *work = NULL, lwork_t;
    magma_int_t *iwork = NULL, liwork_t;
    magma_int_t lwork  = -1;
    magma_int_t lrwork = -1;
    magma_int_t liwork = -1;
    magma_int_t info;
    magma_zheevd_gpu(MagmaVec, MagmaUpper, Dmat, EigenVectors_d, Dmat, EigenValues, wA, Dmat, &lwork_t, lwork, &lrwork_t, lrwork, &liwork_t, liwork, &info);
    lwork  = (magma_int_t)MAGMA_Z_REAL(lwork_t);
    lrwork = (magma_int_t)lrwork_t;
    liwork = liwork_t;

    magma_queue_t queue = NULL;
    magma_int_t   dev   = 0;

    double start, end, T_diag, T_prod;
    // ******************** Allocation ******************** //
    magma_queue_create( dev, &queue );
    magma_zmalloc( &mat_d, Dmat*Dmat );
    magma_zmalloc( &EigenVectors_d, Dmat*Dmat );
    magma_zmalloc( &temp_d, Dmat*Dmat );

    magma_dmalloc_cpu( &EigenValues, Dmat );
    magma_zmalloc_cpu( &wA, Dmat*Dmat );
    magma_zmalloc_cpu( &work, lwork );
    magma_dmalloc_cpu( &rwork, lrwork );
    magma_imalloc_cpu( &iwork, liwork );

    cudaMalloc( (void**)&MTGPStates_d, nBlock*nBlock *sizeof(curandStateMtgp32) );
    cudaMalloc( (void**)&KernelParams_d, sizeof(mtgp32_kernel_params) );
    // ******************** (END)Allocation ******************** //

    curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, KernelParams_d);
    curandMakeMTGP32KernelState(MTGPStates_d, mtgp32dc_params_fast_11213, KernelParams_d, 1, seed);
    setRandomMatrix<<<dimGrid,dimBlock>>>(Dmat, mat_d, MTGPStates_d, nBlock);
    //magma_zprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    magma_zcopymatrix(Dmat, Dmat, mat_d, Dmat, EigenVectors_d, Dmat, queue);

    start = getETtime();
      magma_zheevd_gpu(MagmaVec, MagmaUpper, Dmat, EigenVectors_d, Dmat, EigenValues, wA, Dmat, work, lwork, rwork, lrwork, iwork, liwork, &info);
    end = getETtime();
      if(info != 0) {
        fprintf(stderr, "# Error: magma_zheevd_gpu failed.\n");
        exit(EX_SOFTWARE);
      }
    T_diag = end-start;

    start = getETtime();
      magma_zhemm(MagmaLeft, MagmaUpper, Dmat, Dmat, alpha, mat_d, Dmat, EigenVectors_d, Dmat, beta, temp_d, Dmat, queue);
      magma_zgemm(MagmaConjTrans, MagmaNoTrans, Dmat, Dmat, Dmat, alpha, EigenVectors_d, Dmat, temp_d, Dmat, beta, mat_d, Dmat, queue);
    end = getETtime();
    T_prod = end-start;
    // magma_zprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    // magma_dprint(1, Dmat, w, 1);

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
