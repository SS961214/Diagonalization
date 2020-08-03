#include "mersenne_twister.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <cublas_v2.h>
#include <curand_kernel.h>
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include "magma_v2.h"
#include "magma_lapack.h"

__global__ void SetRandomMatrix(int Dmat, magmaFloatComplex* mat, curandStateMtgp32_t* MTGPStates_d, int nBlock) {
  int i = blockIdx.x*blockDim.x +threadIdx.x;
  int j = blockIdx.y*blockDim.y +threadIdx.y;
  //int blockId = blockIdx.x +nBlock*blockIdx.y;
  if( i<Dmat && j<Dmat ){
    float rand1 = curand_normal(&MTGPStates_d[0]) /sqrt(2.0);
    float rand2 = curand_normal(&MTGPStates_d[0]) /sqrt(2.0);
    //printf("(%d,%d) rand1=%+.4lf, rand2=%+.4lf\n", i,j, rand1, rand2);

    if(i == j) mat[i +Dmat*j] = MAGMA_C_MAKE(sqrt(2.0)*rand1,0);
    else if(i < j) {
      mat[i +Dmat*j] = MAGMA_C_MAKE(rand1,rand2);
      mat[j +Dmat*i] = MAGMA_C_CONJ(mat[i +Dmat*j]);
    }
  }
}

int main(int argc, char **argv) {
    if(argc != 2) {
      fprintf(stderr, "Usage: 1.This 2.Dmat\n");
      exit(1);
    }

    const int Dmat = atoi(argv[1]);
    int nThread = (int)sqrt(1024);
    int nBlock = (int)Dmat/nThread;
    if( Dmat%nThread != 0 ) nBlock += 1;
    printf("nBlock=%d, nThread=%d *%d=%d\n", nBlock, nThread, nThread, nThread*nThread);
  
    dim3 dimGrid(nBlock, nBlock, 1);
    dim3 dimBlock(nThread, nThread, 1);

    const unsigned long long seed = 0;
    curandStateMtgp32 *MTGPStates_d;
    mtgp32_kernel_params *KernelParams_d;
    //curandGenerator_t generator;
    //curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_MTGP32);
    //curandSetPseudoRandomGeneratorSeed(generator, seed);

    void (*po)(int, magmaFloatComplex*, curandStateMtgp32_t*, int);
    po = SetRandomMatrix;
    struct cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr, po);
    printf("constSizeBytes     = %d\n", attr.constSizeBytes);
    printf("maxThreadsPerBlock = %d\n", attr.maxThreadsPerBlock);
    
      
    magma_init();
    magmaFloatComplex *mat_d, *EigenVectors_d, *temp_d;
    float *w, *rwork;
    magmaFloatComplex *wA, *work;
    magma_int_t lwork;
    magma_int_t lrwork = 1+5*Dmat+2*Dmat*Dmat;
    magma_int_t *iwork;
    magma_int_t liwork = 3+5*Dmat;
    magma_int_t info;
    magma_int_t nb = magma_get_zhetrd_nb(Dmat);

    struct timespec start_time, end_time;
    clock_t start, end;

    int temp1 = Dmat +Dmat*nb;
    int temp2 = 2*Dmat+Dmat*Dmat;
    if(temp1 < temp2) lwork = temp2;
    else lwork = temp1;

    magma_queue_t queue=NULL;
    magma_int_t dev = 0;
    magma_queue_create( dev, &queue );
    
    magma_cmalloc( &mat_d, Dmat*Dmat );
    magma_cmalloc( &EigenVectors_d, Dmat*Dmat );
    magma_cmalloc( &temp_d, Dmat*Dmat );

    magma_smalloc_cpu( &w, Dmat );
    magma_cmalloc_cpu( &wA, Dmat*Dmat );
    magma_cmalloc_cpu( &work, lwork );
    magma_smalloc_cpu( &rwork, lrwork );
    magma_imalloc_cpu( &iwork, liwork );
    
    cudaMalloc( (void**)&MTGPStates_d, nBlock*nBlock *sizeof(curandStateMtgp32) );
    cudaMalloc( (void**)&KernelParams_d, sizeof(mtgp32_kernel_params) );

    curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, KernelParams_d);
    curandMakeMTGP32KernelState(MTGPStates_d, mtgp32dc_params_fast_11213, KernelParams_d, 1, seed);
    SetRandomMatrix<<<dimGrid,dimBlock>>>(Dmat, mat_d, MTGPStates_d, nBlock);
    magma_cprint_gpu(Dmat, Dmat, mat_d, Dmat, queue); 
    magma_ccopymatrix(Dmat, Dmat, mat_d, Dmat, EigenVectors_d, Dmat, queue);

    start = clock();
    clock_gettime(CLOCK_REALTIME, &start_time);
    magma_cheevd_gpu(MagmaVec, MagmaUpper, Dmat, EigenVectors_d, Dmat, w, wA, Dmat, work, lwork, rwork, lrwork, iwork, liwork, &info);
    end = clock();
    clock_gettime(CLOCK_REALTIME, &end_time);

    unsigned int sec;
    int nsec;
    double d_sec;
    sec = end_time.tv_sec - start_time.tv_sec;
    nsec = end_time.tv_nsec - start_time.tv_nsec;
    d_sec = (double)sec +(double)nsec / (1000 * 1000 * 1000);
    fprintf(stderr, "%d %f %f\n", Dmat, (double)(end-start)/CLOCKS_PER_SEC, d_sec);

    /*
    magmaFloatComplex alpha = MAGMA_C_MAKE(1,0);
    magmaFloatComplex beta  = MAGMA_C_MAKE(0,0);
    magma_chemm(MagmaLeft, MagmaUpper, Dmat, Dmat, alpha, mat_d, Dmat, EigenVectors_d, Dmat, beta, temp_d, Dmat, queue);
    magma_cgemm(MagmaConjTrans, MagmaNoTrans, Dmat, Dmat, Dmat, alpha, EigenVectors_d, Dmat, temp_d, Dmat, beta, mat_d, Dmat, queue);
    magma_cprint_gpu(Dmat, Dmat, mat_d, Dmat, queue);
    magma_sprint(1, Dmat, w, 1);
    */
    
    magma_queue_destroy(queue);
    magma_free(mat_d);
    magma_free(EigenVectors_d);
    magma_free(temp_d);

    magma_free_cpu(w);
    magma_free_cpu(wA);
    magma_free_cpu(work);
    magma_free_cpu(rwork);
    magma_free_cpu(iwork);

    cudaFree(MTGPStates_d);
    cudaFree(KernelParams_d);

    magma_finalize();
    //curandDestroyGenerator(generator);
    
    return 0;
}
