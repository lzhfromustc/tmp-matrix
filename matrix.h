
// parameters for preparation
#define RAND_MAX 0x7fffffff
#define GEMM_SIMD_ALIGN_SIZE 32

// Tiling
#define DGEMM_KC 64
#define DGEMM_MC 64
#define DGEMM_NC 64
#define DGEMM_MR 4
#define DGEMM_NR 4

double *bl_malloc_aligned(
    int m,
    int n,
    int size);

void fillMatRandom(double *p, int n);

void square_dgemm(int lda, double *A, double *B, double *C);

int intLogIndex();

void *ptrLog();
