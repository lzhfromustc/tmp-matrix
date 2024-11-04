#include "matrix.h"
#include <stdlib.h>
#include "stdio.h"

double *bl_malloc_aligned(
    int m,
    int n,
    int size)
{
    double *ptr;
    int err;

    err = posix_memalign((void **)&ptr, (size_t)GEMM_SIMD_ALIGN_SIZE, size * m * n);

    if (err)
    {
        printf("bl_malloc_aligned(): posix_memalign() failures");
        exit(1);
    }

    return ptr;
}

void fillMatRandom(double *p, int n)
{
    long int Rmax = RAND_MAX;
    long int Rmax_2 = Rmax >> 1;
    long int RM = Rmax_2 + 1;
    for (int i = 0; i < n; ++i)
    {
        long int r = rand(); // Uniformly distributed ints over [0,RAND_MAX]
                             // Typical value of RAND_MAX: 2^31 - 1
        long int R = r - RM;
        p[i] = (double)R / (double)RM; // Uniformly distributed over [-1, 1]
    }
}