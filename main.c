// Code is from our CSE260 homework: https://github.com/cse260-fa22/pa1-zhy025-zil060

#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
    struct timeval stop, start;
    long total_time_in_us;
    gettimeofday(&start, NULL);

#ifdef LOGGING
    // printf("Logging version\n");
#endif
    int n = (long)strtol(argv[1], NULL, 10); // Size of matrix
    // printf("n is %d\n", n);
    double *A, *B, *C;
    A = bl_malloc_aligned(n, n, sizeof(double));
    B = bl_malloc_aligned(n, n, sizeof(double));
    C = bl_malloc_aligned(n, n, sizeof(double));

    // Initialize A, B, C with random numbers
    fillMatRandom(A, n * n);
    fillMatRandom(B, n * n);
    fillMatRandom(C, n * n);

    // Compute C
    square_dgemm(n, A, B, C);

    // Lack of verification

    // Free heap
    free(A);
    free(B);
    free(C);

    gettimeofday(&stop, NULL);
    total_time_in_us = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
    printf("%lu\n", total_time_in_us);

    // for (int i = 0; i < 50; i++)
    // {
    //     printf("%p\t", ptrLog(i));
    // }
}
