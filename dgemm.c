#include "matrix.h"

#define min(i, j) i < j ? i : j
#define a(i, j, ld) a[(i) * (ld) + (j)]
#define b(i, j, ld) b[(i) * (ld) + (j)]
#define c(i, j, ld) c[(i) * (ld) + (j)]

void *logging[100000000];
int logIndex = 0;

int intLogIndex()
{
    return logIndex;
}

void *ptrLog(int i)
{
    if (i < logIndex)
    {
        return logging[i];
    }
    else
    {
        return logging[logIndex];
    }
}

void bl_dgemm_ukr(int k,
                  int m,
                  int n,
                  double *a,
                  double *b,
                  double *c,
                  unsigned long long ldc)
{
    int l, j, i;

    // #ifdef LOGGING
    //     logging[logIndex] = (void *)a;
    //     logIndex++;
    //     logging[logIndex] = (void *)b;
    //     logIndex++;
    //     logging[logIndex] = (void *)c;
    //     logIndex++;
    // #endif

    for (l = 0; l < k; ++l)
    {
        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < m; ++i)
            {
                c(i, j, ldc) += a(i, l, ldc) * b(l, j, ldc);
#ifdef LOGGING
                logging[logIndex] = &a(i, l, ldc);
                logIndex++;
                logging[logIndex] = &b(l, j, ldc);
                logIndex++;
                logging[logIndex] = &c(i, j, ldc);
                logIndex++;
                logging[logIndex] = &c(i, j, ldc);
                logIndex++;
#endif
            }
        }
    }
}

static inline void bl_macro_kernel(
    int m,
    int n,
    int k,
    double *packA,
    double *packB,
    double *C,
    int ldc)
{
    int i, j;

    for (i = 0; i < m; i += DGEMM_MR)
    { // 2-th loop around micro-kernel
        for (j = 0; j < n; j += DGEMM_NR)
        { // 1-th loop around micro-kernel
            bl_dgemm_ukr(
                k,
                min(m - i, DGEMM_MR),
                min(n - j, DGEMM_NR),
                &packA[i * ldc],
                &packB[j],
                &C[i * ldc + j],
                (unsigned long long)ldc);
        } // 1-th loop around micro-kernel
    }     // 2-th loop around micro-kernel
}

void bl_dgemm(
    int m,
    int n,
    int k,
    double *XA,
    int lda,
    double *XB,
    int ldb,
    double *C,
    int ldc)
{
    int ic, ib, jc, jb, pc, pb;
    double *packA, *packB;

    for (ic = 0; ic < m; ic += DGEMM_MC)
    { // 5-th loop around micro-kernel
        ib = min(m - ic, DGEMM_MC);
        for (pc = 0; pc < k; pc += DGEMM_KC)
        { // 4-th loop around micro-kernel
            pb = min(k - pc, DGEMM_KC);

            // Currently no pack
#ifdef LOGGING
            logging[logIndex] = (void *)XA;
            logIndex++;
#endif
            packA = &XA[pc + ic * lda];

            for (jc = 0; jc < n; jc += DGEMM_NC)
            { // 3-rd loop around micro-kernel
                jb = min(m - jc, DGEMM_NC);

                // Currently no pack
#ifdef LOGGING
                logging[logIndex] = (void *)XB;
                logIndex++;
#endif
                packB = &XB[ldb * pc + jc];

                bl_macro_kernel(
                    ib,
                    jb,
                    pb,
                    packA,
                    packB,
                    &C[ic * ldc + jc],
                    ldc);
            } // End 3.rd loop around micro-kernel
        }     // End 4.th loop around micro-kernel
    }         // End 5.th loop around micro-kernel
}

void square_dgemm(int lda, double *A, double *B, double *C)
{
    bl_dgemm(lda, lda, lda, A, lda, B, lda, C, lda);
}
