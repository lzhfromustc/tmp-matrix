# rm ./perf.data*
gcc -O4 -g -o origin main.c matrix.h dgemm.c utils.c
