#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "compress.h"

#define N 40

unsigned char * printBits( int a[][N] , int bits_dim1, int bits_dim2 )
{
  unsigned char * buf = (unsigned char *)malloc((bits_dim1 + 1) * bits_dim2 * sizeof(uint8_t));
  for (int i = 0 ; i < bits_dim1 ; i++) {
    for (int j = 0 ; j < bits_dim2 ; j++) {
      buf[i*(bits_dim1 + 1) + j] = (a[i][j] != 0) ? '#': '-';
    }
    buf[i*(bits_dim1 + 1) + bits_dim2] = '\n';
  }
  return buf;
}

void init_automat(int a[][N], int size)
{
  time_t t;
  srand((unsigned) time(&t));
  for (int i = 0 ; i < size ; i++) {
    for (int j = 0 ; j < size ; j++) {
      if (rand()%2 == 0) { a[i][j] = 1; }
    }
  }
}

const int rule = 2124;

void update_step_totalistic(int base[][N], int size, int rule, int ** A)
{

  int count;
  for (int i = 0 ; i < size ; i++) {
    memcpy(A[i], base[i], size * sizeof(int));

    for (int j = 0 ; j < size ; j++) {
      count = 0;
      for (int k = -1 ; k < 2 ; k ++) {
        for (int l = -1 ; l < 2 ; l++) {
          if (k != 0 | l != 0) {
            count += base[(i + k + size)%size][(j + l + size)%size];
          }
        }
      }
      if (rule & 1 << (2*count + base[i][j])) {
        A[i][j] = 1;
      } else {
        A[i][j] = 0;
      }
    }
  }
  for (int i = 0 ; i < size ; i++) {
    memcpy(base[i], A[i], size*sizeof(int));
  }
}

int main()
{
  int A[N][N]  = {{}};
  init_automat(A, N);
  int ** placeholder = (int **)malloc(N * sizeof(int *));
  for (int i = 0 ; i < N ; i ++) {
    placeholder[i] = (int *)malloc(N * sizeof(int));
  }
  printf("%i\n", compress_memory_size(printBits(A, N, N), (N+1)*N));
  for (int i = 0 ; i < 100 ; i ++) {
    printf("%s", printBits(A, N, N));
    update_step_totalistic(A, N, rule, placeholder);
    printf("%i\n", compress_memory_size(printBits(A, N, N), (N+1)*N));
  }
}
