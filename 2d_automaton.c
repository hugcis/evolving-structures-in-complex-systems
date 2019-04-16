#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define N 100

unsigned char * printBits( int a[][N] , int bits_dim1, int bits_dim2 )
{
  unsigned char * buf = (unsigned char *)malloc((bits_dim1 + 1) * bits_dim2 * sizeof(uint8_t));
  for (int i = 0 ; i < bits_dim1 ; i++) {
    for (int j = 0 ; j < bits_dim2 ; j++) {
      buf[i*(bits_dim1 + 1) + j] = (a[i][j] != 0) ? '1': '0';
    }
    buf[i*(bits_dim1 + 1) + bits_dim2] = '\n';
  }
  return buf;
}


int main()
{
  int A[N][N]  = {{}};
  A[5][5] = 1;
  printf("%s", printBits(A, N, N));
}
