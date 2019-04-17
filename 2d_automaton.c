#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "compress.h"

#define N 100
#define STEPS 100

char * printBits( int a[][N] , int bits_dim1, int bits_dim2 )
{
  char * buf = (char *)malloc(((bits_dim1 + 1) * bits_dim2 + 1) * sizeof(uint8_t));
  for (int i = 0 ; i < bits_dim1 ; i++) {
    for (int j = 0 ; j < bits_dim2 ; j++) {
      buf[i*(bits_dim1 + 1) + j] = (a[i][j] != 0) ? '#': '-';
    }
    buf[i*(bits_dim1 + 1) + bits_dim2] = '\n';
  }
  buf[(bits_dim1 + 1) * bits_dim2] = '\0';
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

void write_to_file(int rule) {
  FILE *out_file;
  char * fname;
  asprintf(&fname, "data_2d/out%i.dat", rule);
  out_file = fopen(fname, "w+");

  int A[N][N]  = {{}};
  init_automat(A, N);
  int ** placeholder = (int **)malloc(N * sizeof(int *));
  for (int i = 0 ; i < N ; i ++) {
    placeholder[i] = (int *)malloc(N * sizeof(int));
  }

  for (int i = 0 ; i < STEPS ; i ++) {
    update_step_totalistic(A, N, rule, placeholder);
    if (i%30 == 0) {
      char * out_string = printBits(A, N, N);
      fprintf(out_file, "%i    %i\n", i,
              compress_memory_size(out_string, (N+1)*N));
      free(out_string);
    }
  }

  for (int i = 0 ; i < N ; i ++) {
    free(placeholder[i]);
  }
  free(placeholder);
  fclose(out_file);
}

int main()
{
  printf("%i\n", RAND_MAX);
  for (int rule = 0 ; rule < 1 ; rule++) {
    printf("\rRule %i", rule);
    write_to_file(rule);
  }
}
