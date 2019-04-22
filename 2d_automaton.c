#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "compress.h"

#define N 100
#define STEPS 1000

char* print_bits(int a[][N], size_t bits_dim1, size_t bits_dim2)
{
  char* buf = (char*)malloc(((bits_dim1 + 1) * bits_dim2 + 1)
      * sizeof(uint8_t));
  for (size_t i = 0; i < bits_dim1; i++) {
    for (size_t j = 0; j < bits_dim2; j++) {
      buf[i * (bits_dim1 + 1) + j] = (a[i][j] != 0) ? '1' : '0';
    }
    buf[i * (bits_dim1 + 1) + bits_dim2] = '\n';
  }
  buf[(bits_dim1 + 1) * bits_dim2] = '\0';
  return buf;
}

void init_automat(int a[][N], size_t size)
{
  assert(size > 20);
  time_t t;
  srand((unsigned)time(&t));
  for (size_t i = size / 2 - 10; i < size / 2 + 10; i++) {
    for (size_t j = size / 2 - 10; j < size / 2 + 10; j++) {
      if (rand() % 2 == 0) {
        a[i % size][j % size] = 1;
      }
    }
  }
}

void update_step_totalistic(int base[][N], size_t size, int rule, int** A)
{
  int count;
  for (size_t i = 0; i < size; i++) {
    memcpy(A[i], base[i], size * sizeof(int));

    for (size_t j = 0; j < size; j++) {
      count = 0;
      for (int k = -1; k < 2; k++) {
        for (int l = -1; l < 2; l++) {
          if (k != 0 | l != 0) {
            count += base[(i + k + size) % size][(j + l + size) % size];
          }
        }
      }
      if (rule & 1 << (2 * count + base[i][j])) {
        A[i][j] = 1;
      } else {
        A[i][j] = 0;
      }
    }
  }
  for (size_t i = 0; i < size; i++) {
    memcpy(base[i], A[i], size * sizeof(int));
  }
}

void write_to_file(int rule, int save_steps)
{
  FILE* out_file;
  char* fname;
  asprintf(&fname, "data_2d/out%i.dat", rule);
  out_file = fopen(fname, "w+");

  int A[N][N] = { {} };
  init_automat(A, N);
  int** placeholder = (int**)malloc(N * sizeof(int*));
  for (int i = 0; i < N; i++) {
    placeholder[i] = (int*)malloc(N * sizeof(int));
  }

  for (int i = 0; i < STEPS; i++) {
    update_step_totalistic(A, N, rule, placeholder);
    if (i % 50 == 0) {
      char* out_string = print_bits(A, N, N);
      int compressed_size = compress_memory_size(out_string, (N + 1) * N);

      fprintf(out_file, "%i    %i\n", i, compressed_size);

      if (save_steps == 1) {
        FILE* out_step_file;
        char* step_fname;
        asprintf(&step_fname, "step_2d/out%i_%i.step", rule, i);
        out_step_file = fopen(step_fname, "w+");
        fprintf(out_step_file, "%s", out_string);

        free(out_string);
        fclose(out_step_file);
      }
    }
  }

  for (int i = 0; i < N; i++) {
    free(placeholder[i]);
  }
  free(placeholder);
  fclose(out_file);
}

void process_rule(int rule, int save_steps)
{
  fflush(stdout);
  write_to_file(rule, save_steps);
}

int main(int argc, char** argv)
{
  /* Write steps for a given rule  */
  if (argc == 2) {
    int rule = atoi(argv[1]);
    process_rule(rule, 1);
    return 0;
  }

  /* Generate compression plots for many rules */
  time_t t;
  srand((unsigned)time(&t));

  for (int i = 0; i < 400; i++) {
    int random_rule = rand();
    while (random_rule >= 1 << 18) {
      random_rule = rand();
    }
    random_rule = i + 122100;

    printf("\r%i: Rule %i", i, random_rule);
    process_rule(random_rule, 0);
  }
}
