#include "compress.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BITS 800
#define STEPS 800

enum InitMode { ONE,
  RANDOM,
  RAND_SMALL };

void set_bit(int A[], int k)
{
  A[k / 32] |= 1 << (k % 32); // Set the bit at the k-th position in A[i]
}

void clear_bit(int A[], int k)
{
  A[k / 32] &= ~(1 << (k % 32));
}

int test_bit(int A[], int k)
{
  return ((A[k / 32] & (1 << (k % 32))) != 0);
}

char* print_bits_spaced(int a[], int bits_to_read, int spaced)
{
  int n_bits = spaced == 1 ? 2 * bits_to_read - 1 : bits_to_read;
  int buf_index;

  char* buf = (char*)malloc((n_bits + 1) * sizeof(char));

  for (int i = 0; i < bits_to_read; i++) {
    buf_index = spaced == 1 ? 2 * i : i;
    buf[buf_index] = (a[i / 32] & (1 << (i % 32))) != 0 ? '1' : '0';
    if (spaced == 1 & buf_index + 1 < n_bits) {
      buf[buf_index + 1] = ' ';
    }
  }
  buf[n_bits] = '\0';
  return buf;
}

char* printBits(int a[], int bits_to_read)
{
  return print_bits_spaced(a, bits_to_read, 0);
}

void update_step(int base[], int len, int rule)
{
  int size = len / 32 + 1;
  int* A = (int*)malloc(size * sizeof(int));

  memcpy(A, base, size * sizeof(int));

  for (int r = 0; r < len; r++) {
    int c = 0;
    if (test_bit(base, (r - 1 + len) % len)) {
      c |= (1 << 2);
    }
    if (test_bit(base, r)) {
      c |= (1 << 1);
    }
    if (test_bit(base, (r + 1) % len)) {
      c |= (1 << 0);
    }
    if (rule & (1 << c)) {
      set_bit(A, r);
    } else {
      clear_bit(A, r);
    }
  }

  memcpy(base, A, size * sizeof(int));
  free(A);
}

void init_automat(int a[], int size, enum InitMode mode)
{
  if (mode == ONE) {
    set_bit(a, size / 2);
  } else if (mode == RANDOM) {
    time_t t;
    srand((unsigned)time(&t));
    for (int i = 0; i < size; i++) {
      if (rand() % 2 == 0) {
        set_bit(a, i);
      }
    }
  } else if (mode == RAND_SMALL) {
    time_t t;
    srand((unsigned)time(&t));
    set_bit(a, 0);
    for (int i = 1; i < size / 5; i++) {
      if (rand() % 2 == 0) {
        set_bit(a, i);
      }
    }
  }
}

void write_step(int A[], int rule, int step)
{
  FILE* steps_file;
  char* steps_fname;
  asprintf(&steps_fname, "steps/out%i_%i.step", rule, step);
  steps_file = fopen(steps_fname, "w+");

  char* output_string = printBits(A, BITS);
  fputs(output_string, steps_file);
  free(output_string);
  fclose(steps_file);
}

void write_to_file(int rule, int print_automaton, int write)
{
  FILE* out_file;
  char* fname;
  asprintf(&fname, "data/out%i.dat", rule);
  out_file = fopen(fname, "w+");

  FILE* out_steps_file;
  char* out_steps_fname;
  asprintf(&out_steps_fname, "steps/out%i.steps", rule);
  out_steps_file = fopen(out_steps_fname, "w+");

  int A[BITS / 32 + 1] = {};
  init_automat(A, BITS, RAND_SMALL);

  for (int i = 0; i < STEPS; i++) {

    if (i % 50 == 0 & write == 1) {
      write_step(A, rule, i);
    }

    char* spaced_output = print_bits_spaced(A, BITS, 1);
    fprintf(out_steps_file, "%s\n", spaced_output);
    free(spaced_output);

    if (print_automaton == 1) {
      printf("%s\n", printBits(A, BITS));
    }

    update_step(A, BITS, rule);

    if (i % 30 == 0) {
      char* output_string = printBits(A, BITS);
      fprintf(out_file, "%i    %i\n", i,
          compress_memory_size(output_string, BITS));
      free(output_string);
    }
  }
  fclose(out_file);
  fclose(out_steps_file);
}

int main()
{
  for (int rule = 0; rule < 256; rule++) {
    printf("\rRule %i", rule);
    write_to_file(rule, 0, 1);
  }
  printf("\n");
  return 0;
}
