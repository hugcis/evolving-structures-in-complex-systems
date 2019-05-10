#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include "wolfram_automaton.h"
#include "compress.h"

#define GRAIN 50
#define NEIGHBORS 3
#define STATES 2


char* print_bits_spaced(int bits_to_read, uint8_t a[bits_to_read])
{
  char* buf = (char*)malloc((bits_to_read + 1) * sizeof(char));

  for (int i = 0; i < bits_to_read; i++) {
    buf[i] = a[i] != 0 ? '1' : '0';
  }
  buf[bits_to_read] = '\0';
  return buf;
}


void update_step(size_t size, size_t rule_size, uint8_t base[size],
                 uint8_t rule[rule_size])
{
  uint8_t* A = (uint8_t*)malloc(size * sizeof(uint8_t));

  memcpy(A, base, size * sizeof(uint8_t));

  for (int r = 0; r < size; r++) {
    int c = 0;
    if (base[(r - 1 + size) % size] == 1) {
      c += (1 << 2);
    }
    if (base[r] == 1) {
      c += (1 << 1);
    }
    if (base[(r + 1) % size] == 1) {
      c += (1 << 0);
    }
    if (rule[c]) {
      A[r] = 1;
    } else {
      A[r] = 0;
    }
  }

  memcpy(base, A, size * sizeof(uint8_t));
  free(A);
}

static void init_automat(size_t size, uint8_t a[size], enum InitMode mode)
{
  memset(a, 0, sizeof(uint8_t) * size);
  if (mode == ONE) {
    a[size/2] = 1;
    return;
  }
  if (mode == RANDOM) {
    time_t t;
    srand((unsigned)time(&t));
    for (int i = 0; i < size; i++) {
      if (rand() % 2 == 0) {
        a[i] = 1;
      }
    }
  } else if (mode == RAND_SMALL) {
    time_t t;
    srand((unsigned)time(&t));
    for (int i = 1; i < size / 5; i++) {
      if (rand() % 2 == 0) {
        a[i] = 1;
      }
    }
  }
}

unsigned long rule_number(size_t rule_size, uint8_t rule[rule_size])
{
  unsigned long count = 0;
  unsigned long mult = 1;
  for (int i = 0; i < rule_size; ++i) {
    count += rule[i] * mult;
    mult *= STATES;
  }
  return count;
}

void write_step(size_t size, size_t rule_size, uint8_t A[size],
                uint8_t rule[rule_size], int step)
{
  FILE* steps_file;
  char* steps_fname;
  asprintf(&steps_fname, "steps/out%lu_%i.step",
           rule_number(rule_size, rule), step);

  steps_file = fopen(steps_fname, "w+");

  char* output_string = print_bits_spaced(size, A);
  fputs(output_string, steps_file);
  free(output_string);
  fclose(steps_file);
}

void write_to_file(size_t size, size_t rule_size, uint8_t rule[rule_size],
                   int print_automaton, int write, struct Options1D* options)
{
  size_t steps = options->timesteps;

  FILE* out_file;
  char* fname;
  asprintf(&fname, "data/out%lu.dat", rule_number(rule_size, rule));
  out_file = fopen(fname, "w+");

  FILE* out_steps_file;
  char* out_steps_fname;
  asprintf(&out_steps_fname, "steps/out%lu.steps",
           rule_number(rule_size, rule));
  out_steps_file = fopen(out_steps_fname, "w+");

  uint8_t A[size];
  init_automat(size, A, options->init);
  int v, comp_size = 0;

  for (v = 0; v < steps; ++v) {

    char* spaced_output = print_bits_spaced(size, A);
    fprintf(out_steps_file, "%s\n", spaced_output);
    free(spaced_output);

    if (print_automaton == 1) {
      printf("%s\n", print_bits_spaced(size, A));
    }

    update_step(size, rule_size, A, rule);

    if (v % GRAIN == 0) {
      char* output_string = print_bits_spaced(size, A);
      comp_size = compress_memory_size(output_string, size);

      fprintf(out_file, "%i    %i    %lu\n", v, comp_size,
              rule_number(rule_size, rule));
      free(output_string);
    }

  }
  printf("%i", comp_size);
  fclose(out_file);
  fclose(out_steps_file);
}
