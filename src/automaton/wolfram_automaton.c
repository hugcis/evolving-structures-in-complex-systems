#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include "wolfram_automaton.h"
#include "compress.h"
#include "utils.h"

#define NEIGHBORS 3


char* print_bits_spaced(int bits_to_read, uint8_t a[bits_to_read])
{
  char* buf = (char*)malloc((bits_to_read + 1) * sizeof(char));

  for (int i = 0; i < bits_to_read; i++) {
    buf[i] = a[i] + '0';
  }
  buf[bits_to_read] = '\0';
  return buf;
}


void update_step(size_t size, size_t rule_size, uint8_t base[size],
                 uint8_t rule[rule_size], int states)
{
  uint8_t* A = (uint8_t*)malloc(size * sizeof(uint8_t));

  memcpy(A, base, size * sizeof(uint8_t));

  for (size_t r = 0; r < size; r++) {
    int c = 0;
    if (base[(r - 1 + size) % size] == 1) {
      c += ipow(states, 2);
    }
    if (base[r] == 1) {
      c += ipow(states, 1);
    }
    if (base[(r + 1) % size] == 1) {
      c += ipow(states, 0);
    }
    A[r] = rule[c];
  }

  memcpy(base, A, size * sizeof(uint8_t));
  free(A);
}

static void init_automat(int states, size_t size, uint8_t a[size],
                         enum InitMode mode)
{
  memset(a, 0, sizeof(uint8_t) * size);
  if (mode == ONE) {
    a[size/2] = 1;
    return;
  }
  if (mode == RANDOM) {
    time_t t;
    srand((unsigned)time(&t));
    for (size_t i = 0; i < size; i++) {
      a[i] = rand() % states;
    }
  } else if (mode == RAND_SMALL) {
    time_t t;
    srand((unsigned)time(&t));
    for (size_t i = 1; i < size / 5; i++) {
      a[i] = rand() % states;
    }
  }
}

unsigned long rule_number(int states, size_t rule_size,
                          uint8_t rule[rule_size])
{
  unsigned long count = 0;
  unsigned long mult = 1;
  for (size_t i = 0; i < rule_size; ++i) {
    count += rule[i] * mult;
    mult *= states;
  }
  return count;
}

void write_step(size_t size, size_t rule_size, uint8_t A[size],
                uint8_t rule[rule_size], int step, int states)
{
  FILE* steps_file;
  char* steps_fname;
  asprintf(&steps_fname, "steps/out%lu_%i.step",
           rule_number(states, rule_size, rule), step);

  steps_file = fopen(steps_fname, "w+");

  char* output_string = print_bits_spaced(size, A);
  fputs(output_string, steps_file);
  free(output_string);
  fclose(steps_file);
}

void write_to_file(size_t size, size_t rule_size, uint8_t rule[rule_size],
                   int print_automaton, struct Options1D* options,
                   int states)
{
  size_t steps = options->timesteps;

  FILE* out_file;
  char* fname;
  asprintf(&fname, "data/out%lu.dat", rule_number(states, rule_size, rule));
  out_file = fopen(fname, "w+");

  FILE* out_steps_file;
  char* out_steps_fname;
  asprintf(&out_steps_fname, "steps/out%lu.steps",
           rule_number(states, rule_size, rule));
  out_steps_file = fopen(out_steps_fname, "w+");

  char* spaced_output;
  char* final_output;
  char* dbl_ouput = (char*) malloc(( 2 * size + 1 ) * sizeof(char));

  uint8_t A[size];
  uint8_t final[size];
  init_automat(states, size, A, options->init);
  memcpy(final, A, size * sizeof(uint8_t));

  size_t v;
  int comp_size = 0, dbl_comp_size = 0;

  for (v = 0; v < steps; ++v) {
    update_step(size, rule_size, final, rule, states);
  }
  final_output = print_bits_spaced(size, A);

  for (v = 0; v < steps; ++v) {

    spaced_output = print_bits_spaced(size, A);
    fprintf(out_steps_file, "%s\n", spaced_output);

    if (print_automaton == 1) {
      printf("%s\n", print_bits_spaced(size, A));
    }

    update_step(size, rule_size, A, rule, states);

    if (v % options->grain == 0) {
      comp_size = compress_memory_size(spaced_output, size);

      fprintf(out_file, "%lu    %i    %i\n", v, comp_size,
              dbl_comp_size);
      sprintf(dbl_ouput, "%s%s", spaced_output, final_output);
      dbl_comp_size = compress_memory_size(dbl_ouput, 2 * size);

      if (options->write == WRITE_STEP) {
        write_step(size, rule_size, A, rule, v, states);
      }
    }
    free(spaced_output);

  }
  printf("%i\t%i", comp_size, dbl_comp_size);
  free(final_output);
  fclose(out_file);
  fclose(out_steps_file);
}
