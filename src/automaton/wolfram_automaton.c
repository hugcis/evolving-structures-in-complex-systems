#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include "wolfram_automaton.h"
#include "utils/compress.h"
#include "utils/utils.h"

#define NEIGHBORS 3


void print_bits_spaced(int bits_to_read, uint8_t a[bits_to_read],
                       char buf[bits_to_read + 1])
{
  for (int i = 0; i < bits_to_read; i++) {
    buf[i] = a[i] + '0';
  }
  buf[bits_to_read] = '\0';
}


void update_step(size_t size, size_t rule_size, uint8_t base[size],
                 uint8_t placeholder[size],
                 uint8_t rule[rule_size],
                 int states, int radius)
{
  memcpy(placeholder, base, size * sizeof(uint8_t));

  for (size_t r = 0; r < size; r++) {
    int c = 0;
    for (int w = -radius; w <= radius; w++) {
      if (base[(r + w + size) % size] == 1) {
        c += ipow(states, 2 * radius - (w + radius));
      }
    }
    placeholder[r] = rule[c];
  }

  memcpy(base, placeholder, size * sizeof(uint8_t));
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
    for (size_t i = 0; i < size; i++) {
      a[i] = rand() % states;
    }
  } else if (mode == RAND_SMALL) {
    for (size_t i = size / 2 - size / 5; i < size / 2 + size / 5; i++) {
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

  char output_string[size + 1];
  print_bits_spaced(size, A, output_string);
  fputs(output_string, steps_file);
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

  char* dbl_ouput = (char*) malloc(( 2 * size + 1 ) * sizeof(char));

  uint8_t A[size];
  uint8_t placeholder[size];
  uint8_t final[size];
  char spaced_output[size + 1];
  char final_output[size + 1];

  init_automat(states, size, A, options->init);
  memcpy(final, A, size * sizeof(uint8_t));

  size_t v;
  int end_comp_size = 0, comp_size = 0, dbl_comp_size = 0;

  for (v = 0; v < steps; ++v) {
    update_step(size, rule_size, final, placeholder, rule,
                states, options->radius);
  }
  print_bits_spaced(size, final, final_output);
  end_comp_size = compress_memory_size(final_output, size);

  for (v = 0; v < steps; ++v) {

    print_bits_spaced(size, A, spaced_output);
    fprintf(out_steps_file, "%s\n", spaced_output);

    if (print_automaton == 1) {
      printf("%s\n", spaced_output);
    }

    update_step(size, rule_size, A, placeholder,
                rule, states, options->radius);

    if (v % options->grain == 0) {
      comp_size = compress_memory_size(spaced_output, size);

      fprintf(out_file, "%lu    %i    %i\n", v, comp_size,
              dbl_comp_size);
      sprintf(dbl_ouput, "%s%s", spaced_output, final_output);
      dbl_comp_size = compress_memory_size(dbl_ouput, 2 * size);

      if (options->write == WRITE_STEP) {
        write_step(size, rule_size, A, rule, v, states);
      }

      if (v == options->grain) {
        printf("%i\t", comp_size);
      }
    }
  }
  printf("%i\t%i\t%i", end_comp_size, comp_size, dbl_comp_size);
  fclose(out_file);
  fclose(out_steps_file);
}
