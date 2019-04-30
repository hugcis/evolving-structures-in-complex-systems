#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include "compress.h"


#define SIZE 50
#define STEPS 800
#define GRAIN 50
#define NEIGHBORS 3
#define STATES 2
#define RULE_SIZE (int)pow(STATES, NEIGHBORS)
#define RUNS 2000

enum InitMode { ONE,
  RANDOM,
  RAND_SMALL };


char* print_bits_spaced(uint8_t a[SIZE], int bits_to_read)
{
  char* buf = (char*)malloc((bits_to_read + 1) * sizeof(char));

  for (int i = 0; i < bits_to_read; i++) {
    buf[i] = a[i] != 0 ? '1' : '0';
  }
  buf[bits_to_read] = '\0';
  return buf;
}


void update_step(uint8_t base[SIZE], size_t size, uint8_t rule[RULE_SIZE])
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

void init_automat(uint8_t a[], int size, enum InitMode mode)
{
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

unsigned long rule_number(uint8_t rule[RULE_SIZE], size_t size)
{
  unsigned long count = 0;
  unsigned long mult = 1;
  for (int i = 0; i < size; ++i) {
    count += rule[i] * mult;
    mult *= STATES;
  }
  return count;
}

void write_step(uint8_t A[SIZE], uint8_t rule[RULE_SIZE], int step)
{
  FILE* steps_file;
  char* steps_fname;
  asprintf(&steps_fname, "steps/out%lu_%i.step",
           rule_number(rule, RULE_SIZE), step);

  steps_file = fopen(steps_fname, "w+");

  char* output_string = print_bits_spaced(A, SIZE);
  fputs(output_string, steps_file);
  free(output_string);
  fclose(steps_file);
}

uint32_t write_to_file(uint8_t rule[RULE_SIZE], int print_automaton, int write)
{
  uint8_t new_rule[RULE_SIZE];
  unsigned long previous_rule = 0;

  FILE* out_file;
  char* fname;
  asprintf(&fname, "data/out%lu.dat", rule_number(rule, RULE_SIZE));
  out_file = fopen(fname, "w+");

  FILE* out_steps_file;
  char* out_steps_fname;
  asprintf(&out_steps_fname, "steps/out%lu.steps",
           rule_number(rule, RULE_SIZE));
  out_steps_file = fopen(out_steps_fname, "w+");

  uint8_t A[SIZE] = {};
  init_automat(A, SIZE, RANDOM);
  int v;

  for (v = 0; v < 100; ++v) {
    for (int i = 0; i < RULE_SIZE; i++) {

      /* if (i % GRAIN == 0 & write == 1) { */
        /* write_step(A, rule, i); */
      /* } */

      char* spaced_output = print_bits_spaced(A, SIZE);
      fprintf(out_steps_file, "%s\n", spaced_output);
      free(spaced_output);

      if (print_automaton == 1) {
        printf("%s\n", print_bits_spaced(A, SIZE));
      }

      update_step(A, SIZE, rule);
      new_rule[i] = A[SIZE/2];

      if (i % 30 == 0) {
        char* output_string = print_bits_spaced(A, SIZE);
        fprintf(out_file, "%i    %i    %lu\n", i,
                compress_memory_size(output_string, SIZE),
                rule_number(rule, RULE_SIZE));
        free(output_string);
      }
    }
    memcpy(rule, new_rule, RULE_SIZE * sizeof(uint8_t));
    /* printf("%lu\n", rule_number(rule, RULE_SIZE)); */

    if (rule_number(rule, RULE_SIZE) == previous_rule) {
      fclose(out_file);
      fclose(out_steps_file);
      return v;
    }

    previous_rule = rule_number(rule, RULE_SIZE);
  }
  fclose(out_file);
  fclose(out_steps_file);
  return v;
}

int main()
{
  uint8_t rule[RULE_SIZE];
  uint32_t length[256] = {};

  time_t t;
  srand((unsigned)time(&t));

  for (int c = 0; c < RUNS; ++c) {
    for (int n = 0; n < 256; n++) {
      for (int i = 0; i < RULE_SIZE; ++i) {
        rule[i] = (n & (1 << i)) ? 1: 0;
        /* printf("%"PRIu8, rule[i]); */
      }
      /* printf("\n"); */
      /* printf("%lu\n", rule_number(rule, RULE_SIZE)); */
      length[n] += write_to_file(rule, 0, 1);
    }
  }

  for (int n = 0; n < 256; n++) {
    printf("%i\t%f\n", n, ((float)length[n])/RUNS);
  }
  /* for (int rule = 0; rule < RULE_SIZE; rule++) { */
    /* printf("\rRule %i", rule); */
    /* write_to_file(rule, 0, 1); */
  /* } */
  return 0;
}
