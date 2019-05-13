#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "2d_automaton.h"
#include "compress.h"

#define GRAIN 50
#define NEIGHBOR ( (STATES - 1) * NEIGH_SIZE )
#define RULE_SIZE ( (NEIGHBOR + 1) * STATES )

typedef void (*ProcessF)(uint64_t, uint8_t[][N] , size_t,
                         uint8_t[] , uint8_t**);

uint32_t ipow(int base, int exp)
{
    uint32_t result = 1;
    for (;;)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

unsigned long hash(char *str)
{
  unsigned long hash = 5381;
  int c;

  while (( c = *str++ )) {
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  }

  return hash;
}

/**
 * Format the 2d automaton as a string in buf.
 */
void print_bits(uint8_t a[][N], size_t dim1, size_t dim2, char* buf)
{
  char out_string = '0';

  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      out_string = '0' + a[i][j];
      buf[i * (dim1 + 1) + j] = out_string;
    }
    buf[i * (dim1 + 1) + dim2] = '\n';
  }
  buf[(dim1 + 1) * dim2] = '\0';
}

void print_aggregated_bits(size_t dim1, size_t dim2, uint8_t a[][dim2],
                           char* buf, int offset)
{
  int inc1 = 0;
  int inc2, counter;
  char out_string = '0';
  for (size_t i = 0; i < dim1 - offset + 1; i+=offset) {
    inc2 = 0;
    for (size_t j = 0; j < dim2 - offset + 1; j+=offset) {
      counter = 0;
      for (int c1 = 0; c1 < offset; ++c1) {
        for (int c2 = 0; c2 < offset; ++c2) {
          counter += a[(i + c1) % dim1][(j + c2) % dim2];
        }
      }
      out_string = '0' + counter;
      buf[inc1 * (dim1/offset + 1) + inc2] = out_string;
      inc2++;
    }
    buf[inc1 * (dim1/offset + 1) + inc2] = '\n';
    inc1++;
  }
  buf[(inc1 - 1) * (dim1/offset + 1) + inc2 + 1] = '\0';
}

void init_automat(uint8_t a[][N], size_t size)
{
  assert(size > 20);
  time_t t;
  srand((unsigned)time(&t));
  for (size_t i = size / 2 - 10; i < size / 2 + 10; i++) {
    for (size_t j = size / 2 - 10; j < size / 2 + 10; j++) {
      a[i % size][j % size] = rand() % STATES;
    }
  }
}

void update_step_general(uint64_t grule_size, uint8_t base[][N], size_t size,
                         uint8_t rule[grule_size], uint8_t** A)
{
  uint64_t position;
  int current_value;
  int increment;

  for (size_t i = 0; i < size; i++) {
    memcpy(A[i], base[i], size * sizeof(uint8_t));

    for (size_t j = 0; j < size; j++) {
      position = 0;
      increment = 0;
      for (int k = - HORIZON; k <= HORIZON; k++) {
        for (int l = - HORIZON; l <= HORIZON; l++) {
          current_value = base[(i + k + size) % size][(j + l + size) % size];
          position += current_value * ipow(STATES, increment);
          ++increment;
        }
      }
      A[i][j] = rule[position];
    }
  }
  for (size_t i = 0; i < size; i++) {
    memcpy(base[i], A[i], size * sizeof(uint8_t));
  }
}

void update_step_totalistic(uint64_t rule_size, uint8_t base[][N], size_t size,
                            uint8_t rule[RULE_SIZE], uint8_t** A)
{
  int count;
  for (size_t i = 0; i < size; i++) {
    memcpy(A[i], base[i], size * sizeof(uint8_t));

    for (size_t j = 0; j < size; j++) {
      count = 0;
      for (int k = - HORIZON; k <= HORIZON; k++) {
        for (int l = - HORIZON; l <= HORIZON; l++) {
          if (k != 0 | l != 0) {
            count += base[(i + k + size) % size][(j + l + size) % size];
          }
        }
      }
      A[i][j] = rule[STATES * count + base[i][j]];
    }
  }
  for (size_t i = 0; i < size; i++) {
    memcpy(base[i], A[i], size * sizeof(uint8_t));
  }
}

int count_cells(uint8_t A[N][N], size_t size)
{
  assert(A && size);
  int* counts = (int*)calloc(STATES, sizeof(int));
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      counts[A[i][j]] += 1;
    }
  }
  int value = counts[0];
  for (int i = 0; i < STATES; ++i) {
    if (counts[i] < value) {
      value = counts[i];
    }
  }
  return value + 1;
}

/**
 * Given a rule, create and simulate the corresponding automaton.
 */
void process_rule(uint64_t grule_size, uint8_t rule[grule_size],
                  char rule_buf[], int save_steps, int joint_complexity,
                  int totalistic, long steps)
{
  FILE* out_file;
  char* fname;
  int last_compressed_size;
  int compressed_size;
  /* int agg_compressed_size; */
  int last_cell_count;
  int cell_count = 0;
  int size_sum;
  float ratio = 0;
  float ratio2 = 0;

  ProcessF process_function = totalistic == 1 ?
    update_step_totalistic: update_step_general;

  asprintf(&fname, "data_2d_%i/out%s.dat", STATES, rule_buf);
  out_file = fopen(fname, "w+");

  uint8_t A[N][N] = { {} };
  init_automat(A, N);

  uint8_t** placeholder = (uint8_t**)malloc(N * sizeof(uint8_t*));
  for (int i = 0; i < N; i++) {
    placeholder[i] = (uint8_t*)malloc(N * sizeof(uint8_t));
  }

  char double_placeholder[2 * ((N + 1) * N + 1)];
  char out_string[(N+1) * N + 1];
  int flag = 0;

  for (int i = 0; i < steps; i++) {
    process_function(grule_size, A, N, rule, placeholder);

    if (i % GRAIN == 0) {
      last_compressed_size = compressed_size;
      last_cell_count = cell_count;

      /* print_aggregated_bits(N, N, A, out_string, 2); */
      /* agg_compressed_size = compress_memory_size(out_string, (N/2 + 1) * N); */

      print_bits(A, N, N, out_string);
      compressed_size = compress_memory_size(out_string, (N + 1) * N);
      cell_count = count_cells(A, N);

      /* Check if state has evolved from last time */
      if (compressed_size == last_compressed_size && flag == 1) {
        printf("\n");
        /* Cleanup before returning */
        for (int i = 0; i < N; i++) {
          free(placeholder[i]);
        }
        free(placeholder);
        fclose(out_file);
        return;
      } else if (compressed_size == last_compressed_size) {
        flag = 1;
      }

      if (joint_complexity == 1) {
        memcpy(&double_placeholder[(N + 1) * N + 1],
               out_string, (N + 1) * N + 1);

        int double_compressed_size =
          compress_memory_size(double_placeholder, 2 * ((N + 1) * N + 1));

        memcpy(double_placeholder, out_string, (N + 1) * N + 1);

        if (i > 0) {
          size_sum = last_compressed_size + compressed_size;
          ratio2 = (size_sum - double_compressed_size)/(float)size_sum;
          ratio = ( (last_compressed_size / (float)last_cell_count) +
                    (compressed_size / (float)cell_count) ) /
            (double_compressed_size / (float)(last_cell_count + cell_count));
        }
        fprintf(out_file, "%i    %i    %f    %f    %i    %i    %i\n",
                i, compressed_size, ratio, ratio2, cell_count, last_cell_count,
                double_compressed_size);
      } else {
        fprintf(out_file, "%i    %i\n", i, compressed_size);
      }

      printf("%i  ", compressed_size);
      fflush(stdout);

      if (save_steps == 1) {
        FILE* out_step_file;
        char* step_fname;
        asprintf(&step_fname, "step_2d_%i/out%s_%i.step",
                 STATES, rule_buf, i);
        out_step_file = fopen(step_fname, "w+");
        fprintf(out_step_file, "%s", out_string);
        fclose(out_step_file);
      }
    }
  }
  printf("\n");

  for (int i = 0; i < N; i++) {
    free(placeholder[i]);
  }
  free(placeholder);
  fclose(out_file);
}

void generate_totalistic_rule(uint8_t rule_array[RULE_SIZE],
                              char rule_buf[RULE_SIZE + 1])
{
  unsigned long rule_number = 0UL;

  for (int s = 0 ; s < RULE_SIZE ; ++s) {
    rule_array[s] = rand() % STATES;
    rule_number += rule_array[s] * ipow(STATES, s);
    if (STATES >= 3) {
      rule_buf[s] = '0' + rule_array[s];
    }
  }

  if (STATES == 2) {
    sprintf(rule_buf, "%lu", rule_number);
  }
}

/**
 * Symmetrize the rule by setting all the states and their symmetries to having
 * the same output.
 */
void symmetrize_rule(uint64_t grule_size, uint8_t rule_array[grule_size])
{
  uint32_t position_180;
  uint32_t position_90;
  uint32_t position_270;
  uint32_t position_vflip;
  uint32_t position_hflip;
  uint32_t position_dflip;
  uint32_t position_aflip;

  /* Keep track of already seen positions with book-keeping */
  uint8_t book_keep[grule_size];
  for (int i = 0; i < grule_size; ++i) {
    book_keep[i] = 0;
  }

  int pos;

  for (int i = 0; i < grule_size; ++i) {
    /* Skip already seen positions when looping through the rule */
    /* if (book_keep[i] == 1) { */
      /* continue; */
    /* } */

    position_180 = 0;
    position_90 = 0;
    position_270 = 0;
    position_vflip = 0;
    position_hflip = 0;
    position_dflip = 0;
    position_aflip = 0;

    /* Create the representation of the symmetrized position by swapping the
       states in its number representation. */
    for (int p = 0; p < NEIGH_SIZE + 1; ++p) {
      /* 180° rotation */
      position_180 += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, NEIGH_SIZE  - p)) % STATES);

      /* 90° rotation */
      pos = (NEIGH_SIZE - SIDE + 1 - (SIDE * (p%SIDE)) + p/SIDE);
      position_90 += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, pos)) % STATES);

      /* 270° rotation */
      position_270 += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, NEIGH_SIZE  - pos)) % STATES);

      /* Vertical flip */
      pos = (SIDE * (p / SIDE)) + (SIDE - 1 - (p % 3));
      position_vflip += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, pos)) % STATES);

      /* Horizontal flip */
      if (p/SIDE < SIDE/2) {
        pos = (p - SIDE + NEIGH_SIZE + 1)%(NEIGH_SIZE + 1);
      } else if (p/SIDE > SIDE/2) {
        pos = (p + SIDE)%(NEIGH_SIZE + 1);
      } else {
        pos = p;
      }
      position_hflip += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, pos)) % STATES);

      /* Diagonal flip */
      pos = NEIGH_SIZE - ((p % SIDE) * SIDE  + (p / SIDE));
      position_dflip += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, pos)) % STATES);

      /* Antidiagonal flip */
      pos = ((p % SIDE) * SIDE  + (p / SIDE));
      position_aflip += (uint32_t)ipow(STATES, p) *
        ((i / (uint32_t)ipow(STATES, pos)) % STATES);

    }

    /* Add all seen positions to the book to not process them again */
    book_keep[i] = 1;
    book_keep[position_180] = 1;
    book_keep[position_90] = 1;
    book_keep[position_270] = 1;
    book_keep[position_vflip] = 1;
    book_keep[position_hflip] = 1;
    book_keep[position_dflip] = 1;
    book_keep[position_aflip] = 1;

    rule_array[position_180] = rule_array[i];
    rule_array[position_90] = rule_array[i];
    rule_array[position_270] = rule_array[i];
    rule_array[position_vflip] = rule_array[i];
    rule_array[position_hflip] = rule_array[i];
    rule_array[position_dflip] = rule_array[i];
    rule_array[position_aflip] = rule_array[i];
  }

}

/**
 * Build a rule from the provided command-line arguments. Depending on the
 * number of states, it either expect a base-( STATES - 1 ) input or a base-10
 * input representing the rule.
 */
void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char* rule_arg)
{
  /* unsigned long rule_number = 0UL; */
  /* if (STATES == 2) { */
  /*   char* eptr; */
  /*   rule_number = strtoul(argv[1], &eptr, 10); */
  /*   for (int s = 0 ; s < RULE_SIZE ; ++s) { */
  /*     rule_array[s] = (uint8_t)((rule_number / (int)ipow(STATES, s)) % STATES); */
  /*   } */

  /*   const int n = snprintf(NULL, 0, "%lu", rule_number); */
  /*   assert(n > 0); */
  /*   int c = snprintf(rule_buf, n+1, "%lu", rule_number); */
  /*   assert(rule_buf[n] == '\0'); */
  /*   assert(c == n); */

  /* } else { */
  /*   for (int s = 0 ; s < RULE_SIZE ; ++s) { */
  /*     rule_array[s] = argv[1][s] - '0'; */
  /*     printf("%"PRIu8, rule_array[s]); */
  /*     rule_buf[s] = argv[1][s]; */
  /*   } */
  /*   rule_buf[RULE_SIZE] = '\0'; */
  /* } */

  /* Rule is given in base-(STATES - 1) format */
  for (int s = 0 ; s < grule_size; ++s) {
    rule_array[s] = rule_arg[s] - '0';
    rule_buf[s] = rule_arg[s];
  }
  rule_buf[grule_size] = '\0';
}

/**
 * Generate a random general rule.
 */
void generate_general_rule(uint64_t grule_size,
                           uint8_t rule_array[grule_size],
                           char rule_buf[grule_size + 1])
{
  int ur = 2 + rand()%7;
  for (int v = 0; v < grule_size; v++) {
    if (rand() % 10 < ur) {
      rule_array[v] = 1;
    } else {
      rule_array[v] = (uint8_t)(rand() % STATES);
    }
    /* rule_array[v] = (uint8_t)(rand() % STATES); */
  }

  symmetrize_rule(grule_size, rule_array);

  for (int v = 0; v < grule_size; v++) {
    sprintf(&rule_buf[v], "%"PRIu8, rule_array[v]);
  }
}
