#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "compress.h"

#define N 256
#define STEPS 1000
#define GRAIN 50
#define STATES 2
#define HORIZON (1)
#define SIDE ( (2 * HORIZON) + 1 )
#define NEIGH_SIZE ( (SIDE * SIDE) - 1 )
#define NEIGHBOR ( (STATES - 1) * NEIGH_SIZE )
#define RULE_SIZE ( (NEIGHBOR + 1) * STATES )
#define GRULE_SIZE  ( (int)pow(STATES, NEIGH_SIZE + 1) )

typedef void (*ProcessF)(uint8_t[][N] , size_t,
                         uint8_t[] , uint8_t**);

unsigned long hash(char *str)
{
  unsigned long hash = 5381;
  int c;

  while (( c = *str++ )) {
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  }

  return hash;
}


char* print_bits(uint8_t a[][N], size_t bits_dim1, size_t bits_dim2)
{
  char* buf = (char*)malloc(((bits_dim1 + 1) * bits_dim2 + 1)
      * sizeof(uint8_t));
  char out_string = '0';

  for (size_t i = 0; i < bits_dim1; i++) {
    for (size_t j = 0; j < bits_dim2; j++) {
      out_string = '0' + a[i][j];
      buf[i * (bits_dim1 + 1) + j] = out_string;
    }
    buf[i * (bits_dim1 + 1) + bits_dim2] = '\n';
  }
  buf[(bits_dim1 + 1) * bits_dim2] = '\0';
  return buf;
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


void update_step_general(uint8_t base[][N], size_t size,
                         uint8_t rule[GRULE_SIZE], uint8_t** A)
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
          position += current_value * (1 << increment);
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

void update_step_totalistic(uint8_t base[][N], size_t size,
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

void process_rule(uint8_t rule[GRULE_SIZE], char rule_buf[],
                  int save_steps, int joint_complexity,
                  int totalistic)
{
  FILE* out_file;
  char* fname;
  int last_compressed_size;
  int compressed_size;
  int size_sum;
  float ratio = 0;

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

  char* last_step[(N + 1) * N + 1];
  char* double_placeholder[2 * ((N + 1) * N + 1)];
  char* out_string;
  int flag = 0;

  for (int i = 0; i < STEPS; i++) {
    process_function(A, N, rule, placeholder);

    if (i % GRAIN == 0) {
      last_compressed_size = compressed_size;
      out_string = print_bits(A, N, N);
      compressed_size = compress_memory_size(out_string, (N + 1) * N);

      if (compressed_size == last_compressed_size && flag == 1) {
        printf("\n");
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

        memcpy(last_step, out_string, (N + 1) * N + 1);
        memcpy(double_placeholder, out_string, (N + 1) * N + 1);
        if (i > 0) {
          size_sum = last_compressed_size + compressed_size;
          ratio = (size_sum - double_compressed_size)/(float)size_sum;
        }
        fprintf(out_file, "%i    %i    %f\n", i, compressed_size, ratio);
      } else {
        fprintf(out_file, "%i    %i\n", i, compressed_size);
      }

      printf("%i  ", compressed_size);

      if (save_steps == 1) {
        FILE* out_step_file;
        char* step_fname;
        asprintf(&step_fname, "step_2d_%i/out%s_%i.step",
                 STATES, rule_buf, i);
        out_step_file = fopen(step_fname, "w+");
        fprintf(out_step_file, "%s", out_string);
        fclose(out_step_file);
      }

      free(out_string);
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
    rule_number += rule_array[s] * pow(STATES, s);
    if (STATES >= 3) {
      rule_buf[s] = '0' + rule_array[s];
    }
  }

  if (STATES == 2) {
    sprintf(rule_buf, "%lu", rule_number);
  }
}

void symmetrize_rule(uint8_t rule_array[GRULE_SIZE])
{
  uint32_t position_180;
  uint32_t position_90;
  uint32_t position_270;
  uint32_t position_vflip;
  /* Keep track of already seen positions */
  uint8_t book_keep[GRULE_SIZE];
  for (int i = 0; i < GRULE_SIZE; ++i) {
    book_keep[i] = 0;
  }

  int pos;

  for (int i = 0; i < GRULE_SIZE; ++i) {
    /* Skip already seen positions */
    /* if (book_keep[i] == 1) { */
      /* continue; */
    /* } */

    position_180 = 0;
    position_90 = 0;
    position_270 = 0;
    position_vflip = 0;
    for (int p = 0; p < NEIGH_SIZE + 1; ++p) {
      /* 180° rotation */
      position_180 += (uint32_t)pow(STATES, p) *
        ((i / (uint32_t)pow(STATES, NEIGH_SIZE  - p)) % STATES);

      /* 90° rotation */
      pos = (NEIGH_SIZE - SIDE + 1 - (SIDE * (p%SIDE)) + p/SIDE);
      position_90 += (uint32_t)pow(STATES, p) *
        ((i / (uint32_t)pow(STATES, pos)) % STATES);

      /* 270° rotation */
      position_270 += (uint32_t)pow(STATES, p) *
        ((i / (uint32_t)pow(STATES, NEIGH_SIZE  - pos)) % STATES);

      /* Vertical flip */
      pos = (SIDE * (p / SIDE)) + (SIDE - 1 - (p % 3));
      if (i==99){
        printf("%i", pos);
        fflush(stdout);
      }
      position_vflip += (uint32_t)pow(STATES, p) *
        ((i / (uint32_t)pow(STATES, pos)) % STATES);

      /* Horizontal flip */
    }

    rule_array[position_180] = rule_array[i];
    rule_array[position_90] = rule_array[i];
    rule_array[position_270] = rule_array[i];
    rule_array[position_vflip] = rule_array[i];
  }

}

void build_rule_from_args(uint8_t rule_array[RULE_SIZE],
                          char rule_buf[RULE_SIZE + 1],
                          char** argv)
{
  unsigned long rule_number = 0UL;
  if (STATES == 2) {
    char* eptr;
    rule_number = strtoul(argv[1], &eptr, 10);
    for (int s = 0 ; s < RULE_SIZE ; ++s) {
      rule_array[s] = (uint8_t)((rule_number / (int)pow(STATES, s)) % STATES);
    }

    const int n = snprintf(NULL, 0, "%lu", rule_number);
    assert(n > 0);
    int c = snprintf(rule_buf, n+1, "%lu", rule_number);
    assert(rule_buf[n] == '\0');
    assert(c == n);

  } else {
    for (int s = 0 ; s < RULE_SIZE ; ++s) {
      rule_array[s] = argv[1][s] - '0';
      printf("%"PRIu8, rule_array[s]);
      rule_buf[s] = argv[1][s];
    }
    rule_buf[RULE_SIZE] = '\0';
  }
}

void generate_general_rule(uint8_t rule_array[GRULE_SIZE],
                           char rule_buf[GRULE_SIZE + 1])
{
  for (int v = 0; v < GRULE_SIZE; v++) {
    rule_array[v] = (uint8_t)(rand() % STATES);
    sprintf(&rule_buf[v], "%"PRIu8, rule_array[v]);
  }

  symmetrize_rule(rule_array);

  sprintf(rule_buf, "%lu", hash(rule_buf));
}

int main(int argc, char** argv)
{
  char rule_buf[GRULE_SIZE + 1];
  uint8_t rule_array[GRULE_SIZE];

  /* Write steps for a given rule  */
  if (argc == 2) {
    build_rule_from_args(rule_array, rule_buf, argv);

    process_rule(rule_array, rule_buf, 1, 1, 1);
    return 0;
  }


  /* Generate compression plots for many rules */
  time_t t;
  srand((unsigned)time(&t));

  fflush(stdout);
  for (int i = 0; i < 1; ++i) {
    /* generate_totalistic_rule(rule_array, rule_buf); */
    for (int v = 0; v < GRULE_SIZE; v++) {
      rule_array[v] = (uint8_t)(rand() % STATES);
      sprintf(&rule_buf[v], "%"PRIu8, rule_array[v]);
    }

    symmetrize_rule(rule_array);
    assert(rule_array[99] = rule_array[270]);
    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("%i: Rule %s\n", i, rule_buf);
    fflush(stdout);

    /* process_rule(rule_array, rule_buf, 1, 1, 0); */
  }
  printf("\n");
  return 0;
}
