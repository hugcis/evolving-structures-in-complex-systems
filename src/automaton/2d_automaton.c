#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "2d_automaton.h"
#include "rule.h"
#include "nn/nn.h"
#include "utils/compress.h"
#include "utils/utils.h"
#include "utils/hashmap.h"

#define KEY_MAX_LENGTH (256)
#define PERT 1E-15

typedef void (*ProcessF)(uint64_t, size_t size, uint8_t [size][size],
                         uint8_t[] , uint8_t [size][size], int, int);

typedef struct entrop_state_ph_s
{
  double entropy;
  uint32_t visited;
  int states;
} entrop_state_ph_t;

typedef struct cross_ent_struct_s
{
  map_t source_map;
  int states;
  double cross_entropy;
} cross_ent_struct_t;

typedef struct data_struct_s
{
  char key_string[KEY_MAX_LENGTH];
  double* number_array;
} data_struct_t;

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
void print_bits(size_t dim1, size_t dim2, uint8_t a[dim1][dim2], char* buf)
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

void init_automat(size_t size, uint8_t a[size][size], int states)
{
  assert(size > 20);
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      /* Initialize center 20x20 square to random */
      if ((i >= size/2 - 10) && (i < size/2 + 10) &&
          (j >= size/2 - 10) && (j < size/2 + 10)) {
        a[i % size][j % size] = (uint8_t)(rand() % states);
      }
      /* Rest is 0 */
      else {
        a[i][j] = (uint8_t)0;
      }
    }
  }
}

void update_step_general(uint64_t grule_size, size_t size,
                         uint8_t A[size][size],
                         uint8_t rule[grule_size],
                         uint8_t base[size][size], int states,
                         int horizon)
{
  uint64_t position;
  uint8_t current_value;
  int increment;

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      position = 0;
      increment = 0;
      for (int k = - horizon; k <= horizon; k++) {
        for (int l = - horizon; l <= horizon; l++) {
          current_value = base[(i + k + size) % size][(j + l + size) % size];
          position += current_value * ipow(states, increment);
          ++increment;
        }
      }

      A[i][j] = rule[position];
    }
  }
  memcpy(base, A, size * size * sizeof(uint8_t));
}

void update_step_totalistic(uint64_t rule_size, size_t size,
                            uint8_t base[size][size],
                            uint8_t rule[rule_size],
                            uint8_t A[size][size], int states,
                            int horizon)
{
  int count;
  for (size_t i = 0; i < size; i++) {
    memcpy(A[i], base[i], size * sizeof(uint8_t));

    for (size_t j = 0; j < size; j++) {
      count = 0;
      for (int k = - horizon; k <= horizon; k++) {
        for (int l = - horizon; l <= horizon; l++) {
          if (k != 0 | l != 0) {
            count += base[(i + k + size) % size][(j + l + size) % size];
          }
        }
      }
      A[i][j] = rule[states * count + base[i][j]];
    }
  }
  for (size_t i = 0; i < size; i++) {
    memcpy(base[i], A[i], size * sizeof(uint8_t));
  }
}

int count_cells(size_t size, uint8_t A[size][size], int states)
{
  int* counts = (int*)calloc(states, sizeof(int));
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      counts[A[i][j]] += 1;
    }
  }
  int value = counts[0];
  for (int i = 0; i < states; ++i) {
    if (counts[i] < value) {
      value = counts[i];
    }
  }
  free(counts);
  return value + 1;
}

int normalize_probs(any_t states, any_t data)
{
  data_struct_t* in = (data_struct_t*) data;
  int st = *(int*)states;
  double sum = 0;
  for (int i = 0; i < st; i++) {
    sum += in->number_array[i];
  }
  for (int i = 0; i < st; i++) {
    in->number_array[i] /= sum;
  }
  return MAP_OK;
}

void build_key_string(int key_len, char key[key_len],
                      size_t size, uint8_t automaton[size][size],
                      int offset, size_t i, size_t j)
{
  int key_index;
  for (int a = -offset; a < offset + 1; a++) {
    for (int b = -offset; b < offset + 1; b++) {
      key_index = (a + offset) * (2 * offset + 1) + (b + offset);
      if (a == 0 && b == 0) {
        key[key_index] = 'x';
      }
      else {
        key[key_index] =
          '0' + automaton[(i + a + size) % size][(j + b + size) % size];
      }
    }
  }
  key[key_len - 1] = '\0';
}

void populate_map(map_t* map, size_t size, uint8_t automaton[size][size],
                  int offset, int states)
{
  data_struct_t* value;
  int key_len = (2*offset + 1) * (2*offset + 1) + 1;
  char key[key_len];
  int error;

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      build_key_string(key_len, key, size, automaton, offset, i, j);
      error = hashmap_get(map, key, (void**)(&value));

      if (error==MAP_MISSING) {
        value = malloc(sizeof(data_struct_t));

        snprintf(value->key_string, key_len, "%s", key);
        value->number_array = (double*)calloc(states, sizeof(double));

        error = hashmap_put(map, value->key_string, value);
        assert(error==MAP_OK);
      }
      value->number_array[automaton[i][j]] += 1;
    }
  }
  hashmap_iterate(map, (PFany)&normalize_probs, (any_t)&states);
}

int reduce_entropy(any_t state_entropy, any_t in)
{
  entrop_state_ph_t* state_ent = (entrop_state_ph_t*) state_entropy;
  data_struct_t* data = (data_struct_t*) in;

  for (int i = 0; i < state_ent->states; i ++) {
    if (data->number_array[i] > 0) {
      state_ent->entropy += -log(data->number_array[i]) *
        (data->number_array[i]);
    }
  }
  state_ent->visited += 1;
  return MAP_OK;
}

entrop_state_ph_t* entropy_score(map_t map, int states)
{
  entrop_state_ph_t* ph = calloc(1, sizeof(entrop_state_ph_t));
  ph->states = states;
  hashmap_iterate(map, (PFany)&reduce_entropy, (any_t)ph);

  return ph;
}

double predictive_score(map_t map, int states, size_t size,
                        uint8_t automaton[size][size], int offset)
{
  data_struct_t* value;
  int key_len = (2*offset + 1) * (2*offset + 1) + 1;
  char key[key_len];
  int error;
  double result = 0;

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      build_key_string(key_len, key, size, automaton, offset, i, j);
      error = hashmap_get(map, key, (void**)(&value));
      if (error==MAP_MISSING) {
        result += - log2(1/(double)states);
      }
      else {
        result +=  - log2((value->number_array[automaton[i][j]] > PERT) ?
                          value->number_array[automaton[i][j]]: PERT);
      }
    }
  }
  return result / (double)(size * size);
}

double cross_entropy_between_arrays(size_t size, double* source, double* target)
{
  double result = 0;
  for (size_t i = 0; i < size; ++i) {
    result += -target[i] * log2((source[i] > PERT) ? source[i]: PERT);
  }
  return result;
}

int compute_cross_ent(any_t in_item, any_t in)
{
  int error;
  data_struct_t* source;
  data_struct_t* target = (data_struct_t*) in;
  cross_ent_struct_t* item = (cross_ent_struct_t*) in_item;

  /* Placeholder with uniform distribution over the states */
  double* base_array = calloc(item->states, sizeof(double));
  for (int i = 0; i < item->states; ++i) {
    base_array[i] = 1/(double)item->states;
  }

  double* ph_array;

  error = hashmap_get(item->source_map, target->key_string, (void**)(&source));
  /*  By default the prediction is same proability for all states */
  if (error==MAP_MISSING) {
    ph_array = base_array;
  }
  else {
    ph_array = source->number_array;
  }
  item->cross_entropy +=
    cross_entropy_between_arrays(item->states,
                                 ph_array,
                                 target->number_array);
  free(base_array);
  return MAP_OK;
}

double map_cross_entropy(map_t map_source, size_t size,
                         uint8_t automaton[size][size], int offset, int states)
{
  double cross_ent;
  int key_len = (2*offset + 1) * (2*offset + 1) + 1;
  char key[key_len];
  cross_ent_struct_t* item = calloc(1, sizeof(cross_ent_struct_t));
  item->states = states;
  item->source_map = map_source;

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      build_key_string(key_len, key, size, automaton, offset, i, j);
    }
  }

  cross_ent = item->cross_entropy;
  free(item);
  return cross_ent;
}

int free_data(any_t _, any_t in)
{
  (void)(_); /* Unused parameter */
  data_struct_t* data = (data_struct_t*) in;
  free(data->number_array);
  free(data);
  data = NULL;

  return MAP_OK;
}

void free_map(map_t map)
{
  hashmap_iterate(map, (PFany)&free_data, NULL);
  hashmap_free(map);
}

void add_results_to_file(map_t map_source, size_t size,
                         uint8_t automaton[size][size], int states,
                         int offset, FILE* file)
{
  entrop_state_ph_t* result = entropy_score(map_source, states);
  fprintf(file, "%f    %f    %"PRIu32"    ",
          predictive_score(map_source, states, size, automaton, offset),
          result->entropy, result->visited);
  free(result);
}


/**
 * Given a rule, create and simulate the corresponding automaton.
 */
void process_rule(uint64_t grule_size, uint8_t rule[grule_size],
                  char rule_buf[], int save_steps, int joint_complexity,
                  int totalistic, long steps, int states, int horizon,
                  size_t size, int grain, enum WriteStepMode save_flag)
{
  FILE* out_file = NULL;
  char* fname;

  FILE* mult_time_file;
  char* mult_time_fname;

  FILE* entrop_file = NULL;
  char* entrop_fname;

  FILE* nn_file = NULL;
  char* nn_fname;

  int last_compressed_size;
  int compressed_size;
  /* int agg_compressed_size; */
  int last_cell_count;
  int cell_count = 0;
  int size_sum;
  int offset_a = 2, offset_b = 1;
  float ratio = 0;
  float ratio2 = 0;

  ProcessF process_function = totalistic == 1 ?
    update_step_totalistic: update_step_general;

  asprintf(&fname, "data_2d_%i/out/out%s.dat", states, rule_buf);
  out_file = fopen(fname, "w+");

  asprintf(&mult_time_fname, "data_2d_%i/mult/mult%s.dat", states, rule_buf);
  mult_time_file = fopen(mult_time_fname, "w+");


  uint8_t A[size][size];
  uint8_t base[size][size];

  uint8_t (*automat5) [size] =
    (uint8_t (*) [size])malloc(size * size * sizeof(uint8_t));
  uint8_t (*automat50) [size] =
    (uint8_t (*) [size])malloc(size * size * sizeof(uint8_t));
  uint8_t (*automat300) [size] =
    (uint8_t (*) [size])malloc(size * size * sizeof(uint8_t));

  init_automat(size, A, states);
  memcpy(base, A, size * size * sizeof(uint8_t));

  char dbl_pholder[2 * ((size + 1) * size + 1)];
  char out_string[(size + 1) * size + 1];
  char out_string300[(size + 1) * size + 1];
  char out_string200[(size + 1) * size + 1];
  char out_string100[(size + 1) * size + 1];
  char out_string50[(size + 1) * size + 1];
  char out_string10[(size + 1) * size + 1];
  char out_string5[(size + 1) * size + 1];
  /* int dbl_cmp300, dbl_cmp200, dbl_cmp100, dbl_cmp50, dbl_cmp10, dbl_cmp5; */

  map_t map300 = hashmap_new();
  map_t map50 = hashmap_new();
  map_t map5 = hashmap_new();
  map_t map300b = hashmap_new();
  map_t map50b = hashmap_new();
  map_t map5b = hashmap_new();

  int flag = 0;

  for (int i = 0; i < steps; i++) {

    process_function(grule_size, size, base, rule, A, states, horizon);

    if (i % grain == 0) {
      last_compressed_size = compressed_size;
      last_cell_count = cell_count;

      print_bits(size, size, A, out_string);
      compressed_size = compress_memory_size(out_string, (size + 1) * size);
      cell_count = count_cells(size, A, states);

      /* Check if state has evolved from last time */
      if (compressed_size == last_compressed_size && flag == 1) {
        printf("\n");
        /* Cleanup before returning */
        break;
      } else if (compressed_size == last_compressed_size) {
        flag = 1;
      }

      if (joint_complexity == 1) {
        memcpy(&dbl_pholder[(size + 1) * size + 1],
               out_string, (size + 1) * size + 1);

        int dbl_comp_size =
          compress_memory_size(dbl_pholder, 2 * ((size + 1) * size + 1));

        memcpy(dbl_pholder, out_string, (size + 1) * size + 1);

        if (i > 0) {
          size_sum = last_compressed_size + compressed_size;
          ratio2 = (size_sum - dbl_comp_size)/(float)size_sum;
          ratio = ( (last_compressed_size / (float)last_cell_count) +
                    (compressed_size / (float)cell_count) ) /
            (dbl_comp_size / (float)(last_cell_count + cell_count));
        }
        fprintf(out_file, "%i    %i    %f    %f    "
                          "%i    %i    %i\n",
                i, compressed_size, ratio, ratio2, cell_count, last_cell_count,
                dbl_comp_size);
      } else {
        fprintf(out_file, "%i    %i\n", i, compressed_size);
      }

      printf("%i  ", compressed_size);
      fflush(stdout);

      if (save_steps == 1) {
        FILE* out_step_file;
        char* step_fname;
        if (save_flag == TMP_FILE) {
          asprintf(&step_fname, "rule_gif/tmp_%i.step", i);
        }
        else {
          asprintf(&step_fname, "step_2d_%i/out%s_%i.step",
                   states, rule_buf, i);
        }
        out_step_file = fopen(step_fname, "w+");
        free(step_fname);
        fprintf(out_step_file, "%s", out_string);
        fclose(out_step_file);
      }
    }

    if (i == steps - 301) {
      print_bits(size, size, A, out_string300);
      populate_map(map300, size, A, offset_a, states);
      populate_map(map300b, size, A, offset_b, states);

      memcpy(automat300, A, size * size * sizeof(uint8_t));
    }
    if (i == steps - 200) {
      print_bits(size, size, A, out_string200);
    }
    if (i == steps - 100) {
      print_bits(size, size, A, out_string100);
    }
    if (i == steps - 51) {
      print_bits(size, size, A, out_string50);
      populate_map(map50, size, A, offset_a, states);
      populate_map(map50b, size, A, offset_b, states);

      memcpy(automat50, A, size * size * sizeof(uint8_t));
    }
    if (i == steps - 10) {
      print_bits(size, size, A, out_string10);
    }
    if (i == steps - 5) {
      print_bits(size, size, A, out_string5);
      populate_map(map5, size, A, offset_a, states);
      populate_map(map5b, size, A, offset_b, states);

      memcpy(automat5, A, size * size * sizeof(uint8_t));
    }
    if (i == steps - 1) {
      /* print_bits(size, size, A, out_string); */
      /* memcpy(&dbl_pholder[(size + 1) * size + 1], */
      /*        out_string, (size + 1) * size + 1); */

      /* memcpy(dbl_pholder, out_string300, (size + 1) * size + 1); */
      /* dbl_cmp300 = compress_memory_size(dbl_pholder, */
      /*                                   2 * ((size + 1) * size + 1)); */
      /* memcpy(dbl_pholder, out_string200, (size + 1) * size + 1); */
      /* dbl_cmp200 = compress_memory_size(dbl_pholder, */
      /*                                   2 * ((size + 1) * size + 1)); */
      /* memcpy(dbl_pholder, out_string100, (size + 1) * size + 1); */
      /* dbl_cmp100 = compress_memory_size(dbl_pholder, */
      /*                                   2 * ((size + 1) * size + 1)); */
      /* memcpy(dbl_pholder, out_string50, (size + 1) * size + 1); */
      /* dbl_cmp50 = compress_memory_size(dbl_pholder, */
      /*                                   2 * ((size + 1) * size + 1)); */
      /* memcpy(dbl_pholder, out_string10, (size + 1) * size + 1); */
      /* dbl_cmp10 = compress_memory_size(dbl_pholder, */
      /*                                   2 * ((size + 1) * size + 1)); */
      /* memcpy(dbl_pholder, out_string5, (size + 1) * size + 1); */
      /* dbl_cmp5 = compress_memory_size(dbl_pholder, */
      /*                                  2 * ((size + 1) * size + 1)); */

      /* fprintf(mult_time_file, "%i    %i    %i    %i    " */
      /*                         "%i    %i    %i    %i    " */
      /*                         "%i    %i    %i    %i    " */
      /*                         "%i\n", */
      /*         compress_memory_size(out_string, (size + 1) * size), */
      /*         dbl_cmp5, */
      /*         compress_memory_size(out_string5, (size + 1) * size), */
      /*         dbl_cmp10, */
      /*         compress_memory_size(out_string10, (size + 1) * size), */
      /*         dbl_cmp50, */
      /*         compress_memory_size(out_string50, (size + 1) * size), */
      /*         dbl_cmp100, */
      /*         compress_memory_size(out_string100, (size + 1) * size), */
      /*         dbl_cmp200, */
      /*         compress_memory_size(out_string200, (size + 1) * size), */
      /*         dbl_cmp300, */
      /*         compress_memory_size(out_string300, (size + 1) * size)); */

      asprintf(&entrop_fname, "data_2d_%i/ent/ent%s.dat", states, rule_buf);
      entrop_file = fopen(entrop_fname, "w+");

      add_results_to_file(map300, size, A, states, offset_a, entrop_file);
      add_results_to_file(map50, size, A, states, offset_a, entrop_file);
      add_results_to_file(map5, size, A, states, offset_a, entrop_file);
      fprintf(entrop_file, "\n");


      add_results_to_file(map300b, size, A, states, offset_b, entrop_file);
      add_results_to_file(map50b, size, A, states, offset_b, entrop_file);
      add_results_to_file(map5b, size, A, states, offset_b, entrop_file);
      fprintf(entrop_file, "\n");


      asprintf(&nn_fname, "data_2d_%i/nn/nn%s.dat", states, rule_buf);
      nn_file = fopen(nn_fname, "w+");
      train_nn_on_automaton(size, states, automat300, A, 7, nn_file);
      train_nn_on_automaton(size, states, automat50, A, 7, nn_file);
      train_nn_on_automaton(size, states, automat5, A, 7, nn_file);
    }
  }
  printf("\n");

  /* Cleanup before finishing */

  free_map(map5);
  free_map(map300);
  free_map(map50);
  free_map(map5b);
  free_map(map300b);
  free_map(map50b);

  free(automat5);
  free(automat50);
  free(automat300);

  free(fname);
  free(mult_time_fname);

  if (entrop_file) {
    free(entrop_fname);
    fclose(entrop_file);
  }
  if (nn_file) {
    free(nn_fname);
    fclose(nn_file);
  }
  fclose(mult_time_file);
  fclose(out_file);
}

