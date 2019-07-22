#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "2d_automaton.h"
#include "automaton/rule.h"
#include "nn/nn.h"
#include "utils/compress.h"
#include "utils/utils.h"
#include "utils/hashmap.h"

#define KEY_MAX_LENGTH (256)
#define PERT 1E-15

typedef void (*ProcessF)(uint64_t, size_t size, uint8_t*,
                         uint8_t[] , uint8_t* , int, int);

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
void print_bits(size_t dim1, size_t dim2, uint8_t* a, char* buf)
{
  char out_string = '0';

  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      out_string = '0' + a[i * dim2 + j];
      buf[i * (dim1 + 1) + j] = out_string;
    }
    buf[i * (dim1 + 1) + dim2] = '\n';
  }
  buf[(dim1 + 1) * dim2] = '\0';
}

/**
 * @brief Automaton initializer function.
 *
 */
void init_automat(size_t size, uint8_t* a, int states)
{
  assert(size > 20);
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      /* Initialize center 20x20 square to random */
      if ((i >= size/2 - 10) && (i < size/2 + 10) &&
          (j >= size/2 - 10) && (j < size/2 + 10)) {
        a[(i % size) * size + (j % size)] = (uint8_t)(rand() % states);
      }
      /* Rest is 0 */
      else {
        a[i * size + j] = (uint8_t)0;
      }
    }
  }
}

void update_step_general(uint64_t grule_size, size_t size,
                         uint8_t* A,
                         uint8_t rule[grule_size],
                         uint8_t* base, int states,
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
          current_value = base[((i + k + size) % size) * size
                               + ((j + l + size) % size)];
          position += current_value * ipow(states, increment);
          ++increment;
        }
      }

      A[i * size + j] = rule[position];
    }
  }
  memcpy(base, A, size * size * sizeof(uint8_t));
}

void update_step_totalistic(uint64_t rule_size, size_t size,
                            uint8_t* base,
                            uint8_t rule[rule_size],
                            uint8_t* A, int states,
                            int horizon)
{
  int count;
  memcpy(A, base, size * size * sizeof(uint8_t));

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      count = 0;
      for (int k = - horizon; k <= horizon; k++) {
        for (int l = - horizon; l <= horizon; l++) {
          if (k != 0 | l != 0) {
            count += base[((i + k + size) % size) * size
                          + ((j + l + size) % size)];
          }
        }
      }
      A[i * size + j] = rule[states * count + base[i * size + j]];
    }
  }
  memcpy(base, A, size * size * sizeof(uint8_t));
}

int count_cells(size_t size, uint8_t* A, int states)
{
  int* counts = (int*)calloc(states, sizeof(int));
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      counts[A[i* size + j]] += 1;
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
    if (in->number_array[i] == 0.0) {
      in->number_array[i] += 1E-3;
    }
    sum += in->number_array[i];
  }
  for (int i = 0; i < st; i++) {
    in->number_array[i] /= sum;
  }
  return MAP_OK;
}

void build_key_string(int key_len, char key[key_len],
                      size_t size, uint8_t* automaton,
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
          '0' + automaton[((i + a + size) % size) * size
                          + ((j + b + size) % size)];
      }
    }
  }
  key[key_len - 1] = '\0';
}

entrop_state_ph_t* populate_map(map_t* map, size_t size,
                                uint8_t* automaton,
                                int offset, int states)
{
  data_struct_t* value;
  entrop_state_ph_t* out_data =
    (entrop_state_ph_t *) calloc(1, sizeof(entrop_state_ph_t));

  int key_len = (2*offset + 1) * (2*offset + 1) + 1;
  char key[key_len];
  int error;
  out_data->entropy = 0.0;

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
      value->number_array[automaton[i * size + j]] += 1;
    }
  }
  hashmap_iterate(map, (PFany)&normalize_probs, (any_t)&states);

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      build_key_string(key_len, key, size, automaton, offset, i, j);
      error = hashmap_get(map, key, (void**)(&value));
      if (error==MAP_MISSING) {
        out_data->entropy += -log( 1 / (double) states );
      }
      else {
        out_data->entropy +=
          -log(value->number_array[automaton[i * size + j]]);
      }
    }
  }
  out_data->entropy /= (double) (size * size);
  return out_data;
}

int reduce_entropy(any_t state_entropy, any_t in)
{
  entrop_state_ph_t* state_ent = (entrop_state_ph_t*) state_entropy;
  data_struct_t* data = (data_struct_t*) in;

  for (int i = 0; i < state_ent->states; i ++) {
    if (data->number_array[i] > 0) {
      state_ent->entropy += -log(data->number_array[i]);
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
                        uint8_t* automaton, int offset)
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
        result += - log(1/(double)states);
      }
      else {
        result +=  - log((value->number_array[automaton[i * size + j]] > PERT) ?
                          value->number_array[automaton[i * size + j]]: PERT);
      }
    }
  }
  return result / (double)(size * size);
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
                         uint8_t* automaton, int states,
                         int offset, FILE* file, entrop_state_ph_t* result)
{
  fprintf(file, "%f    %f    %"PRIu32"    ",
          predictive_score(map_source, states, size, automaton, offset),
          result->entropy, result->visited);
  free(result);
}

void add_nn_results_to_file(FILE* file, network_result_t* res)
{
  fprintf(file, "%f    %f\n", res->train_error, res->test_error);
}

void process_rule(uint64_t grule_size, uint8_t rule[grule_size],
                  char rule_buf[],
                  int totalistic, long steps,
                  struct Options2D* opts, results_nn_t* results)
{
  int states = opts->states;
  size_t size = opts->size;

  FILE* out_file = NULL;
  char* fname;

  FILE* mult_time_file;
  char* mult_time_fname;

  FILE* entrop_file = NULL;
  char* entrop_fname;

  FILE* nn_file = NULL;
  char* nn_fname = NULL;

  int last_compressed_size;
  int compressed_size;
  int last_cell_count;
  int cell_count = 0;
  int size_sum;
  float ratio = 0;
  float ratio2 = 0;

  ProcessF process_function = totalistic == 1 ?
    update_step_totalistic: update_step_general;

  asprintf(&fname, "data_2d_%i/out/out%s.dat", states, rule_buf);
  out_file = fopen(fname, "w+");

  asprintf(&mult_time_fname, "data_2d_%i/mult/mult%s.dat", states, rule_buf);
  mult_time_file = fopen(mult_time_fname, "w+");


  uint8_t* A = (uint8_t*) calloc(size * size, sizeof(uint8_t));
  uint8_t* base = (uint8_t*) calloc(size * size, sizeof(uint8_t));

  uint8_t* automat5 = (uint8_t*) calloc(size * size, sizeof(uint8_t));
  uint8_t* automat50 = (uint8_t*) calloc(size * size, sizeof(uint8_t));
  uint8_t* automat300 = (uint8_t*) calloc(size * size, sizeof(uint8_t));

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
  entrop_state_ph_t *res300, *res50, *res5, *res300b, *res50b, *res5b;

  int flag = 0;

  for (int i = 0; i < steps; i++) {

    process_function(grule_size, size, base, rule, A, states, opts->horizon);

    if (i % opts->grain == 0) {
      last_compressed_size = compressed_size;
      last_cell_count = cell_count;

      print_bits(size, size, A, out_string);
      compressed_size = compress_memory_size(out_string, (size + 1) * size);
      cell_count = count_cells(size, A, states);

      /* Check if state has evolved from last time (stop mechanism) */
      if (compressed_size == last_compressed_size
          && flag == 1 && opts->early == EARLY) {
        printf("\n");
        break;
      }
      else if (compressed_size == last_compressed_size) {
        flag = 1;
      }

      if (opts->joint_complexity == 1) {
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

    }

    /* Save steps every grain_write */
    if (opts->grain_write > 0
        && i % opts->grain_write == 0
        && opts->save_steps == 1) {

      FILE* out_step_file;
      char* step_fname;

      print_bits(size, size, A, out_string);
      if (opts->save_flag == TMP_FILE) {
        asprintf(&step_fname, "%s/tmp_%i.step", opts->out_step_dir, i);
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

    int offset_a = 2, offset_b = 1;

    if (i == steps - 301) {
      print_bits(size, size, A, out_string300);
      res300 = populate_map(map300, size, A, offset_a, states);
      res300b = populate_map(map300b, size, A, offset_b, states);

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
      res50 = populate_map(map50, size, A, offset_a, states);
      res50b = populate_map(map50b, size, A, offset_b, states);

      memcpy(automat50, A, size * size * sizeof(uint8_t));
    }
    if (i == steps - 10) {
      print_bits(size, size, A, out_string10);
    }
    if (i == steps - 5) {
      print_bits(size, size, A, out_string5);
      res5 = populate_map(map5, size, A, offset_a, states);
      res5b = populate_map(map5b, size, A, offset_b, states);

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

      add_results_to_file(map300, size, A, states, offset_a,
                          entrop_file, res300);
      add_results_to_file(map50, size, A, states, offset_a,
                          entrop_file, res50);
      add_results_to_file(map5, size, A, states, offset_a,
                          entrop_file, res5);
      fprintf(entrop_file, "\n");


      add_results_to_file(map300b, size, A, states, offset_b,
                          entrop_file, res300b);
      add_results_to_file(map50b, size, A, states, offset_b,
                          entrop_file, res50b);
      add_results_to_file(map5b, size, A, states, offset_b,
                          entrop_file, res5b);
      fprintf(entrop_file, "\n");


      asprintf(&nn_fname, "data_2d_%i/nn/nn%s.dat", states, rule_buf);
      nn_file = fopen(nn_fname, "w+");

      network_result_t res = {1.};
      network_opts_t n_opts = {10, 30, 3};
      train_nn_on_automaton(size, states, automat300, A, &n_opts, &res);
      add_nn_results_to_file(nn_file, &res);
      results->nn_tr_300 = res.train_error;
      results->nn_te_300 = res.test_error;

      /* train_nn_on_automaton(size, states, automat50, A, &n_opts, &res); */
      /* add_nn_results_to_file(nn_file, &res); */
      /* results->nn_tr_50 = res.train_error; */
      /* results->nn_te_50 = res.test_error; */

      /* train_nn_on_automaton(size, states, automat5, A, &n_opts, &res); */
      /* add_nn_results_to_file(nn_file, &res); */
      /* results->nn_tr_5 = res.train_error; */
      /* results->nn_te_5 = res.test_error; */

    /*   n_opts.offset = 1; */
    /*   train_nn_on_automaton(size, states, automat300, A, &n_opts, &res); */
    /*   add_nn_results_to_file(nn_file, &res); */
    /*   train_nn_on_automaton(size, states, automat50, A, &n_opts, &res); */
    /*   add_nn_results_to_file(nn_file, &res); */
    /*   train_nn_on_automaton(size, states, automat5, A, &n_opts, &res); */
    /*   add_nn_results_to_file(nn_file, &res); */

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

