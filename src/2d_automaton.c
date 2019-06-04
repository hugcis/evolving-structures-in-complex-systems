#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "2d_automaton.h"
#include "nn.h"
#include "compress.h"
#include "utils.h"
#include "hashmap.h"

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
                         uint8_t base[size][size],
                         uint8_t rule[grule_size],
                         uint8_t** A, int states,
                         int horizon)
{
  uint64_t position;
  uint8_t current_value;
  int increment;

  for (size_t i = 0; i < size; i++) {
    memcpy(A[i], base[i], size * sizeof(uint8_t));

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
  for (size_t i = 0; i < size; i++) {
    memcpy(base[i], A[i], size * sizeof(uint8_t));
  }
}

void update_step_totalistic(uint64_t rule_size, size_t size,
                            uint8_t base[size][size],
                            uint8_t rule[rule_size],
                            uint8_t** A, int states,
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
  assert(A && size);
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
        result += 1/(double)states;
      }
      else {
        result += value->number_array[automaton[i][j]];
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
  double* base_array = calloc(item->states, sizeof(double));
  double* ph_array;

  error = hashmap_get(item->source_map, target->key_string, (void**)(&source));
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
                  size_t size, int grain)
{
  FILE* out_file = NULL;
  char* fname;

  FILE* mult_time_file;
  char* mult_time_fname;

  FILE* entrop_file = NULL;
  char* entrop_fname;

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
  int dbl_cmp300, dbl_cmp200, dbl_cmp100, dbl_cmp50, dbl_cmp10, dbl_cmp5;

  map_t map300 = hashmap_new();
  map_t map50 = hashmap_new();
  map_t map5 = hashmap_new();
  map_t map300b = hashmap_new();
  map_t map50b = hashmap_new();
  map_t map5b = hashmap_new();

  map_t map_a = hashmap_new();
  map_t map_b = hashmap_new();

  int flag = 0;

  for (int i = 0; i < steps; i++) {
    process_function(grule_size, size, A, rule, placeholder, states, horizon);

    if (i % grain == 0) {
      last_compressed_size = compressed_size;
      last_cell_count = cell_count;

      /* print_aggregated_bits(size, size, A, out_string, 2); */
      /* agg_compressed_size = compress_memory_size(out_string,
         (size/2 + 1) * size); */

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
        fprintf(out_file, "%i    %i    %f    %f"
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
        asprintf(&step_fname, "step_2d_%i/out%s_%i.step",
                 states, rule_buf, i);
        out_step_file = fopen(step_fname, "w+");
        free(step_fname);
        fprintf(out_step_file, "%s", out_string);
        fclose(out_step_file);
      }
    }

    if (i == steps - 300) {
      print_bits(size, size, A, out_string300);
      populate_map(map300, size, A, offset_a, states);
      populate_map(map300b, size, A, offset_b, states);
    }
    if (i == steps - 200) {
      print_bits(size, size, A, out_string200);
    }
    if (i == steps - 100) {
      print_bits(size, size, A, out_string100);
    }
    if (i == steps - 50) {
      print_bits(size, size, A, out_string50);
      populate_map(map50, size, A, offset_a, states);
      populate_map(map50b, size, A, offset_b, states);
    }
    if (i == steps - 10) {
      print_bits(size, size, A, out_string10);
    }
    if (i == steps - 5) {
      print_bits(size, size, A, out_string5);
      populate_map(map5, size, A, offset_a, states);
      populate_map(map5b, size, A, offset_b, states);
    }
    if (i == steps - 1) {
      print_bits(size, size, A, out_string);
      memcpy(&dbl_pholder[(size + 1) * size + 1],
             out_string, (size + 1) * size + 1);

      memcpy(dbl_pholder, out_string300, (size + 1) * size + 1);
      dbl_cmp300 = compress_memory_size(dbl_pholder,
                                        2 * ((size + 1) * size + 1));
      memcpy(dbl_pholder, out_string200, (size + 1) * size + 1);
      dbl_cmp200 = compress_memory_size(dbl_pholder,
                                        2 * ((size + 1) * size + 1));
      memcpy(dbl_pholder, out_string100, (size + 1) * size + 1);
      dbl_cmp100 = compress_memory_size(dbl_pholder,
                                        2 * ((size + 1) * size + 1));
      memcpy(dbl_pholder, out_string50, (size + 1) * size + 1);
      dbl_cmp50 = compress_memory_size(dbl_pholder,
                                        2 * ((size + 1) * size + 1));
      memcpy(dbl_pholder, out_string10, (size + 1) * size + 1);
      dbl_cmp10 = compress_memory_size(dbl_pholder,
                                        2 * ((size + 1) * size + 1));
      memcpy(dbl_pholder, out_string5, (size + 1) * size + 1);
      dbl_cmp5 = compress_memory_size(dbl_pholder,
                                       2 * ((size + 1) * size + 1));

      fprintf(mult_time_file, "%i    %i    %i    %i    "
                              "%i    %i    %i    %i    "
                              "%i    %i    %i    %i    "
                              "%i\n",
              compress_memory_size(out_string, (size + 1) * size),
              dbl_cmp5,
              compress_memory_size(out_string5, (size + 1) * size),
              dbl_cmp10,
              compress_memory_size(out_string10, (size + 1) * size),
              dbl_cmp50,
              compress_memory_size(out_string50, (size + 1) * size),
              dbl_cmp100,
              compress_memory_size(out_string100, (size + 1) * size),
              dbl_cmp200,
              compress_memory_size(out_string200, (size + 1) * size),
              dbl_cmp300,
              compress_memory_size(out_string300, (size + 1) * size));

      asprintf(&entrop_fname, "data_2d_%i/ent/ent%s.dat", states, rule_buf);
      entrop_file = fopen(entrop_fname, "w+");

      populate_map(map_a, size, A, offset_a, states);
      populate_map(map_b, size, A, offset_b, states);

      add_results_to_file(map300, map_a, states, entrop_file);
      add_results_to_file(map50, map_a, states, entrop_file);
      add_results_to_file(map5, map_a, states, entrop_file);
      fprintf(entrop_file, "\n");


      add_results_to_file(map300b, map_b, states, entrop_file);
      add_results_to_file(map50b, map_b, states, entrop_file);
      add_results_to_file(map5b, map_b, states, entrop_file);
      fprintf(entrop_file, "\n");
    }
  }
  printf("\n");

  /* Cleanup before finishing */
  for (size_t i = 0; i < size; i++) {
    free(placeholder[i]);
  }
  free(placeholder);

  free_map(map5);
  free_map(map300);
  free_map(map50);
  free_map(map5b);
  free_map(map300b);
  free_map(map50b);

  free_map(map_a);
  free_map(map_b);

  free(fname);
  free(mult_time_fname);

  if (entrop_file) {
    free(entrop_fname);
    fclose(entrop_file);
  }
  fclose(mult_time_file);
  fclose(out_file);
}

void generate_totalistic_rule(uint64_t rule_size, uint8_t rule_array[rule_size],
                              char rule_buf[rule_size + 1], int states)
{
  unsigned long rule_number = 0UL;

  for (uint64_t s = 0 ; s < rule_size; ++s) {
    rule_array[s] = rand() % states;
    rule_number += rule_array[s] * ipow(states, s);
    if (states >= 3) {
      rule_buf[s] = '0' + rule_array[s];
    }
  }

  if (states == 2) {
    sprintf(rule_buf, "%lu", rule_number);
  }
}

/**
 * Symmetrize the rule by setting all the states and their symmetries to having
 * the same output.
 */
void symmetrize_rule(uint64_t grule_size,
                     uint8_t rule_array[grule_size],
                     int states, int horizon)
{
  int side = 2 * horizon + 1;
  int neigh_size = side * side - 1;

  uint32_t position_180;
  uint32_t position_90;
  uint32_t position_270;
  uint32_t position_vflip;
  uint32_t position_hflip;
  uint32_t position_dflip;
  uint32_t position_aflip;

  /* grule_size can be very big, this array is better on the heap */
  /* Keep track of already seen positions with book-keeping */
  uint8_t* book_keep = calloc(grule_size, sizeof(uint8_t));

  int pos;

  for (uint64_t i = 0; i < grule_size; ++i) {
    /* Skip already seen positions when looping through the rule */
    if (book_keep[i] == 1) {
      continue;
    }

    position_180 = 0;
    position_90 = 0;
    position_270 = 0;
    position_vflip = 0;
    position_hflip = 0;
    position_dflip = 0;
    position_aflip = 0;

    /* Create the representation of the symmetrized position by swapping the
       states in its number representation. */
    for (int p = 0; p < neigh_size + 1; ++p) {
      /* 180° rotation */
      position_180 += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, neigh_size  - p)) % states);

      /* 90° rotation */
      pos = (neigh_size - side + 1 - (side * (p%side)) + p/side);
      position_90 += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, pos)) % states);

      /* 270° rotation */
      position_270 += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, neigh_size  - pos)) % states);

      /* Vertical flip */
      pos = (side * (p / side)) + (side - 1 - (p % 3));
      position_vflip += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, pos)) % states);

      /* Horizontal flip */
      if (p/side < side/2) {
        pos = (p - side + neigh_size + 1)%(neigh_size + 1);
      } else if (p/side > side/2) {
        pos = (p + side)%(neigh_size + 1);
      } else {
        pos = p;
      }
      position_hflip += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, pos)) % states);

      /* Diagonal flip */
      pos = neigh_size - ((p % side) * side  + (p / side));
      position_dflip += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, pos)) % states);

      /* Antidiagonal flip */
      pos = ((p % side) * side  + (p / side));
      position_aflip += (uint32_t)ipow(states, p) *
        ((i / (uint32_t)ipow(states, pos)) % states);

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
  free(book_keep);
}

/**
 * Build a rule from the provided command-line arguments.
 */
void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char* rule_arg, int states)
{
  /* Rule is given in base-(states - 1) format */
  for (uint64_t s = 0 ; s < grule_size; ++s) {
    rule_array[s] = rule_arg[s] - '0';
    rule_buf[s] = rule_arg[s];
    assert(rule_array[s] < states);
  }
  rule_buf[grule_size] = '\0';
}

/**
 * Generate a random general rule. The rule is written to rule_array and a
 * string  representation is also saved in rule_buf.
 */
void generate_general_rule(uint64_t grule_size,
                           uint8_t rule_array[grule_size],
                           char rule_buf[grule_size + 1],
                           int states, int horizon)
{
  int inc;
  /* The following implements a Uniform-lambda sampling strategy */

  /* Choose lambda parameter at random as well as the proportion of */
  /* transitions to other states*/
  float lambda[states - 1];
  for (int i = 0; i < states - 1; ++i) {
    lambda[i] = ((double)rand() / (double)((unsigned)RAND_MAX + 1));
  }

  for (uint64_t v = 0; v < grule_size; v++) {
    inc = 0;
    /* Assign the rule to the first state that passes the test */
    while ((double)rand() / (double)((unsigned)RAND_MAX + 1) > lambda[inc]
           && inc < states - 1) {
      inc++;
    }
    rule_array[v] = (uint8_t)inc;
  }

  symmetrize_rule(grule_size, rule_array, states, horizon);

  /* Populate buffer with newly created rule */
  for (uint64_t v = 0; v < grule_size; v++) {
    sprintf(&rule_buf[v], "%"PRIu8, rule_array[v]);
  }
}
