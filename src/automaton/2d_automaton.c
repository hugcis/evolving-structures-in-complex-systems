#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <ctype.h>
#include "automaton/2d_automaton.h"
#include "automaton/rule.h"
#include "nn/nn.h"
#include "utils/compress.h"
#include "utils/utils.h"
#include "utils/hashmap.h"

#if PROFILE
#define PROF(x) {\
    clock_t t1=clock();\
    x;\
    clock_t t2=clock();\
    t += t2 - t1;\
    printf("\rLoop per second %g ", ((double) i / ((double) t * 1E-6)));\
  }
#else
#define PROF(x) x;
#endif

#define KEY_MAX_LENGTH (256)
#define PERT 1E-15
#define WINDOW 501
#define W_STEP 13

typedef struct masking_element
{
  uint8_t state;
  size_t pos_x;
  size_t pos_y;
} masking_element_t;

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
void init_automat(size_t size, uint8_t* a, int states, long init_type)
{
  assert((long) size > init_type);
  size_t border = (init_type > 0) ? init_type/2: size/2;

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      /* Very small automata: initialize at random */
      if (size <= 20) {
        a[(i % size) * size + (j % size)] = (uint8_t)(rand() % states);
      }
      else {
        /* Initialize center (border x border) square to random */
        if ((i >= size/2 - border) && (i < size/2 + border) &&
            (j >= size/2 - border) && (j < size/2 + border)) {
          a[(i % size) * size + (j % size)] = (uint8_t)(rand() % states);
        }
        /* Rest is 0 */
        else {
          a[i * size + j] = (uint8_t)0;
        }
      }
    }
  }
}

/**
 * @brief A function that will read a pattern file and initialize the automaton
 * with the specified pattern.
 */
int parse_pattern(size_t size, uint8_t* a, FILE* pattern_file)
{
  char ch; /* Character placeholder for the file */
  int pattern_start = 0; /* Flag indicating wether we are reading the header
                            or the pattern */
  size_t x = size/2, y = size/2; /* Pattern's upper left corner is the center
                                    of the automaton */

  while ((ch = fgetc(pattern_file)) != EOF && pattern_start != -1)  {
    switch (ch) {
    case '#':
        if (pattern_start == 0) {
          pattern_start = 1; /* Pattern is starting */
        }
        else if (pattern_start == 1) {
          pattern_start = 2;  /* Pattern ends */
        }
        break;
    case '\n':
      if (pattern_start == 1) {
        /* Line breaks indicate new line in the pattern */
        x++;
        y = size/2;
      }
      break;
    case 'B': /* Beginning of BG flag to fill background with a state */
      /* We need to read the next 3 characters since the format is BG=X */
      if ((ch = fgetc(pattern_file)) != EOF
          && ch == 'G'
          && (ch = fgetc(pattern_file)) != EOF
          && ch == '='
          && (ch = fgetc(pattern_file)) != EOF
          && isdigit(ch) != 0) {
        for (size_t i = 0; i < size; ++i) {
          for (size_t j = 0; j < size; ++j) {
            a[i * size + j] = (uint8_t)(ch - '0');
          }
        }
      }
      break;
    default:
      /* If we are reading the pattern and a digit comes up, add it
         to the automaton and increment coordinates */
      if (pattern_start == 1 && isdigit(ch) != 0) {
        a[x * size + y] = (uint8_t)(ch - '0');
        y++;
      }
      break;
    }
  }

  fclose(pattern_file);
  return pattern_start;
}

/**
 * Function that updates the automaton state from the buffer last_autom to the
 * buffer autom. It uses a rule stored as a lookup table.
 */
void update_step_general(size_t size,
                         uint8_t* autom,
                         uint8_t* rule,
                         uint8_t* last_autom,
                         int horizon,
                         uint32_t* pows)
{
  uint32_t position;
  uint8_t current_value;
  int increment;

  /* Take care of cells not on boundary first */
  for (size_t i = horizon; i < size - horizon; i++) {
    for (size_t j = horizon; j < size - horizon; j++) {
      position = 0;
      increment = 0;
      for (int k = - horizon; k <= horizon; k++) {
        for (int l = - horizon; l <= horizon; l++) {
          current_value = last_autom[(i + k) * size
                                     + (j + l)];
          position += current_value * pows[increment];
          ++increment;
        }
      }
      autom[i * size + j] = rule[position];
    }
  }

  /* Deal with boundary conditions separately */
  /* TODO: These boundary conditions updates only work with horizon=1 */
  /* First and last line */
  for (size_t j = 0; j < size; ++j) {
    position = 0;
    increment = 0;
    for (int k = - horizon; k <= horizon; k++) {
      for (int l = - horizon; l <= horizon; l++) {
        current_value = last_autom[((k + size) % size) * size
                                   + ((j + l + size) % size)];
        position += current_value * pows[increment];
        ++increment;
      }
    }
    autom[j] = rule[position];

    position = 0;
    increment = 0;
    for (int k = - horizon; k <= horizon; k++) {
      for (int l = - horizon; l <= horizon; l++) {
        current_value = last_autom[((size - 1 + k + size) % size) * size
                                   + ((j + l + size) % size)];
        position += current_value * pows[increment];
        ++increment;
      }
    }
    autom[(size - 1) * size + j] = rule[position];
  }

  /* First and last column */
  for (size_t i = 0; i < size; ++i) {
    position = 0;
    increment = 0;
    for (int k = - horizon; k <= horizon; k++) {
      for (int l = - horizon; l <= horizon; l++) {
        current_value = last_autom[((i + k + size) % size) * size
                                   + ((l + size) % size)];
        position += current_value * pows[increment];
        ++increment;
      }
    }
    autom[i * size] = rule[position];

    position = 0;
    increment = 0;
    for (int k = - horizon; k <= horizon; k++) {
      for (int l = - horizon; l <= horizon; l++) {
        current_value = last_autom[((i + k + size) % size) * size
                                   + ((size - 1 + l + size) % size)];
        position += current_value * pows[increment];
        ++increment;
      }
    }
    autom[i * size + size - 1] = rule[position];
  }
}


/**
 * Count the cells in each state and return the minority state cell count
 */
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

void add_entropy_results_to_file(map_t map_source, size_t size,
                         uint8_t* automaton, int states,
                         int offset, FILE* file, entrop_state_ph_t* result)
{
  if (result) {
    fprintf(file, "%f    %f    %"PRIu32"    ",
            predictive_score(map_source, states, size, automaton, offset),
            result->entropy, result->visited);
    free(result);
  }
}

/**
 * Add results computed with the neural network to a file. Some parameters are
 * also written to have tracability of the results.
 */
void add_nn_results_to_file(FILE* file, network_opts_t* opts,
                            network_result_t* res, int timesteps)
{
  fprintf(file, "%i %i %i %i:%f %f\n", opts->num_hid, opts->max_epoch,
          opts->offset, timesteps, res->train_error, res->test_error);
}

void add_fisher_results_to_file(FILE* file, network_opts_t* opts,
                                network_result_t* res, int timesteps)
{
  fprintf(file, "%i %i %i %i:%f\n", opts->num_hid, opts->max_epoch,
          opts->offset, timesteps, res->fisher_info);
}

void mask_autom(int pert, size_t size, masking_element_t* mask,
                uint8_t* automaton)
{
  for (int i = 0; i < pert; ++i) {
    automaton[mask[i].pos_x * size + mask[i].pos_y] = mask[i].state;
  }
}

void make_mask(int pert, masking_element_t* mask, size_t size, int states)
{
  for (int i = 0; i < pert; ++i) {
    mask[i].pos_x = size
      * ((double)rand() / (double)((unsigned)RAND_MAX + 1));
    mask[i].pos_y = size
      * ((double)rand() / (double)((unsigned)RAND_MAX + 1));
    mask[i].state = rand() % states;
  }
}

/**
 * Function to add random noise to an automaton with probability `rate`. Noise
 * consists in switching a cell to a random different state.
 */
void random_noise(size_t size, uint8_t* automaton, int states, double rate)
{
  size_t rand_x, rand_y;
  int pert = (int) (rate * (double) size * (double) size);
  for (int i = 0; i < pert; ++i) {
    rand_x = size * ((double)rand() / (double)((unsigned)RAND_MAX + 1));
    rand_y = size * ((double)rand() / (double)((unsigned)RAND_MAX + 1));

    automaton[(rand_x % size) * size + (rand_y % size)]
      += (1 + rand() % (states - 1));
    automaton[(rand_x % size) * size + (rand_y % size)]
      %= states;
  }
}

/**
 * Main entrypoint that handles all the automaton processing.
 */
void process_rule(uint64_t grule_size, uint8_t rule[grule_size],
                  char rule_buf[],
                  long steps,
                  struct Options2D* opts, results_nn_t* results)
{
  int states = opts->states;
  size_t size = opts->size;

  FILE* out_file = NULL;
  char* fname = NULL;

  FILE* mult_time_file = NULL;
  char* mult_time_fname = NULL;

  FILE* entrop_file = NULL;
  char* entrop_fname = NULL;

  FILE* nn_file = NULL;
  char* nn_fname = NULL;

  FILE* fisher_file = NULL;
  char* fisher_fname = NULL;

  int last_compressed_size;
  int compressed_size;
  int last_cell_count;
  int cell_count = 0;
  /* int size_sum; */
  /* float ratio = 0; */
  /* float ratio2 = 0; */

  ProcessF process_function = update_step_general;


  uint8_t* (*frame1) = malloc(sizeof(uint8_t* (*)));
  *frame1 = (uint8_t*) malloc(size * size * sizeof(uint8_t));
  uint8_t* (*frame2) = malloc(sizeof(uint8_t* (*)));
  *frame2 = (uint8_t*) malloc(size * size * sizeof(uint8_t));
  uint8_t* (*temp_frame);

  int pert = (int) (opts->noise_rate * (double) size * (double) size);

  masking_element_t* mask = NULL;
  if (opts->mask == MASK) {
      mask = (masking_element_t*) calloc(pert, sizeof(masking_element_t));
  }


  if (opts->init_pattern_file == NULL) {
    init_automat(size, *frame1, states, opts->init_type);
  }
  else if (parse_pattern(size, *frame1, opts->init_pattern_file) <= 0){
    fprintf(stderr, "Error parsing initialization pattern\n");
    exit(EXIT_FAILURE);
  }
  memcpy(*frame2, *frame1, size * size * sizeof(uint8_t));
  uint8_t* automat5 = NULL;
  uint8_t* automat50 = NULL;
  uint8_t* automat300 = NULL;

  uint8_t** test_automata =
    (uint8_t**) malloc(sizeof(uint8_t*) * (WINDOW / W_STEP));

  if (opts->output_data != NO_OUTPUT) {
    /* asprintf(&fname, "data_2d_%i/out/out%s.dat", states, rule_buf); */
    /* out_file = fopen(fname, "w+"); */

    automat5 = (uint8_t*) calloc(size * size, sizeof(uint8_t));
    automat50 = (uint8_t*) calloc(size * size, sizeof(uint8_t));
    automat300 = (uint8_t*) calloc(size * size, sizeof(uint8_t));

    for (int i = 0; i < (WINDOW / W_STEP); ++i) {
      test_automata[i] = (uint8_t*) calloc(size * size, sizeof(uint8_t));
    }
  }

  /* char* dbl_pholder = (char*) */
    /* malloc(sizeof(char) * (steps  * ((size + 1) * size + 1))); */
  char out_string[(size + 1) * size + 1];
  char out_string300[(size + 1) * size + 1];
  char out_string50[(size + 1) * size + 1];
  char out_string5[(size + 1) * size + 1];

  map_t map300 = hashmap_new();
  map_t map50 = hashmap_new();
  map_t map5 = hashmap_new();
  map_t map300b = hashmap_new();
  map_t map50b = hashmap_new();
  map_t map5b = hashmap_new();
  entrop_state_ph_t *res300 = NULL, *res50 = NULL, *res5 = NULL,
    *res300b = NULL, *res50b = NULL, *res5b = NULL;

  int flag = 0;
  int neigs = (2 * opts->horizon + 1) * (2 * opts->horizon + 1);
  #if PROFILE
  clock_t t = 0;
  #endif

  uint32_t pows[neigs];
  for (int i = 0; i < neigs; ++i) {
    pows[i] = ipow(states, i);
  }

  for (int i = 0; i < steps; ++i) {

    /* if ((i + 1) % opts->noise_step == 0) { */
    /*   random_noise(size, frame1, states, opts->noise_rate); */
    /* } */
    if (i % 20 == 0 && opts->mask == MASK) {
      make_mask(pert, mask, size, opts->states);
    }

    if (opts->mask == MASK) {
      mask_autom(pert, size, mask, *frame1);
    }
    /* Make update from frame1 to frame2 */

    /* Macro to profile if flag is on */
    PROF(
      process_function(size, *frame2, rule, *frame1, opts->horizon, pows);
      /* Swap pointers */
      temp_frame = frame1;
      frame1 = frame2;
      frame2 = temp_frame;
    )

    if (opts->mask == MASK) {
      mask_autom(pert, size, mask, *frame1);
    }

    /* Save steps every grain_write */
    if (opts->grain_write > 0
        && i % opts->grain_write == 0
        && opts->save_steps == 1) {

      FILE* out_step_file;
      char* step_fname;

      print_bits(size, size, *frame1, out_string);
      if (opts->save_flag == TMP_FILE) {
        asprintf(&step_fname, "%s/tmp_%i.step", opts->out_step_dir, i);
      }
      else {
        asprintf(&step_fname, "data_2d_%i/steps/out%s_%i.step",
                 states, rule_buf, i);
      }

      out_step_file = fopen(step_fname, "w+");
      free(step_fname);
      fprintf(out_step_file, "%s", out_string);
      fclose(out_step_file);
    }

    if (opts -> output_data == NO_OUTPUT) {
      continue;
    }

    print_bits(size, size, *frame1, out_string);
    /* memcpy(&dbl_pholder[i * ((size + 1) * size + 1)], */
           /* out_string, (size + 1) * size + 1); */

    if (i % opts->grain == 0) {
      last_compressed_size = compressed_size;
      last_cell_count = cell_count;

      print_bits(size, size, *frame1, out_string);
      compressed_size = compress_memory_size(out_string, (size + 1) * size);
      cell_count = count_cells(size, *frame1, states);

      /* Check if state has evolved from last time (stop mechanism) */
      if (compressed_size == last_compressed_size
          && flag == 1 && opts->early == EARLY) {
        printf("\n");
        break;
      }
      else if (compressed_size == last_compressed_size) {
        flag = 1;
      }

      /* if (opts->joint_complexity == 1) { */
      /*   /\* memcpy(&dbl_pholder[size * size + 1], *\/ */
      /*          /\* out_string, size * size + 1); *\/ */

      /*   int dbl_comp_size = */
      /*     compress_memory_size(dbl_pholder, i * ((size + 1) * size + 1)); */
      /*   /\* memcpy(dbl_pholder, out_string, (size + 1) * size + 1); *\/ */

      /*   if (i > 0) { */
      /*     size_sum = last_compressed_size + compressed_size; */
      /*     ratio2 = (size_sum - dbl_comp_size)/(float)size_sum; */
      /*     ratio = ( (last_compressed_size / (float)last_cell_count) + */
      /*               (compressed_size / (float)cell_count) ) / */
      /*       (dbl_comp_size / (float)(last_cell_count + cell_count)); */
      /*   } */
      /*   fprintf(out_file, "%i    %i    %f    %f    " */
      /*                     "%i    %i    %i\n", */
      /*           i, compressed_size, ratio, ratio2, cell_count, last_cell_count, */
      /*           dbl_comp_size); */
      /* } else { */
      /*   fprintf(out_file, "%i    %i\n", i, compressed_size); */
      /* } */

      printf("%i  ", compressed_size);
      fflush(stdout);
    }

    int offset_a = 2, offset_b = 1;


    if (i == steps - WINDOW) {
      print_bits(size, size, *frame1, out_string300);
      res300 = populate_map(map300, size, *frame1, offset_a, states);
      res300b = populate_map(map300b, size, *frame1, offset_b, states);

      memcpy(automat300, *frame1, size * size * sizeof(uint8_t));
    }
    if (i > (steps - WINDOW) && (i - (steps - WINDOW)) % W_STEP == 0) {
      memcpy(test_automata[((i - (steps - WINDOW)) / W_STEP) - 1],
             *frame1, size * size * sizeof(uint8_t));
    }

    if (i == steps - 51) {
      print_bits(size, size, *frame1, out_string50);
      res50 = populate_map(map50, size, *frame1, offset_a, states);
      res50b = populate_map(map50b, size, *frame1, offset_b, states);

      memcpy(automat50, *frame1, size * size * sizeof(uint8_t));
    }
    if (i == steps - 5) {
      print_bits(size, size, *frame1, out_string5);
      res5 = populate_map(map5, size, *frame1, offset_a, states);
      res5b = populate_map(map5b, size, *frame1, offset_b, states);

      memcpy(automat5, *frame1, size * size * sizeof(uint8_t));
    }
    if (i == steps - 1) {

      asprintf(&entrop_fname, "data_2d_%i/ent/ent%s.dat", states, rule_buf);
      entrop_file = fopen(entrop_fname, "w+");

      add_entropy_results_to_file(map300, size, *frame1, states, offset_a,
                          entrop_file, res300);
      add_entropy_results_to_file(map50, size, *frame1, states, offset_a,
                          entrop_file, res50);
      add_entropy_results_to_file(map5, size, *frame1, states, offset_a,
                          entrop_file, res5);
      fprintf(entrop_file, "\n");


      add_entropy_results_to_file(map300b, size, *frame1, states, offset_b,
                          entrop_file, res300b);
      add_entropy_results_to_file(map50b, size, *frame1, states, offset_b,
                          entrop_file, res50b);
      add_entropy_results_to_file(map5b, size, *frame1, states, offset_b,
                          entrop_file, res5b);
      fprintf(entrop_file, "\n");


      /* asprintf(&nn_fname, "data_2d_%i/nn/nn%s.dat", states, rule_buf); */
      /* nn_file = fopen(nn_fname, "w+"); */
      asprintf(&fisher_fname, "data_2d_%i/nn/fisher%s.dat", states, rule_buf);
      fisher_file = fopen(fisher_fname, "w+");

      network_result_t res = {1., 1., 1., 0.};
      network_opts_t n_opts = {10, 40, 3, MOMENTUM, DECAY, NO_FISHER, 1};

      for (int i = 4; i < 5; ++i) {
        n_opts.num_hid = 10;
        n_opts.offset = i;
        n_opts.fisher = NO_FISHER;

        train_nn_on_automaton(size, states, automat300, test_automata,
                              WINDOW / W_STEP, &n_opts, &res);

        /*add_nn_results_to_file(nn_file, &n_opts, &res, 50);*/
        fprintf(fisher_file, "%f", res.error_var);

        results->nn_tr_50 = res.error_var;
        results->nn_te_50 = 1.;
      }
    }
  }
  printf("\n");

  /* Cleanup before finishing */
  /* free(dbl_pholder); */

  free_map(map5);
  free_map(map300);
  free_map(map50);
  free_map(map5b);
  free_map(map300b);
  free_map(map50b);

  if (test_automata) {
    for (int i = 0; i < WINDOW / W_STEP; ++i) {
      free(test_automata[i]);
    }
    free(test_automata);
  }

  if (automat5)
    free(automat5);
  if (automat50)
    free(automat50);
  if (automat300)
    free(automat300);

  if (fname) {
    free(fname);
    fclose(out_file);
  }

  if (mult_time_file) {
    free(mult_time_fname);
    fclose(mult_time_file);
  }

  free(*frame1);
  free(frame1);
  free(*frame2);
  free(frame2);
  free(mask);

  if (entrop_file) {
    free(entrop_fname);
    fclose(entrop_file);
  }
  if (nn_file) {
    free(nn_fname);
    fclose(nn_file);
  }
  if (fisher_file) {
    free(fisher_fname);
    fclose(fisher_file);
  }
}
