/**
 * @file
 * @brief All 2D automaton functions.
 */
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>

#ifndef TWOD_AUTOMATON_H /* Include guard */
#define TWOD_AUTOMATON_H

enum WriteStepMode { TMP_FILE, STEP_FILE };
enum EarlyStop { EARLY, NO_STOP };
enum MaskEnum { MASK, NO_MASK };
enum DataOutput { OUTPUT, NO_OUTPUT };

/** A set of options to pass for generating and processing an automaton from a
 *  rule.
 */
struct Options2D {
  char data_dir_name[1024];
  size_t size; /**< The size of the square grid. */
  int grain_write; /**< The grain at which write operations are executed. */
  int grain; /**< The grain at which compression is done. */
  enum WriteStepMode save_flag; /**< Wether to write output to a temporary
                                   file or step file */
  char out_step_dir[1024]; /**< Output dir for the temporary file */
  int states;
  int horizon;
  int joint_complexity;
  int save_steps;
  enum EarlyStop early;
  int noise_step;
  double noise_rate;
  enum MaskEnum mask;
  enum DataOutput output_data; /**< Wether to compute metrics on the
                                  simulations */
  char init_pattern_fname[1024];
  FILE* init_pattern_file;
  long init_type; /**< Size of the random initialization zone (-1 for full) */
};

typedef struct results_nn_s
{
  double nn_tr_300;
  double nn_te_300;
  double nn_tr_50;
  double nn_te_50;
  double nn_tr_5;
  double nn_te_5;
} results_nn_t;


typedef void (*ProcessF)(size_t size, uint8_t*,
                         uint8_t[], uint8_t*,
                         int, uint32_t*);

unsigned long hash(char*);

/**
 * @brief Main 2D rule processing function.
 *
 * Given a rule, create and simulate the corresponding automaton.
 */
void process_rule(uint64_t grule_size,
                  uint8_t rule[grule_size],
                  char rule_buf[],
                  long,
                  struct Options2D*,
                  results_nn_t*);


void generate_general_rule(uint64_t grule_size,
                           uint8_t rule_array[grule_size],
                           char rule_buf[grule_size + 1], int, int);

#endif // TWOD_AUTOMATON_H
