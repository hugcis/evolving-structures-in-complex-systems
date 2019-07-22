#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct network_opts_s
{
  int num_hid;
  int max_epoch;
  int offset;
} network_opts_t;

typedef struct network_result_s
{
  double train_error;
  double test_error;
} network_result_t;

void train_nn_on_automaton(size_t size, int, uint8_t*,
                           uint8_t*,
                           network_opts_t*,
                           network_result_t*);
