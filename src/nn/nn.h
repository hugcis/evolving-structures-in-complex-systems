#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct network_result_s
{
  double train_error;
  double test_error;
} network_result_t;

void train_nn_on_automaton(size_t size, int, uint8_t[size][size],
                           uint8_t[size][size], int, FILE*,
                           network_result_t*);
