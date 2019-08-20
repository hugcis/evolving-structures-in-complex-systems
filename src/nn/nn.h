#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

enum OptimType { MOMENTUM, ADAM, NESTEROV, SGD };
enum LRDecay { NO_DECAY, DECAY };
enum FisherInfo { FISHER, NO_FISHER };

typedef struct network_opts_s
{
  int num_hid;
  int max_epoch;
  int offset;
  enum OptimType optim_type;
  enum LRDecay decay;
  enum FisherInfo fisher;
} network_opts_t;

typedef struct network_result_s
{
  double train_error;
  double test_error;
  double fisher_info;
} network_result_t;

void train_nn_on_automaton(size_t, int,
                           uint8_t*,
                           uint8_t*,
                           network_opts_t*,
                           network_result_t*);

