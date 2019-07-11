#include <stdint.h>
#include "automaton/2d_automaton.h"

void iterative_search(int n_simulations, int input_flag,
                      long timesteps, uint64_t grule_size,
                      uint8_t* rule_array, char* rule_buf,
                      struct Options2D* opts);
