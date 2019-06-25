#include <inttypes.h>
#include "automaton/2d_automaton.h"

#ifndef RULE_H
#define RULE_H

void populate_buf(uint64_t, uint8_t*, char*);

void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char*, int);

void symmetrize_rule(uint64_t grule_size,
                     uint8_t rule_array[grule_size], int, int);

void perturb_rule(uint64_t grule_size,
                  uint8_t rule_array[grule_size],
                  char rule_buf[grule_size + 1],
                  int, int, double);

void make_map(struct Options2D*, char*, int);

#endif /* RULE_H */
