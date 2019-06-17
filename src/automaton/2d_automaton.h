#include <math.h>
#include <stdlib.h>
#include <inttypes.h>

#ifndef TWOD_AUTOMATON_H /* Include guard */
#define TWOD_AUTOMATON_H

enum WriteStepMode { TMP_FILE, STEP_FILE };

unsigned long hash(char*);

void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char*, int);

void symmetrize_rule(uint64_t grule_size,
                     uint8_t rule_array[grule_size], int, int);

void process_rule(uint64_t grule_size,
                  uint8_t rule[grule_size],
                  char rule_buf[], int, int,
                  int, long, int, int, size_t,
                  int, enum WriteStepMode);

void update_step_general(uint64_t grule_size, size_t size,
                         uint8_t [size][size],
                         uint8_t rule[grule_size],
                         uint8_t [size][size], int, int);

void generate_general_rule(uint64_t grule_size,
                           uint8_t rule_array[grule_size],
                           char rule_buf[grule_size + 1], int, int);

#endif // TWOD_AUTOMATON_H
