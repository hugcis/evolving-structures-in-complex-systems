#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "utils/utils.h"
#include "rule.h"

#define DIRICHLET 1
#define D_CONST 0.8

void populate_buf(uint64_t grule_size, uint8_t* rule_array, char* rule_buf)
{
  /* Populate buffer with newly created rule */
  for (uint64_t v = 0; v < grule_size; v++) {
    sprintf(&rule_buf[v], "%"PRIu8, rule_array[v]);
  }
  rule_buf[grule_size] = '\0';
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
 * Simple double comparison function used for sorting in the rule generation
 * function.
 */
int comp (const void * elem1, const void * elem2)
{
  float f = *((float*)elem1);
  float s = *((float*)elem2);
  return (f > s) - (f < s);
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

  /* Choose lambda parameter at random as well as the proportion of */
  /* transitions to other states */

  double alphas[states];
  double theta[states], lambda[states - 1], rand_num;

  /* This method samples the transition probability simplex according to a */
  /* Dirichlet distribution with parameter D_CONST */
  if (DIRICHLET == 1) {
    const gsl_rng_type * T;

    gsl_rng * r;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    time_t t;
    gsl_rng_set(r, (unsigned long)time(&t));

    for (int i = 0; i < states; ++i) {
      alphas[i] = D_CONST;
    }
    gsl_ran_dirichlet(r, states, alphas, theta);
    lambda[0] = theta[0];
    for (int i = 1; i < states - 1; ++i) {
      lambda[i] = theta[i] + lambda[i - 1];
    }

    gsl_rng_free(r);
  }
  /* Second method that samples the simplex uniformly */
  else {
    for (int i = 0; i < states - 1; ++i) {
      lambda[i] = ((double)rand() / (double)((unsigned)RAND_MAX + 1));
    }
    qsort(lambda, sizeof(lambda)/sizeof(*lambda), sizeof(*lambda), comp);
  }

  for (uint64_t v = 0; v < grule_size; v++) {
    /* Assign the rule to the first state that passes the test */
    rand_num = (double)rand() / (double)((unsigned)RAND_MAX + 1);
    inc = 0;
    while (lambda[inc] < rand_num && inc < states - 1) {
      ++inc;
    }
    rule_array[v] = (uint8_t)inc;
  }

  symmetrize_rule(grule_size, rule_array, states, horizon);

  populate_buf(grule_size, rule_array, rule_buf);
}


void perturb_rule(uint64_t grule_size,
                  uint8_t rule_array[grule_size],
                  char rule_buf[grule_size + 1],
                  int states, int horizon, double rate)
{
  for (uint64_t v = 0; v < grule_size; ++v) {
    /* Perturb transisition outcome with probability rate */
    double rand_num = (double)rand() / (double)((unsigned)RAND_MAX + 1);
    if (rand_num < rate) {
      rule_array[v] += (1 + rand() % (states - 1));
      rule_array[v] %= states;
    }
  }

  symmetrize_rule(grule_size, rule_array, states, horizon);
  populate_buf(grule_size, rule_array, rule_buf);
}

void make_map(struct Options2D* opts, char* rule_buf, int step)
{
  FILE* dic_file;
  char* fname;

  asprintf(&fname, "data_2d_%i/map/%lu.map", opts->states, hash(rule_buf));
  dic_file = fopen(fname, "w+");
  fprintf(dic_file, "%s", rule_buf);
  fclose(dic_file);
  free(fname);

  sprintf(rule_buf, "%lu", hash(rule_buf));
  printf("%i: Rule %s\n", step, rule_buf);

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
