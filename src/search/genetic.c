#include <string.h>
#include "nn/nn.h"
#include "automaton/2d_automaton.h"
#include "automaton/rule.h"

#define RATE 0.005

typedef struct val_idx_s
{
  double value;
  int index;
} val_idx_t;

int cmp(const void *a, const void *b)
{
  val_idx_t *a1 = (val_idx_t *)a;
  val_idx_t *a2 = (val_idx_t *)b;
  if ((*a1).value > (*a2).value)
    return -1;
  else if ((*a1).value < (*a2).value)
    return 1;
  else
    return 0;
}

double compute_score(results_nn_t* res)
{
  /* For now, hardcoded weighting coefs for score */
  /* double a = 7./10., b = 7./10., c = -4./10.; */

  /* Either premature stop or some training did go to 0 */
  if (res->nn_te_300 * res->nn_te_50 * res->nn_te_5 == 0) {
    return 0.;
  }
  /* Everything went well */
  else {
    return res->nn_te_300 / res->nn_tr_300;
  }
}


void iterative_search(int n_simulations, int input_flag,
                      long timesteps, uint64_t grule_size,
                      uint8_t* rule_array, char* rule_buf,
                      struct Options2D* opts)
{
  FILE* genealogy_file;
  char* genealogy_fname;

  int population_size = 5;
  int n_children = 10;
  results_nn_t res = {0, 0, 0, 0, 0, 0};

  uint8_t** population = (uint8_t**) malloc(sizeof(uint8_t *)
                                            * population_size);
  uint8_t** tmp_pop = (uint8_t**) malloc(sizeof(uint8_t *)
                                         * population_size);

  uint8_t** children = (uint8_t**) malloc(sizeof(uint8_t *)
                                          * population_size
                                          * n_children);
  val_idx_t* results =
    (val_idx_t*) malloc(sizeof(val_idx_t) *
                         population_size * (n_children + 1));

  for (int i = 0; i < n_simulations; ++i) {
    /* Initialize rule */
    if (i == 0 && input_flag == 0) {
      generate_general_rule(grule_size, rule_array, rule_buf,
                            opts->states, opts->horizon);
    }
    /* Initialize search */
    if (i == 0) {
      asprintf(&genealogy_fname, "data_2d_%i/nn/%lu.gen", opts->states,
               hash(rule_buf));
      genealogy_file = fopen(genealogy_fname, "w+");

      /* Initialization done in 2 times */
      /* First: add the base rule in the initial population */
      population[0] = (uint8_t*) malloc(sizeof(uint8_t) * grule_size);
      tmp_pop[0] = (uint8_t*) malloc(sizeof(uint8_t) * grule_size);

      memcpy(population[0], rule_array, sizeof(uint8_t) * grule_size);
      for (int d = 0; d < n_children; ++d) {
        children[d] =
          (uint8_t *) malloc(sizeof(uint8_t) * grule_size);
      }

      /* Second: add some children of the base rule for the rest of the
         population */
      for (int k = 1; k < population_size; ++k) {
        /* Initialize a population from the base rule */
        population[k] = (uint8_t*) malloc(sizeof(uint8_t) * grule_size);
        tmp_pop[k] = (uint8_t*) malloc(sizeof(uint8_t) * grule_size);

        memcpy(population[k], rule_array, sizeof(uint8_t) * grule_size);


        perturb_rule(grule_size, population[k], rule_buf, opts->states,
                     opts->horizon, 10 * RATE);

        /* Allocate space for childrenrules */
        for (int d = 0; d < n_children; ++d) {
          children[k * n_children + d] =
            (uint8_t *) malloc(sizeof(uint8_t) * grule_size);
        }
      }
    }

    for (int k = 0; k < population_size; ++k) {
      /* Compute results for the parents */

      populate_buf(grule_size, population[k], rule_buf);
      make_map(opts, rule_buf, i);

      res.nn_tr_5 = 1.;
      res.nn_tr_50 = 1.;
      res.nn_tr_300 = 1.;
      res.nn_te_5 = 1.;
      res.nn_te_50 = 1.;
      res.nn_te_300 = 1.;

      process_rule(grule_size, population[k],
                   rule_buf, 0, timesteps,
                   opts, &res);

      results[k * (n_children + 1) + n_children].index =
        k * (n_children + 1) + n_children;
      results[k * (n_children + 1) + n_children].value = compute_score(&res);

      /* Compute results for children */
      for (int d = 0; d < n_children; ++d) {
        int rule_A = rand() % population_size;
        int rule_B = rand() % population_size;

        cross_breed(grule_size, population[rule_A], population[rule_B],
                    children[k * n_children + d], rule_buf, .5, opts->horizon,
                    opts->states);

        make_map(opts, rule_buf, i);

        res.nn_tr_5 = 1.;
        res.nn_tr_50 = 1.;
        res.nn_tr_300 = 1.;
        res.nn_te_5 = 1.;
        res.nn_te_50 = 1.;
        res.nn_te_300 = 1.;

        process_rule(grule_size, children[k * n_children + d],
                     rule_buf, 0, timesteps, opts, &res);

        results[k * (n_children + 1) + d].index = k * (n_children + 1) + d;
        results[k * (n_children + 1) + d].value = compute_score(&res);
      }
    }

    /* Keep the best `population_size` results */
    qsort(results,
          (n_children + 1) * population_size,
          sizeof(results[0]), cmp);

    fprintf(stdout, "Best scores: ");
    for (int k = 0; k < population_size; ++k) {


      int index = results[k].index;

      if ((index % (n_children + 1)) == n_children) {
        populate_buf(grule_size,
                     population[index / (n_children + 1)],
                     rule_buf);

        fprintf(stdout, "p:%lu %f  ",
                hash(rule_buf),
                results[k].value);

        memcpy(tmp_pop[k],
               population[index / (n_children + 1)],
               sizeof(uint8_t) * grule_size);
      }
      else {
        int child_index = (index % (n_children + 1)) +
          n_children * (index / (n_children + 1));
        populate_buf(grule_size,
                     children[child_index],
                     rule_buf);

        fprintf(stdout, "c:%lu %f  ",
                hash(rule_buf),
                results[k].value);

        memcpy(tmp_pop[k],
               children[child_index],
               sizeof(uint8_t) * grule_size);
      }
    }
    fprintf(stdout, "\n");

    for (int k = 0; k < population_size; ++ k) {
      memcpy(population[k], tmp_pop[k], sizeof(uint8_t) * grule_size);

      populate_buf(grule_size, population[k], rule_buf);
      fprintf(genealogy_file, "%lu:%f\t",
              hash(rule_buf),
              results[k].value);
    }
    fprintf(genealogy_file, "\n");
    fflush(genealogy_file);
  }

  /* Free data structures */
  for (int k = 0; k < population_size; ++k) {
    free(population[k]);
    free(tmp_pop[k]);
    for (int d = 0; d < n_children; ++d) {
      free(children[k * n_children + d]);
    }
  }
  fclose(genealogy_file);
  free(results);
  free(tmp_pop);
  free(genealogy_fname);
  free(population);
  free(children);

}
