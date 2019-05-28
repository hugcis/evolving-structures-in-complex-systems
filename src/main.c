#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "2d_automaton.h"
#include "wolfram_automaton.h"
#include "utils.h"

const char help[] = "Use with either 2d or 1d as first argument";

int main_2d(int argc, char** argv)
{
  char usage[] = "%s 2d [-n n_states] [-i input_rule] [-s size] [-t timesteps]"
    "[-g grain] [-c compress]";


  extern char *optarg;
  extern int optind;
  int c, err = 0, input_flag = 0;
  int states = 2, grain = 100;
  int timesteps = 1000, horizon = 1, n_simulations = 1000;
  int compress_flag = 1;
  size_t size = 256;
  char* input_rule;

  /* Optional arguments processing */
  while ((c = getopt(argc - 1, &argv[1], "n:i:s:t:g:c:z:")) != -1)
    switch (c) {
    case 'n':
      states = atoi(optarg);
      break;
    case 'i':
      input_flag = 1;
      input_rule = calloc(strlen(optarg), sizeof(char));
      strcpy(input_rule, optarg);
      break;
    case 's':
      size = atol(optarg);
      break;
    case 't':
      timesteps = atoi(optarg);
      break;
    case 'g':
      grain = atoi(optarg);
      break;
    case 'c':
      compress_flag = 0;
      break;
    case 'z':
      n_simulations = atoi(optarg);
      break;
    case '?':
      err = 1;
      break;
    }

  if (err) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  time_t t;
  srand((unsigned)time(&t));

  const int side = 2 * horizon + 1; /* Sidelength of neighborhood */
  const int neigh_size = side * side - 1; /* Total number of neighbors */
  /* Size of rule */
  const uint64_t grule_size = (int) pow(states, neigh_size + 1);

  char rule_buf[grule_size + 1];
  uint8_t rule_array[grule_size];

  FILE* dic_file;
  char* fname;


  /* If input rule was provided, write steps for a given rule  */
  if (input_flag == 1) {
    build_rule_from_args(grule_size, rule_array, rule_buf, input_rule, states);
    free(input_rule);

    symmetrize_rule(grule_size, rule_array, states, horizon);

    asprintf(&fname, "data_2d_%i/map/%lu.map", states, hash(rule_buf));
    dic_file = fopen(fname, "w+");
    fprintf(dic_file, "%s", rule_buf);
    fclose(dic_file);

    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("Rule %s\n", rule_buf);

    process_rule(grule_size, rule_array, rule_buf,
                 1, 1, 0, timesteps, states, horizon, size, grain);
    return 0;
  }


  /* Generate compression plots for many rules */
  for (int i = 0; i < n_simulations; ++i) {
    /* generate_totalistic_rule(rule_array, rule_buf); */
    generate_general_rule(grule_size, rule_array, rule_buf, states, horizon);

    asprintf(&fname, "data_2d_%i/map/%lu.map", states, hash(rule_buf));
    dic_file = fopen(fname, "w+");
    fprintf(dic_file, "%s", rule_buf);
    fclose(dic_file);
    free(fname);

    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("%i: Rule %s\n", i, rule_buf);
    fflush(stdout);

    process_rule(grule_size, rule_array, rule_buf,
                 0, 1, 0, timesteps, states, horizon, size, grain);
  }
  printf("\n");
  return 0;
}

int main_1d(int argc, char** argv)
{
  extern char *optarg;
  extern int optind;
  int c, err = 0;
  size_t size = 256, neighbor = 3;
  int states = 2;
  struct Options1D options;
  options.init = ONE;
  options.timesteps = 256;
  options.write = NO_WRITE;

  char usage[] = "%s 1d [-s size] [-t timesteps]"
                 "[-n n_states] [-w neighborhood width]"
                 "[-i (one|random|random_small)] [-o]";
  char invalid_init[] = "Invalid value \"%s\" for init option."
                        " Must be one of \"one\", \"random\","
                        " \"random_small\"\n";

  while ((c = getopt(argc - 1, &argv[1], "s:t:n:w:i:o")) != -1)
    switch (c) {
    case 's':
      size = atol(optarg);
      break;
    case 't':
      options.timesteps = atoi(optarg);
      break;
    case 'n':
      states = atoi(optarg);
      break;
    case 'w':
      neighbor = atoi(optarg);
      break;
    case 'o':
      options.write = WRITE_STEP;
      break;
    case 'i':
      if (strcmp("one", optarg) == 0) {
        options.init = ONE;
      }
      else if (strcmp("random", optarg) == 0) {
        options.init = RANDOM;
      }
      else if (strcmp("random_small", optarg) == 0) {
        options.init = RAND_SMALL;
      }
      else {
        printf(invalid_init, optarg);
        exit(1);
      }
      break;
    case '?':
      err = 1;
      break;
    }

  if (err) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  size_t rule_size = (int)pow(states, neighbor);
  uint8_t rule[rule_size];

  time_t t;
  srand((unsigned)time(&t));
  fprintf(stderr, "%"PRIu64, (uint64_t)pow(states, rule_size));
  for (uint64_t n = 0; n < (uint64_t)pow(states, rule_size); n++) {
    for (size_t i = 0; i < rule_size; ++i) {
      rule[i] = (uint8_t)((n / ipow(states, i)) % states);
    }
    printf("%lu\t", rule_number(states, rule_size, rule));
    write_to_file(size, rule_size, rule, 0, &options, states);
    printf("\n");
  }

  return 0;
}

int main(int argc, char** argv)
{
  char base_usage[] = "%s [\"1d\" or \"2d\"]";
  /* No first argument was passed */
  if (argc < 2) {
    printf("%s\n", help);
    return 0;
  }

  if (strcmp("2d", argv[1]) == 0) {
    return main_2d(argc, argv);
  }
  else if (strcmp("1d", argv[1]) == 0) {
    return main_1d(argc, argv);
  }
  else {
    printf("%s\n", help);
    fprintf(stderr, base_usage, argv[0]);
  }
}
