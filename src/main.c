#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "automaton/2d_automaton.h"
#include "automaton/rule.h"
#include "automaton/wolfram_automaton.h"
#include "utils/utils.h"
#include "search/genetic.h"


const char help[] = "Use with either 2d or 1d as first argument";

int main_2d(int argc, char** argv)
{
  char usage[] = "%s 2d [-n n_states] [-i input_rule] [-s size] [-t timesteps]"
    "[-g grain] [-c compress]";
  char one_input[] = "Give only one input, either -i rule or -f rule_file";

  extern char *optarg;
  extern int optind;
  int c, err = 0, input_flag = 0, search = 0;
  int timesteps = 1000, n_simulations = 1000;
  int compress_flag = 1;
  char* input_rule;
  char* input_fname = NULL;
  struct Options2D opts;
  opts.size = 256;
  opts.states = 2;
  opts.grain = 200;
  opts.horizon = 1;
  opts.grain_write = 200;
  opts.save_flag = STEP_FILE;
  opts.joint_complexity = 1;
  opts.early = EARLY;
  sprintf(opts.out_step_dir, "%s", "rule_gif");
  opts.noise_rate = 0.01;
  opts.noise_step = 133;
  opts.mask = NO_MASK;
  opts.output_data = OUTPUT;

  /* Optional arguments processing */
  while ((c = getopt(argc - 1, &argv[1], "n:i:s:t:g:c:z:f:mw:ero:q")) != -1)
    switch (c) {
    case 'r':
      search = 1;
      break;
    case 'n':
      opts.states = atoi(optarg);
      break;
    case 'i':
    case 'f':
      if (input_flag == 1) {
        fprintf(stderr, one_input, NULL);
        exit(1);
      }
      input_flag = 1;
      if (c == 'i') {
        input_rule = (char*) calloc(strlen(optarg), sizeof(char));
        strcpy(input_rule, optarg);
      }
      else {
        input_fname = (char*) calloc(strlen(optarg), sizeof(char));
        strcpy(input_fname, optarg);
      }
      break;
    case 's':
      opts.size = atol(optarg);
      break;
    case 't':
      timesteps = atoi(optarg);
      break;
    case 'g':
      opts.grain = atoi(optarg);
      break;
    case 'c':
      compress_flag = 0;
      break;
    case 'z':
      n_simulations = atoi(optarg);
      break;
    case 'm':
      opts.save_flag = TMP_FILE;
      opts.output_data = NO_OUTPUT;
      break;
    case 'w':
      opts.grain_write = atoi(optarg);
      break;
    case 'e':
      opts.early = NO_STOP;
      break;
    case 'o':
      sprintf(opts.out_step_dir, "%s", optarg);
      break;
    case 'q':
      opts.mask = MASK;
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

  const int side = 2 * opts.horizon + 1; /* Sidelength of neighborhood */
  const int neigh_size = side * side - 1; /* Total number of neighbors */
  /* Size of rule */
  const uint64_t grule_size = (int) pow(opts.states, neigh_size + 1);

  /* These arrays can be very big and are therefore allocated on the heap */
  char* rule_buf = malloc((grule_size + 1) * sizeof(char));
  uint8_t* rule_array = malloc(grule_size * sizeof(uint8_t));


  /* If input was given, it was either directly inline or via a file */
  if (input_fname) {
    FILE* rule_file;
    if ((rule_file = fopen(input_fname, "r"))) {
      uint64_t count = 0;
      while ((c = getc(rule_file)) != EOF && count < grule_size) {
        rule_array[count] = (uint8_t)(c - '0');
        rule_buf[count] = c;
        count++;
      }
      rule_buf[count] = '\0';
      if (count != grule_size) {
        fprintf(stderr, "Incorrect rule in file %s\n", input_fname);
        exit(1);
      }
    }
    else {
      fprintf(stderr, "Error while opening file %s\n", input_fname);
      exit(1);
    }
    free(input_fname);
  }
  /* Input given inline */
  else if (input_flag == 1) {
    build_rule_from_args(grule_size, rule_array, rule_buf, input_rule,
                         opts.states);
    free(input_rule);
  }

  /* If input rule was provided, and not starting a search: write steps for a
     given rule  */
  if (input_flag == 1 && search == 0) {

    results_nn_t res;
    /* This should not be necessary if the provided rule is already symmetric */
    /* TODO: Add possibility to work with either type */
    symmetrize_rule(grule_size, rule_array, opts.states, opts.horizon);

    make_map(&opts, rule_buf, 0);

    opts.save_steps = 1;
    process_rule(grule_size, rule_array, rule_buf, 0, timesteps, &opts, &res);
    free(rule_array);
    free(rule_buf);
    return 0;
  }

  /* The user wants to do a search */
  /* Iterative search */
  if (search == 1) {
    iterative_search(n_simulations, input_flag, timesteps, grule_size,
                     rule_array, rule_buf, &opts);
  }
  /* Generate compression plots for many rules */
  /* Regular search */
  else {
    for (int i = 0; i < n_simulations; ++i) {
      generate_general_rule(grule_size, rule_array, rule_buf,
                            opts.states, opts.horizon);

      make_map(&opts, rule_buf, i);

      results_nn_t res;
      opts.save_steps = 0;
      process_rule(grule_size, rule_array, rule_buf, 0, timesteps, &opts, &res);
    }
  }

  printf("\n");
  free(rule_array);
  free(rule_buf);
  return 0;
}

int main_1d(int argc, char** argv)
{
  extern char *optarg;
  extern int optind;
  int c, err = 0;
  size_t size = 256;
  int states = 2;
  struct Options1D options;
  options.init = ONE;
  options.timesteps = 256;
  options.grain = 50;
  options.write = NO_WRITE;
  options.radius = 1;

  char usage[] = "%s 1d [-s size] [-t timesteps]"
                 "[-n n_states] [-w neighborhood width]"
                 "[-i (one|random|random_small)] [-o]";
  char invalid_init[] = "Invalid value \"%s\" for init option."
                        " Must be one of \"one\", \"random\","
                        " \"random_small\"\n";

  while ((c = getopt(argc - 1, &argv[1], "s:t:n:r:i:og:")) != -1)
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
    case 'r':
      options.radius = atoi(optarg);
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
    case 'g':
      options.grain = atoi(optarg);
      break;
    case '?':
      err = 1;
      break;
    }

  if (err) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  size_t rule_size = (int)pow(states, 2 * options.radius + 1);
  uint8_t* rule = (uint8_t *)malloc(sizeof(uint8_t) * rule_size);

  time_t t;
  srand((unsigned)time(&t));

  uint64_t maxi = ((uint64_t)ipow(states, rule_size) > 300) ? 300:
    (uint64_t)pow(states, rule_size);

  for (uint64_t n = 0; n < maxi; n++) {
    for (size_t i = 0; i < rule_size; ++i) {
      if (maxi == 300) {
        rule[i] = (uint8_t)(rand() % states);
      }
      else {
        rule[i] = (uint8_t)((n / ipow(states, i)) % states);
      }
    }
    printf("%lu\t", rule_number(states, rule_size, rule));
    write_to_file(size, rule_size, rule, 0, &options, states);
    printf("\n");
  }

  free(rule);
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
