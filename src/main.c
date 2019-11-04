#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include "automaton/2d_automaton.h"
#include "automaton/rule.h"
#include "automaton/wolfram_automaton.h"
#include "utils/utils.h"
#include "search/genetic.h"

#define MAJOR 0
#define MINOR 1

extern int errno ;

const char help[] = "Use with either 2d or 1d as first argument";

int main_2d(int argc, char** argv)
{
  char usage[] = "2D Automaton.\n\n\
Usage:\n\
    %s 2d\n\
\n\
Options:\n\
    -h --help               Print this message.\n\
    -v --version            Print version number.\n\
    -n --n_states=<n>       Number of states [default: 2].\n\
    -s --size=<s>           Size of the automaton [default: 256].\n\
    -t --timesteps=<ts>     Number of timesteps of the simulation.\n\
    -g --grain=<g>          Grain of compression operations [default: 200].\n\
    -w --write_grain=<w>    Interval between state writes [default: 200].\n\
    -z --n_simulations=<n>  Number of simulations to run [default: 1000].\n\
    -o --output=<dir>       Output dirname for steps [default: 'rule_gif'].\n\
    -b --init_type=<i>      Size of random disordered initial zone \n\
                            (-1 for fully random) [default: -1].\n\
    -i --input_rule=<rule>  Inline rule input.\n\
    -f --input_file=<fname> File containing the desired rule.\n\
    -j --pattern=<fname>    Pattern filename.\n\
    -r --search             Start a search.\n\
    -m --temp_output        Do a step-only output for visualization.\n\
    -e --no_early_stopping  Disable stopping when periodic.\n\
    -q --masking            Enable masking of the input.\n\
    -c --compress           Disable compression of outputs.\n";

  char one_input[] = "Provide only one input, either -i rule (for inline) or -f"
    " rule_file (for a file).\n";
  char wrong_transitions[] = "Incorrect rule in file %s: %"PRIu64
    " transitions found but %"PRIu64" were expected.\n";
  char too_large_init[] = "Initialization zone size is too large: %l was given"
    " but size is %lu.\n";

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
  opts.init_type = -1;
  opts.init_pattern_file = NULL;

  while (1) {
    static struct option long_options[] = {
       {"help", no_argument, 0, 'h'},
       {"version", no_argument, 0, 'v'},
       {"search", no_argument, 0, 'r'},
       {"temp_output", no_argument, 0, 'm'},
       {"early_stopping", no_argument, 0, 'e'},
       {"masking", no_argument, 0, 'q'},
       {"compress", no_argument, 0, 'c'},
       {"n_states", required_argument, 0, 'n'},
       {"input_rule", required_argument, 0, 'i'},
       {"input_file", required_argument, 0, 'f'},
       {"size", required_argument, 0, 's'},
       {"timesteps", required_argument, 0, 't'},
       {"grain", required_argument, 0, 'g'},
       {"n_simulations", required_argument, 0, 'z'},
       {"write_grain", required_argument, 0, 'w'},
       {"output", required_argument, 0, 'o'},
       {"pattern", required_argument, 0, 'j'},
       {"init_type", required_argument, 0, 'b'},
       {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc - 1, &argv[1],
                     "hvn:i:s:t:g:cz:f:mw:ero:qj:",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

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
        exit(EXIT_FAILURE);
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
    case 'j':
      opts.init_pattern_file = fopen(optarg, "r");
      if (opts.init_pattern_file == NULL) {
        fprintf(stderr, "Error opening file %s\n", optarg);
        err = 1;
      }
      break;
    case 'b':
      opts.init_type = atol(optarg);
      break;
    case 'h':
      fprintf(stdout, usage, argv[0]);
      exit(EXIT_SUCCESS);
      break;
    case 'v':
      fprintf(stdout, "Version %i.%i\n", MAJOR, MINOR);
      exit(EXIT_SUCCESS);
      break;
    case '?':
      err = 1;
      break;
    }

    if (err) {
      fprintf(stderr, usage, argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  if (opts.init_type >= (long) opts.size) {
    fprintf(stderr, too_large_init, opts.init_type, opts.size);
    exit(EXIT_FAILURE);
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
      /* Read the file and put each number in rule_array */
      while ((c = getc(rule_file)) != EOF && count < grule_size) {
        rule_array[count] = (uint8_t)(c - '0');
        rule_buf[count] = c;
        count++;
      }
      rule_buf[count] = '\0'; /* Terminate string buffer */

      /* Incorrect number of transitions */
      if (count != grule_size) {
        fprintf(stderr, wrong_transitions, input_fname,
                count, grule_size);
        exit(EXIT_FAILURE);
      }
    }
    else {
      fprintf(stderr, "Error while opening file %s: %s\n", input_fname,
              strerror(errno));
      exit(EXIT_FAILURE);
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
  return EXIT_SUCCESS;
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
        exit(EXIT_FAILURE);
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
    exit(EXIT_FAILURE);
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
