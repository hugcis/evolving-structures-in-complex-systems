#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "nn.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#define AVERAGE_FISHER 100

/**
 * A function for generating random numbers according to a N(mu, sigma) Gaussian
 * distribution.
 */
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1) {
      call = !call;
      return (mu + sigma * (double) X2);
  }
  do {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
  } while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

/**
 * Shorthand to generate normally distributed random numbers.
 */
double rand_normal()
{
  return randn(0.0, 1.0);
}

/**
 * @brief Shuffle the index array.
 */
void shuffle_index(size_t num_pattern, size_t random_idx[num_pattern])
{
  size_t p, np, op;
  for (p = 0; p < num_pattern; ++p) {
    random_idx[p] = p;
  }
  for (p = 0; p < num_pattern; ++p) {
    np = p + ((double)rand()/((double)RAND_MAX+1)) * (num_pattern - 1 - p);
    op = random_idx[p];
    random_idx[p] = random_idx[np];
    random_idx[np] = op;
  }
}

/**
 * This function fills the input and target vector with the list of training
 * examples from automaton.
 */
void fill_input_target(size_t size, double* input, uint8_t* target,
                       uint8_t* automaton, int offset, int states)
{
  size_t index;
  int counter;
  int side = 2 * offset + 1;
  int num_input = states * (side * side - 1);
  uint8_t val;

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      index = i * size + j;
      counter = 1;

      /* Add bias in the main vector */
      input[index * (num_input + 1)] = 1.0;

      for (int a = -offset; a < offset + 1; ++a) {
        for (int b = -offset; b < offset + 1; ++b) {
          if (a != 0 || b != 0) {  /* Don't take index i,j */
            val = automaton[((i + a + size) % size) * size
                            + ((j + b + size) % size)];
            for (uint8_t s = 0; s < states; ++s) {
              input[index * (num_input + 1) + counter] = (val == s) ? 1.: 0.;
              counter++;
            }
          }
        }
      }

      target[index] = automaton[i * size + j];
    }
  }
}

void init_weights(int num_hidden, int num_input, int num_output,
                  double* delta_w_ih, double* weight_ih,
                  double* delta_w_ho, double* weight_ho)
{
  int i, j, k;
  /* Initialize weights input -> hidden */
  for (j = 0; j < num_hidden; ++j) {
    delta_w_ih[j] = 0.0;
    weight_ih[j] = 0.0;
    for (i = 1; i < num_input + 1; ++i) {
      delta_w_ih[i * num_hidden + j] = 0.0;
      weight_ih[i * num_hidden + j] = rand_normal() *
        sqrt(1 / (double)(num_input));
    }
  }

  /* Initialize weights hidden -> output */
  for (k = 0; k < num_output; ++k) {
    delta_w_ho[k] = 0.0;
    weight_ho[k] = 0.0;
    for (j = 1; j < num_hidden + 1; ++j) {
      delta_w_ho[j * num_output + k] = 0.0;
      weight_ho[j * num_output + k] = rand_normal() *
        sqrt(1 / (double)(num_hidden));
    }
  }
}

void forward(double* input, double* output,
             int num_hidden, int num_pattern,
             int num_input, int num_output,
             double* hidden,
             double* hidden_bias,
             double* weight_ih,
             double* weight_ho)
{
  int j, k, p;
  double max_out, agg_out = 0.0;

  /* Compute hidden activations */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              num_pattern, num_hidden, num_input + 1, 1.0,
              input, num_input + 1, weight_ih, num_hidden,
              0.0, hidden, num_hidden);

  /* ReLU non-linearity */
  for (p = 0; p < num_pattern; ++p) {
    hidden_bias[p * (num_hidden + 1)] = 1.0;
    for (j = 1; j < num_hidden + 1; ++j) {
      hidden_bias[p * (num_hidden + 1) + j] =
        (hidden[p * num_hidden + j - 1] > 0.0) ?
        hidden[p * num_hidden + j - 1]: 0.0;
    }
  }

  /* Compute output unit activations */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              num_pattern, num_output, num_hidden + 1, 1.0,
              hidden_bias, num_hidden + 1, weight_ho, num_output,
              0.0, output, num_output);

  /* Compute softmax of output */
  for (p = 0; p < num_pattern; ++p) {
    max_out = output[p * num_output];
    for (k = 1; k < num_output; ++k) {
      if (output[p * num_output + k] > max_out) {
        max_out = output[p * num_output + k];
      }
    }
    agg_out = 0.0;
    for (k = 0; k < num_output; ++k) {
      output[p * num_output + k] = exp(output[p * num_output + k] - max_out);
      agg_out += output[p * num_output + k];
    }
    for (k = 0; k < num_output; ++k) {
      output[p * num_output + k] /= agg_out;
    }
  }
}

void update_weights(int num_input, int num_hidden, int num_output,
                    double eta,
                    double* batch_error, double reg, double alpha,
                    double* weight_ih, double* delta_w_ih,
                    double* delta_w_ih_prev,
                    double* weight_ho, double* delta_w_ho,
                    double* delta_w_ho_prev)
{
  int i, j, k;
  /* Update weights with average gradient */
  for (i = 0; i < num_input + 1; ++i) {
    for (j = 0; j < num_hidden; ++j) {
      if (reg > 0.) {
        *batch_error += 0.5 * reg *
          weight_ih[i * num_hidden + j] * weight_ih[i * num_hidden + j] ;
        weight_ih[i * num_hidden + j] += (1 + alpha) *
          reg * weight_ih[i * num_hidden + j];
      }

      weight_ih[i * num_hidden + j] -= (1 + alpha) *
        (eta * delta_w_ih[i * num_hidden + j]);

      if (alpha > 0.) {
        weight_ih[i * num_hidden + j] += (-alpha) *
          delta_w_ih_prev[i * num_hidden + j];
      }
    }
  }
  for (j = 0; j < num_hidden + 1; ++j) {
    for (k = 0; k < num_output; ++k) {
      if (reg > 0.) {
        *batch_error += 0.5 * reg *
          weight_ho[j * num_output + k] * weight_ho[j * num_output + k];
        weight_ho[j * num_output + k] += (1 + alpha) *
        reg * weight_ho[j * num_output + k];
      }

      weight_ho[j * num_output + k] -= (1 + alpha) *
        (eta * delta_w_ho[j * num_output + k]);

      if (alpha > 0.) {
        weight_ho[j * num_output + k] += (-alpha) *
          delta_w_ho_prev[j * num_output + k];
      }
    }
  }
}

void compute_batch_gradients(int base_index, double alpha,
                             int batch_size, int num_input, int num_output,
                             int num_hidden, size_t* random_idx,
                             double* output,
                             uint8_t* target, double* delta_output,
                             double* delta_w_ho, double* hidden_bias,
                             double* weight_ho, double* delta_h,
                             double* delta_w_ih, double* input,
                             double* delta_w_ih_prev,
                             double* delta_w_ho_prev,
                             network_opts_t* opts)
{
  int i, j, k, p;

  /* Gradient initialization (Nesterov momentum) */
  if (opts->optim_type == NESTEROV) {
    memcpy(delta_w_ih_prev, delta_w_ih,
           sizeof(double) * num_hidden * (num_input + 1));
    for (i = 0; i < num_input + 1; ++i) {
      for (j = 0; j < num_hidden; ++j) {
        delta_w_ih[i * num_hidden + j] *= alpha;
      }
    }

    memcpy(delta_w_ho_prev, delta_w_ho,
           sizeof(double) * num_output * (num_hidden + 1));
    for (j = 0; j < num_hidden + 1; ++j) {
      for (k = 0; k < num_output; ++k) {
        delta_w_ho[j * num_output + k] *= alpha;
      }
    }
  }

  for (int b = 0; b < batch_size; ++b) {
    p = random_idx[(base_index + b)];

    /* Backpropagation */
    for (k = 0; k < num_output; ++k) {
      /* Output gradients */
      delta_output[b * num_output + k] =
        (output[b * num_output + k] - ((k == target[p])? 1.0: 0.0));
      delta_output[b * num_output + k] /= (double) batch_size;
    }
  }

  /* Dot product hidden_bias x delta_output stored in delta_w_ho */
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
              num_hidden + 1, num_output, batch_size,
              1 , hidden_bias, num_hidden + 1,
              delta_output, num_output,
              (opts->optim_type == NESTEROV)? 1.0: 0.0,
              delta_w_ho, num_output);

  /* Dot product weight_ho x delta_output (we don't count the bias in
     weight_ho) stored in delta_h */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
              batch_size, num_hidden, num_output,
              1.0, delta_output, num_output,
              &weight_ho[num_output], num_output,
              0.0, delta_h, num_hidden);

  /* Hidden layer non linearity gradients */
  for (int b = 0; b < batch_size; ++b) {
    for (j = 0; j < num_hidden; ++j) {
      delta_h[b * num_hidden + j] *=
        ((hidden_bias[b * (num_hidden + 1) + j + 1] > 0) ? 1.0: 0.0);
    }
  }

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
              num_input + 1, num_hidden, batch_size,
              1 , input, num_input + 1,
              delta_h, num_hidden,
              (opts->optim_type == NESTEROV)? 1.0: 0.0,
              delta_w_ih, num_hidden);

}

double compute_loss(int base_index, int batch_size,
                    size_t* random_idx,
                    int num_output, double* output,
                    uint8_t* target)
{
  double batch_error = 0.0;
  int p;
  for (int b = 0; b < batch_size; ++b) {
    p = random_idx[(base_index + b)];
    /* Compute loss */
    batch_error +=
      - log((output[b * num_output + target[p]] > 0) ?
            output[b * num_output + target[p]]: DBL_MIN);
  }
  batch_error /= (double) batch_size;
  return batch_error;
}


double compute_fisher(int states, int neighbors,
                      int num_pattern, int num_output,
                      int num_hidden, int num_input,
                      double* input,
                      double* weight_ih, double* weight_ho)
{
  double fisher_information = 0.0;
  int neighbor_n, base, index, delta;

  /* Build the secondary input that will be perturbed */
  double* perturbed_input =
    (double*) malloc(sizeof(double) * num_pattern * (num_input + 1));

  /* Allocate placeholders in the function to make it more
     adaptable to the input dataset. */
  double* output =
    (double*) malloc(sizeof(double) * num_pattern * num_output);
  double* output_pert =
    (double*) malloc(sizeof(double) * num_pattern * num_output);
  double* hidden =
    (double*) malloc(sizeof(double) * num_pattern * num_hidden);
  double* hidden_bias =
    (double*) malloc(sizeof(double) * num_pattern * (num_hidden + 1));

  /* Compute the output of the network on the base dataset */
  forward(input, output, num_hidden, num_pattern, num_input, num_output,
          hidden, hidden_bias, weight_ih, weight_ho);

  for (int n = 0; n < AVERAGE_FISHER; ++n) {

    /* Perturb second input */
    memcpy(perturbed_input, input,
          sizeof(double) * num_pattern * (num_input + 1));
    for (int i = 0; i < num_pattern; ++i) {
      /* Choose index of the input cell to perturb */
      neighbor_n = rand() % neighbors;
      index = i * (num_input + 1) + neighbor_n * states;

      for (int t = 0; t < states; ++t) {
        base = 1 + (rand() % (states - 1));
        perturbed_input[index + ((base + t) % states)] =
          input[index + t];
      }
    }
    /* Compute the output of the network on the perturbed dataset */
    forward(perturbed_input, output_pert, num_hidden, num_pattern,
            num_input, num_output, hidden, hidden_bias, weight_ih, weight_ho);

    for (int i = 0; i < num_pattern; ++i) {
      for (int j = 0; j < num_output; ++j) {
        if (output_pert[i * num_output + j] > 0
            && output[i * num_output + j] > 0) {
          delta = (log(output[i * num_output + j]) -
                  log(output_pert[i * num_output + j]));
          fisher_information += output[i * num_output + j]
            * delta * delta;
        }
      }
    }

  }
  free(output);
  free(hidden);
  free(hidden_bias);
  free(perturbed_input);
  free(output_pert);

  return fisher_information / (AVERAGE_FISHER * num_pattern);
}

double compute_error(int num_pattern, int num_output,
                     int num_hidden, int num_input,
                     double* input, uint8_t* target,
                     double* weight_ih, double* weight_ho)
{
  double test_error = 0.0;
  double val;

  /* Allocate placeholders in the function to make it more
     adaptable to the input dataset. */
  double* output =
    (double*) malloc(sizeof(double) * num_pattern * num_output);
  double* hidden =
    (double*) malloc(sizeof(double) * num_pattern * num_hidden);
  double* hidden_bias =
    (double*) malloc(sizeof(double) * num_pattern * (num_hidden + 1));

  /* Compute the output of the network */
  forward(input, output, num_hidden, num_pattern, num_input, num_output,
          hidden, hidden_bias, weight_ih, weight_ho);

  /* Compute loss */
  for (int p = 0; p < num_pattern; ++p) {
    val = output[p * num_output + target[p]];
    test_error += - log((val > 0) ? val: DBL_MIN);
  }
  test_error /= num_pattern;

  free(output);
  free(hidden);
  free(hidden_bias);

  return test_error;
}

void train_nn_on_automaton(size_t size, int states,
                           uint8_t* train_automaton,
                           uint8_t* test_automaton,
                           network_opts_t* opts,
                           network_result_t* res)
{
  size_t num_pattern = size * size;
  int side = 2 * opts->offset + 1;
  int num_input = states * (side * side - 1);
  int num_hidden = opts->num_hid;
  int num_output = states;

  double batch_error, error, eta = 1, alpha = 0.9;
  if (opts->optim_type != NESTEROV) {
    alpha = 0.;
  }

  int batch_size = 8;
  /* Regression parameter */
  double reg = 0.;

  /* Network and training variables declaration */
  /* num_pattern x (num_input + 1) array that holds all the training set */
  double* base_input =
    (double *) malloc(num_pattern * (num_input + 1) * sizeof(double));
  /* Array that holds the training labels */
  uint8_t* target = (uint8_t *) malloc(num_pattern * sizeof(uint8_t));
  /* Fill those arrays with the automaton's content */
  fill_input_target(size, base_input, target,
                    train_automaton, opts->offset, states);

  double* test_input =
    (double *) malloc(num_pattern * (num_input + 1) * sizeof(double));
  uint8_t* test_target = (uint8_t *) malloc(num_pattern * sizeof(uint8_t));

  /* Weights of the network */
  double weight_ih[(num_input + 1) * num_hidden];
  double weight_ho[(num_hidden + 1) * num_output];

  /* Gradients of the output and hidden layer */
  double delta_output[batch_size * num_output];
  double delta_h[batch_size * num_hidden];

  /* Weight gradients and previous gradient for Nesterov momentum */
  double delta_w_ih[(num_input + 1) * num_hidden];
  double delta_w_ho[(num_input + 1) * num_hidden];
  double* delta_w_ih_prev = NULL;
  double* delta_w_ho_prev = NULL;

  /* Allocate only if being used later */
  if (opts->optim_type == NESTEROV) {
    delta_w_ih_prev = malloc(sizeof(double) * (num_input + 1) * num_hidden);
    delta_w_ho_prev =  malloc(sizeof(double) * (num_hidden + 1) * num_output);
  }

  /* Allocate the arrays that will hold data for each batch */
  double input[(batch_size) * (num_input + 1)];
  double hidden[batch_size * num_hidden];
  double hidden_bias[batch_size * (num_hidden + 1)];
  double output[batch_size * num_output];

  int epoch;
  size_t random_idx[num_pattern];
  double test_error, fish_info = 0.0;

  /* Initialize the weights of the network */
  init_weights(num_hidden, num_input, num_output,
               delta_w_ih, weight_ih,
               delta_w_ho, weight_ho);

  for (epoch = 0; epoch < opts->max_epoch; ++epoch) {
    /* Initialize error */
    error = 0.0;

    /* Learning rate decay */
    if (epoch > 0 && epoch%3 == 0 && opts->decay == DECAY) {
      eta -= .5 * eta;
    }

    /* Random ordering of input patterns done at each epoch */
    shuffle_index(num_pattern, random_idx);

    /* Loop through every batch in the dataset */
    for (size_t s = 0; s < num_pattern; s += batch_size) {
      for (int b = 0; b < batch_size; ++b) {
        if (s + b >= num_pattern) {
          break;
        }
        /* Copy batch elements to the input array for processing */
        memcpy(input + b * (num_input + 1),
               base_input + (random_idx[(s + b)] *
                             (num_input + 1)),
               sizeof(double) * (num_input + 1));
      }

      /* Forward pass */
      forward(input, output, num_hidden, batch_size, num_input, num_output,
              hidden, hidden_bias, weight_ih, weight_ho);

      batch_error = compute_loss(s, batch_size, random_idx, num_output,
                                 output, target);

      /* Compute the gradients for each weight matrix */
      compute_batch_gradients(s, alpha, batch_size, num_input, num_output,
                              num_hidden,
                              random_idx, output, target,
                              delta_output, delta_w_ho, hidden_bias, weight_ho,
                              delta_h, delta_w_ih, input,
                              delta_w_ih_prev, delta_w_ho_prev, opts);

      /* Update the weights of the network */
      update_weights(num_input, num_hidden, num_output, eta, &batch_error, reg,
                     alpha, weight_ih, delta_w_ih, delta_w_ih_prev, weight_ho,
                     delta_w_ho, delta_w_ho_prev);

      /* Error is only updated here because it might be changed by
         `update_weights` */
      error += batch_error * batch_size;
    }
    error /= (double)(num_pattern);

    if (epoch%5 == 0) {
      fprintf(stdout, "\nEpoch %d: Error = %f", epoch, error);
    }
  }

  /* Compute error on the training set */
  error = compute_error(num_pattern, num_output, num_hidden,
                        num_input, base_input, target,
                        weight_ih, weight_ho);

  /* Compute Fisher information if needed by the option */
  if (opts->fisher == FISHER) {
    res->fisher_info = compute_fisher(states, (side*side - 1), num_pattern,
                                      num_output, num_hidden, num_input,
                                      base_input, weight_ih, weight_ho);
  }

  /* Fill the placeholders with test data */
  fill_input_target(size, test_input, test_target,
                    test_automaton, opts->offset, states);
  /* Compute error on the test set */
  test_error = compute_error(num_pattern, num_output, num_hidden,
                             num_input, test_input, test_target,
                             weight_ih, weight_ho);

  /* Log results */
  fprintf(stdout, "\nTrain error: %f\tTest error: %f\tRatio: %f\tDiff: %f\n",
          error, test_error, test_error/error, test_error - error);
  if (opts->fisher== FISHER) {
    fprintf(stdout, "\tFisher: %f\n", fish_info);
  } else {
    fprintf(stdout, "\n");
  }

  /* Output data */
  res->train_error = error;
  res->test_error = test_error;

  /* Cleanup allocated arrays */
  free(base_input);
  free(test_input);
  free(test_target);
  free(delta_w_ih_prev);
  free(delta_w_ho_prev);
}
