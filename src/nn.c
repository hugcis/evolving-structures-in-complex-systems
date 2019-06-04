#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "nn.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#define NUMHID 10
#define MAX_EPOCH 10

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

double rand_normal()
{
  return randn(0.0, 1.0);
}

void fill_input_target(size_t size, double* input, uint8_t* target,
                       uint8_t automaton[size][size], int offset)
{
  size_t index;
  int counter;
  int side = 2 * offset + 1;
  int num_input = side * side - 1;

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      index = i * size + j;
      counter = 1;

      /* Add bias in the main vector */
      input[index * (num_input + 1)] = 1.0;

      for (int a = -offset; a < offset + 1; ++a) {
        for (int b = -offset; b < offset + 1; ++b) {
          if (a != 0 || b != 0) {  /* Don't take index i,j */
            input[index * (num_input + 1) + counter] =
              (double)automaton[(i + a + size) % size][(j + b + size) % size];
            counter++;
          }
        }
      }

      target[index] = automaton[i][j];
    }
  }
}

void init_weights(int num_hidden, int num_input, int num_output,
                  double* delta_w_ih,
                  double* weight_ih,
                  double* delta_w_ho,
                  double* weight_ho)
{
  int i, j, k;
  /* Initialize weights input -> hidden */
  for (j = 0; j < num_hidden; ++j) {
    delta_w_ih[j] = 0.0;
    weight_ih[j] = 0.0;
    for (i = 1; i < num_input + 1; ++i) {
      delta_w_ih[i * num_hidden + j] = 0.0;
      weight_ih[i * num_hidden + j] = rand_normal() *
        sqrt(1 / (num_input + num_hidden));
    }
  }

  /* Initialize weights hidden -> output */
  for (k = 0; k < num_output; ++k) {
    delta_w_ho[k] = 0.0;
    weight_ho[k] = 0.0;
    for (j = 1; j < num_hidden + 1; ++j) {
      delta_w_ho[j * num_output + k] = 0.0;
      weight_ho[j * num_output + k] = rand_normal() *
        sqrt(1 / (num_output + num_hidden));
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

void write_results_to_file(FILE* file, double train_error, double test_error)
{
  fprintf(file, "%f    %f\n", train_error, test_error);
}

void train_nn_on_automaton(size_t size, int states,
                           uint8_t train_automaton[size][size],
                           uint8_t test_automaton[size][size], int offset,
                           FILE* output_file)
{
  size_t num_pattern = size * size;
  int side = 2 * offset + 1;
  int num_input = side * side - 1;
  int num_hidden = NUMHID;
  int num_output = states;

  /* Network and training variables declaration */
  double* base_input =
    (double *) malloc(num_pattern * (num_input + 1) * sizeof(double));
  uint8_t* target = (uint8_t *) malloc(num_pattern * sizeof(uint8_t));
  fill_input_target(size, base_input, target, train_automaton, offset);

  double weight_ih[(num_input + 1) * num_hidden];

  double weight_ho[(num_hidden + 1) * num_output];

  double delta_output[num_output];
  double delta_h[num_hidden];

  double* delta_w_ih =
    malloc(sizeof(double) * (num_input + 1) * num_hidden);
  double* delta_w_ho =
    malloc(sizeof(double) * (num_hidden + 1) * num_output);
  double batch_error, error, eta = 0.0000001, alpha = 0.02, batch_size = 8;

  double* input =
    (double*) malloc (sizeof(double) * (batch_size) * (num_input + 1));
  double* hidden =
    malloc(sizeof(double) * batch_size * num_hidden);
  double* hidden_bias =
    malloc(sizeof(double) * batch_size * (num_hidden + 1));
  double* output =
    malloc(sizeof(double) * batch_size * num_output);

  int i, j, k, epoch;
  size_t np, op, p;
  size_t random_idx[num_pattern];

  init_weights(num_hidden, num_input, num_output,
               delta_w_ih, weight_ih,
               delta_w_ho, weight_ho);

  for (epoch = 0; epoch < MAX_EPOCH; ++epoch) {
    error = 0.0;
    /* Random ordering of input patterns */
    for (p = 0; p < num_pattern; ++p) {
      random_idx[p] = p;
    }

    for (p = 0; p < num_pattern; ++p) {
      np = p + ((double)rand()/((double)RAND_MAX+1)) * (num_pattern - 1 - p);
      op = random_idx[p];
      random_idx[p] = random_idx[np];
      random_idx[np] = op;
    }

    batch_error = 0.0;
    for (size_t s = 0; s < num_pattern / batch_size; s += batch_size) {
      for (int b = 0; b < batch_size; ++b) {
        if (s + b >= num_pattern) {
          break;
        }
        memcpy(input + b * (num_input + 1),
               &base_input[random_idx[(s + b)] * (num_input + 1)],
               sizeof(double) * (num_input + 1));
      }
      /* Forward pass */
      forward(input, output, num_hidden, batch_size, num_input, num_output,
              hidden, hidden_bias, weight_ih, weight_ho);

      for (int b = 0; b < batch_size; ++b) {
        p = random_idx[(s + b)];
        /* Compute loss */
        batch_error +=
          - log((output[b * num_output + target[p]] > 0) ?
                output[b * num_output + target[p]]: 1e-15);

        /* Backpropagation */
        for (k = 0; k < num_output; ++k) {
          /* Output gradients */
          delta_output[k] = output[b * num_output + k] -
            ((k == target[p])? 1.0: 0.0);

          for (j = 0; j < num_hidden + 1; ++j) {
            delta_w_ho[j * num_output + k] +=
              (eta * hidden_bias[b * (num_hidden + 1) + j] *
              delta_output[k] + alpha * delta_w_ho[j * num_output + k]) /
              batch_size;
          }
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, num_hidden, num_output,
                    1.0, (weight_ho + num_output), num_hidden,
                    delta_output, 1, 0.0, delta_h, 1);
        /* Hidden layer gradients */
        for (j = 0; j < num_hidden; ++j) {
          delta_h[j] *= (hidden_bias[b * (num_hidden + 1) + j + 1] > 0) ?
            1.0: 0.0;

          for (i = 0; i < num_input + 1; ++i) {
            delta_w_ih[i * num_hidden + j] +=
              (eta * input[b * (num_input + 1) + i] *
              delta_h[j] + alpha * delta_w_ih[i * num_hidden + j]) / batch_size;
          }
        }
      }
      batch_error /= batch_size;
      error += batch_error;

      /* Update weights with average gradient */
      for (i = 0; i < num_input + 1; ++i) {
        for (j = 0; j < num_hidden; ++j) {
          weight_ih[i * num_hidden + j] -= delta_w_ih[i * num_hidden + j];
        }
      }
      for (j = 0; j < num_hidden + 1; ++j) {
        for (k = 0; k < num_output; ++k) {
          weight_ho[j * num_output + k] -= delta_w_ho[j * num_output + k];
        }
      }
    }
    error /= (double)(num_pattern/batch_size);

    if (epoch%5 == 0) fprintf(stderr, "\nEpoch %d: Error = %f", epoch, error);
    /* if( error < 0.0004 ) break ;  /\* stop learning when 'near enough' *\/ */
  }

  fill_input_target(size, base_input, target, test_automaton, offset);

  output =
    (double *) realloc(output, sizeof(double) * num_pattern * num_output);
  hidden =
    (double *) realloc(hidden, sizeof(double) * num_pattern * num_hidden);
  hidden_bias =
    (double *) realloc(hidden_bias,
                       sizeof(double) * num_pattern * (num_hidden + 1));
  forward(base_input, output, num_hidden, num_pattern, num_input, num_output,
          hidden, hidden_bias, weight_ih, weight_ho);

  double test_error = 0.0;
  for (p = 0; p < num_pattern; ++p) {
    test_error += - log((output[p * num_output + target[p]] > 0) ?
                        output[p * num_output + target[p]]: 1e-15);
  }
  test_error /= num_pattern;

  fprintf(stderr, "\nTrain error: %f\tTest error: %f\n", error, test_error);
  write_results_to_file(output_file, error, test_error);

  free(base_input);
  free(input);
  free(target);
  free(hidden);
  free(hidden_bias);
  free(output);
  free(delta_w_ih);
  free(delta_w_ho);
}

