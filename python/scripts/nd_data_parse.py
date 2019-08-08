"""
Read 2 states data.
"""
import os
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np

ST = 2

RESULTS_DIC = {}
CPL = []
NN_NAMES = []

def get_results(results_dic, neigh, t_step,
                ratio=False, exp=None, head=None):
    """ Parse data from result dic. """

    static = []
    if ratio:
        fun = lambda x, b: (x[1]/(x[0] if x[0] > 0
                                  else x[0] + 1e-12)) if not b else 1
    else:
        fun = lambda x, b: x[1] - x[0] if not b else 0

    if exp == 'prod':
        fun = lambda x, b: (x[1] * x[0]) if not b else 0
    elif exp == 'sum':
        fun = lambda x, b: -(x[1] + x[0]) if not b else -5000
    elif exp == 'test':
        fun = lambda x, b: (-x[1]) if not b else -100
    elif exp == 'train':
        fun = lambda x, b: (-x[0]) if not b else -100
    elif exp == "harmonic":
        fun = lambda x, b: -(x[1] + x[0]) if not b else 0

    if head is not None:
        return ([rule for rule in head if rule in results_dic],
                np.array([fun(results_dic[rule][neigh][t_step],
                              rule in static)
                          for rule in head if rule in results_dic]))

    return ([rule for rule in results_dic],
            np.array([fun(results_dic[rule][neigh][t_step], rule in static)
                      for rule in results_dic]))

for i in [t for t in os.listdir('data_2d_{}/nn/'.format(ST))
          if t.endswith('dat')]:
    with open(os.path.join('data_2d_{}/nn/'.format(ST), i)) as f:
        x = f.read().strip()
        if not x:
            continue
        with open(os.path.join('data_2d_{}/out/'.format(ST),
                               i.replace('nn', 'out'))) as f:

            s = f.read().split('\n')
            try:
                CPL.append(int(s[-2].split('    ')[1]))
            except IndexError:
                continue
            except ValueError:
                continue

        name = i.split('.')[0].strip('nn')
        NN_NAMES.append(name)
        v = list(map(lambda x: [float(i) for i in x.split('    ')], x.split('\n')))

        RESULTS_DIC[name] = {}

        for radius in range(1, 8):
            RESULTS_DIC[name][radius] = {
                300: [v[radius - 1][0], v[radius - 1][1]]
            }

pkl.dump(RESULTS_DIC, open('python/data/2d_dataset.pkl', 'wb+'))
print(" ".join(NN_NAMES))
plt.plot(get_results(RESULTS_DIC, 3, 300, exp="sum")[1])
plt.show()
