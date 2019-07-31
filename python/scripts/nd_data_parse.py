import os
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np

ST = 2

results_dic = {}
cpl = []
nn_names = []

def get_results(results_dic, neigh, ts, ratio=False, exp=None, head=None):
    static = []
    if ratio:
        fun = lambda x, b: (x[1]/(x[0] if x[0] > 0 else x[0] + 1e-12)) if not b else 1
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
                np.array([fun(results_dic[rule][neigh][ts], rule in static)
                          for rule in head if rule in results_dic]))

    return ([rule for rule in results_dic],
            np.array([fun(results_dic[rule][neigh][ts], rule in static)
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
                cpl.append(int(s[-2].split('    ')[1]))
            except Exception:
                continue
        name = i.split('.')[0].strip('nn')
        nn_names.append(name)
        v = list(map(lambda x: [float(i) for i in x.split('    ')], x.split('\n')))

        results_dic[name] = {}

        for radius in range(1, 8):
            results_dic[name][radius] = {
                300: [v[radius - 1][0], v[radius - 1][1]]
            }

pkl.dump(results_dic, open('python/data/2d_dataset.pkl', 'wb+'))
print(" ".join(nn_names))
plt.plot(get_results(results_dic, 3, 300, exp="sum")[1])
plt.show()
