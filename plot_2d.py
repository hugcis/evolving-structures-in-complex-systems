import sys
import numpy as np
import matplotlib.pyplot as plt

num = sys.argv[1]

for p in range(0, 850, 50):
    f = open('step_2d/out{}_{}.step'.format(num, p))
    for i, line in enumerate(f):
        if i == 0:
            a = np.zeros((len(line.strip()), len(line.strip())))
        a[i,:] = list(map(lambda x: 1 if x == '#' or x == '1' else 0, list(line.strip())))
    plt.matshow(-a, cmap='gray')
    plt.axis('off')
    plt.savefig('test{}.png'.format(p))
    plt.close()
