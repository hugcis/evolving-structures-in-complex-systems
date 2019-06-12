import os
import sys
import numpy as np
import matplotlib.pyplot as plt

STEP_DIR = "steps/"

RULE = sys.argv[1]

with open(os.path.join(STEP_DIR, "out{}.steps".format(RULE))) as f:
    arr = np.array([[int(i) for i in list(t)] for t in f.read().strip().split('\n')])
plt.matshow(-arr, cmap='gray')
plt.axis('off')
plt.show()
