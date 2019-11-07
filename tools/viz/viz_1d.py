import os
import sys
import numpy as np
import matplotlib.pyplot as plt

STEP_DIR = "steps/"

RULE = sys.argv[1:]

COLS = 3 if len(RULE) > 3 else len(RULE)
ROWS = len(RULE)//3 + (1 if not len(RULE)//3 else 0)
fig, axs = plt.subplots(ROWS,
                        COLS, sharex=True, sharey=True,
                        figsize=(6 * COLS, 4 * ROWS))
if len(RULE) == 1:
    axs = [[axs]]
elif len(RULE) <= 3:
    axs = [axs]

for n, rule in enumerate(RULE):
    ax = axs[n//3][n%3]
    with open(os.path.join(STEP_DIR, "out{}.steps".format(rule))) as f:
        arr = np.array([[int(i) for i in list(t)]
                        for t in f.read().strip().split('\n')])
    ax.matshow(-arr, cmap='gray')
    ax.set_xlim(arr.shape[1]/2 - 30.5, arr.shape[1]/2 + 30.5)
    ax.set_ylim(30.5, -.5)
    ax.set_title("Rule {}".format(rule))
    ax.axis('off')

fig.tight_layout()
fig.savefig('figures/output_viz.pdf')
plt.show()
