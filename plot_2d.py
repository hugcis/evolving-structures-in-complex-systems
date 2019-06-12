import sys
import os
import imageio
from python.nb.utils import plot_from_file

NUM = sys.argv[1]
STATE = sys.argv[2]
GIF = False
CROP = False
if len(sys.argv) >= 4 and sys.argv[3] == '-g':
    GIF = True
if len(sys.argv) >= 5 and sys.argv[4] == '-c':
    CROP = True

ALL_NAMES = [name.split('out{}_'.format(NUM))[1].strip('.step')
             for name in os.listdir('step_2d_{}'.format(STATE))
             if 'out{}'.format(NUM) in name]
STEPS = list(map(int, ALL_NAMES))
STEPS = sorted(STEPS)
N_STEPS = max(STEPS)
DIFF = STEPS[-1] - STEPS[-2]

for p in range(0, N_STEPS + 1, DIFF):
    filename = 'step_2d_{}/out{}_{}.step'.format(STATE, NUM, p)
    plot_from_file(filename, crop=CROP,
                   save_as='rule_gif/tmp{}.png'.format(p))

images = []

for t in range(0, N_STEPS + 1, DIFF):
    filename = 'rule_gif/tmp{}.png'.format(t)
    if not filename.split('/')[1] in os.listdir('rule_gif'):
        continue
    images.append(imageio.imread(filename))

if GIF:
    imageio.mimsave('rule_gif/{}.gif'.format(NUM),
                    images, duration=0.2)
else:
    imageio.mimsave('rule_gif/tmp.gif', images, duration=0.2)
