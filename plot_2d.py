import sys
import os
import imageio
import numpy as np
import matplotlib.pyplot as plt

NUM = sys.argv[1]
STATE = sys.argv[2]
GIF = False
if len(sys.argv) >= 4 and sys.argv[3] == '-g':
    GIF = True

for p in range(0, 850, 50):
    filename = 'step_2d_{}/out{}_{}.step'.format(STATE, NUM, p)
    if not filename.split('/')[1] in os.listdir('step_2d_{}'.format(STATE)):
        print("Filename {} not found".format(filename))
        continue


    f = open(filename)
    for i, line in enumerate(f):
        if i == 0:
            a = np.zeros((len(line.strip()), len(line.strip())))
        a[i, :] = list(map(int,
                           list(line.strip())))
    plt.matshow(-a, cmap='gray')
    plt.axis('off')
    plt.savefig('rule_gif/tmp{}.png'.format(p))
    plt.close()

if GIF:
    images = []
    for p in range(0, 850, 50):
        filename = 'rule_gif/tmp{}.png'.format(p)
        if not filename.split('/')[1] in os.listdir('rule_gif'):
            continue
        images.append(imageio.imread(filename))
    imageio.mimsave('rule_gif/{}.gif'.format(NUM), images, duration=0.2)
