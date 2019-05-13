import os

DATA_DIR = 'data'
STEP_DIR = 'steps'

if __name__ == "__main__":
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('_paq.log'):
            out_file = open(os.path.join(DATA_DIR, fname.split('.')[0] + '.dat'), 'w')
            in_file = open(os.path.join(DATA_DIR, fname))
            all_tuple = []
            for line in in_file:
                if '.step' in line and '->' in line:
                    left = line.split(STEP_DIR)[1].split('_')[1].split('.')[0]
                    right = line.split('\x08')[-1].split('->')[1].strip()
                    all_tuple.append((left, right))

            all_tuple.sort(key=lambda x: int(x[0]))
            out_file.write('\n'.join(map('    '.join, all_tuple)))
            out_file.close()
            in_file.close()
