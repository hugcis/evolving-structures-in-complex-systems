import os

if __name__ == "__main__":
    for fname in os.listdir('data_2d'):
        if fname.endswith('_paq.log'):
            out_file = open(os.path.join('data_2d', fname.split('.')[0] + '.dat'), 'w')
            in_file = open(os.path.join('data_2d', fname))
            all_tuple = []
            for line in in_file:
                if '.step' in line and '->' in line:
                    left = line.split('_')[2].split('.')[0]
                    right = line.split('\x08')[-1].strip()
                    all_tuple.append((left, right))

            all_tuple.sort(key=lambda x: int(x[0]))
            out_file.write('\n'.join(map('    '.join, all_tuple)))
            out_file.close()
            in_file.close()
