import numpy as np
import pylab as plt
import os

def plot_data(directory):
    print('Plotting data...')
    for file in os.listdir(directory):
        if(not file.endswith('.txt')):
            continue
        fig = plt.figure()
        print(file)
        values = np.loadtxt(f'{directory}/{file}')
        plt.rcParams['lines.linewidth'] = 1
        plt.rcParams['lines.color'] = 'b'
        plt.rcParams['lines.markersize'] = 0

        plt.plot(values[:,0], values[:,1], '.-')
        plt.xlabel('z')
        plt.ylabel('Intensity')
        plt.title('Intensity vs z, oil thickness = ' + file[5:-4])
        plt.savefig(f'plots/plot_{file[5:-4]}.png')

if __name__ == '__main__':
    plot_data('')