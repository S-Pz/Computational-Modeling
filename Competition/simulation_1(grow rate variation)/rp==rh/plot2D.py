import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys 
import os
import glob
import pandas as pd

def getMaxMinValueFrom2DMatrix(u):
    maxv = np.max(np.max(u))
    minv = np.min(np.min(u))
    return [maxv,minv]

data_path = 'results/'
fig_path = 'figs/'

x = np.loadtxt(data_path + 'x.csv')
y = np.loadtxt(data_path + 'y.csv')
times = np.loadtxt(data_path + 't.csv')

values = [("P", "blue"), ("H", "red")]

for name, color in values: 
    file_name = name + "*.csv"
    all_files = glob.glob(os.path.join(data_path, file_name))
    
    for f in all_files:
        t = f.split("_")[1].split(".csv")[0]
        table = pd.read_csv(f)
        fig_name = f.split("_")[0].split("/")[1]

        Z = np.array(table["value"])
        Z = np.reshape(Z,(np.shape(x)[0],np.shape(y)[0]))
        
        cm = plt.get_cmap('jet')
        fig = plt.figure(figsize=(10,7))
        plt.title(fig_name)
        pcm = plt.pcolormesh(x, y, Z, cmap=cm) #, vmin=minDensity, vmax=maxDensity)
        plt.colorbar(pcm)        
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        fig.savefig(fig_path + fig_name + '_' + t + '.png', format='png', bbox_inches='tight')
        plt.close()