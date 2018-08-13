import matplotlib.pyplot as plt
import numpy as np

from utilities import read_mooring

def plot_panel(x, ax, ylabel='', q=5):

    xsrt = x.dropna().sort_values(ascending=False)
    n=xsrt.size
    p=np.ceil(n*q/100.).astype(int)
    xsrt = xsrt[:p]
    
    ax.plot(x.index.month, x, 'o')
    ax.plot(xsrt.index.month, xsrt, 'o', c='r')
    ax.set_ylabel(ylabel)
    
    return

def main():

    heat = read_mooring(type='HEAT', column='MeanCorr')
    transport = read_mooring(type='TRANSPORT', column='MeanCorr')
    temperature = read_mooring(type='TEMPERATURE', column='Mean')

    fig, ax = plt.subplots(3,1, figsize=(7,15))

    plot_panel(heat, ax[0], ylabel='Heat Flux')
    plot_panel(transport, ax[1], ylabel='Transport')
    plot_panel(temperature, ax[2], ylabel='Temperature')

    plt.show()
    
    return

if __name__ == "__main__":
    main()

