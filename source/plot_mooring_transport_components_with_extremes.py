#----------------------------------------------------------------------
# Plots wind and pressure head components of water transport through
# Berig strait as measured at the A3 mooring.  Dataprovided by Rebecca
# Woodgate
#
# 2019-05-01 A.P.Barrett <apbarret@nsidc.org>
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd

from utilities import read_mooring

def main():

    ph = read_mooring('PH', 'MeanCorr').loc['1991':'2016']
    transport = read_mooring('TRANSPORT', 'MeanCorr').loc['1991':'2016']
    wind = transport - ph

    # Merge values to allow sorting
    merge = pd.concat([transport, ph, wind], axis=1,
                      keys=['TRANSPORT', 'PH', 'WIND']).dropna()
    highest = merge.sort_values('TRANSPORT', ascending=False)[0:10]
    lowest = merge.sort_values('TRANSPORT')[0:10]

    print (highest)
    print (lowest)
    
    fig, ax = plt.subplots(3, 1, figsize=(10,6))

    transport.plot(ax=ax[0], label='Total')
    highest['TRANSPORT'].plot(ax=ax[0], marker='.', linestyle='',
                              markersize=10)
    lowest['TRANSPORT'].plot(ax=ax[0], marker='.', linestyle='',
                             markersize=10)
    ax[0].axhline(0., color='k')
    ax[0].set_ylabel('Sv')
    ax[0].set_title('Total Transport')
    ax[0].set_xlim('1991','2017')
    
    ph.plot(ax=ax[1], label='Pressure head')
    highest['PH'].plot(ax=ax[1], marker='.', linestyle='',
                              markersize=10)
    lowest['PH'].plot(ax=ax[1], marker='.', linestyle='',
                             markersize=10)
    ax[1].axhline(0., color='k')
    ax[1].set_ylabel('Sv')
    ax[1].set_title('Pressure head')
    ax[1].set_xlim('1991','2017')
                    
    wind.plot(ax=ax[2], label='Wind')
    highest['WIND'].plot(ax=ax[2], marker='.', linestyle='',
                              markersize=10)
    lowest['WIND'].plot(ax=ax[2], marker='.', linestyle='',
                             markersize=10)
    ax[2].axhline(0., color='k')
    ax[2].set_ylabel('Sv')
    ax[2].set_title('Wind')
    ax[2].set_xlim('1991','2017')

    plt.tight_layout()
    
#    plt.show()
    plt.savefig('a3_mooring_transport_components_with_extremes.png')


if __name__ == "__main__":
    main()
    
    
