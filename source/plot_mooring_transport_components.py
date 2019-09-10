#----------------------------------------------------------------------
# Plots wind and pressure head components of water transport through
# Berig strait as measured at the A3 mooring.  Dataprovided by Rebecca
# Woodgate
#
# 2019-05-01 A.P.Barrett <apbarret@nsidc.org>
#----------------------------------------------------------------------

import matplotlib.pyplot as plt

from utilities import read_mooring

def main():

    ph = read_mooring('PH', 'MeanCorr')
    transport = read_mooring('TRANSPORT', 'MeanCorr')
    wind = transport - ph

    fig, ax = plt.subplots(figsize=(10,6))

    transport.plot(ax=ax, label='Total')
    ph.plot(ax=ax, label='Pressure head')
    wind.plot(ax=ax, label='Wind')

    ax.axhline(0., color='k')

    ax.legend()
    ax.set_ylabel('Sv')

    plt.show()
#    plt.savefig('a3_mooring_transport_components.png')


if __name__ == "__main__":
    main()
    
    
