# Compare wind from MERRA2 and NCEP for the Bering Strait A3 mooring

from extract_winds import extract_merra2_winds, extract_merra2_winds_avg, extract_ncep_winds
import matplotlib.pyplot as plt

date_begin = '19960101'
date_end = '19961231'

def main():

    ncep_uwnd, ncep_vwnd = extract_ncep_winds(date_begin, date_end, verbose=True)

    merra2_uwnd, merra2_vwnd = extract_merra2_winds(date_begin, date_end, verbose=True)

    merra2_uwnd_avg, merra2_vwnd_avg = extract_merra2_winds_avg(date_begin, date_end, verbose=True)

    fig, ax = plt.subplots(2, 1, figsize=(10, 8))

    ncep_uwnd.plot(ax=ax[0], label='NCEP - nearest cell')
    merra2_uwnd.plot(ax=ax[0], label='MERRA2 - nearest cell')
    merra2_uwnd_avg.plot(ax=ax[0], label='MERRA2 - NCEP cell')
    ax[0].legend()

    ncep_vwnd.plot(ax=ax[1], label='NCEP - nearest cell')
    merra2_vwnd.plot(ax=ax[1], label='MERRA2 - nearest cell')
    merra2_vwnd_avg.plot(ax=ax[1], label='MERRA2 - NCEP cell')
    ax[0].legend()

    plt.tight_layout()

    fig.savefig('compare_bering_strait_winds_1996_timeseries.png')

    plt.show()

if __name__ == "__main__":
    main()
    
