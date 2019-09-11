import matplotlib.pyplot as plt
import calendar

import xarray as xr
import pandas as pd
import numpy as np

from affine import Affine

from compare_bering_strait_winds_scatter import get_data

def main():

    df = get_data()

    vwnd_diff = df.NCEP.vwnd - df.MERRA2.vwnd
    uwnd_diff = df.NCEP.uwnd - df.MERRA2.uwnd
    vwnd_avg_diff = df.NCEP.vwnd - df.MERRA2_CELLAVG.vwnd
    uwnd_avg_diff = df.NCEP.uwnd - df.MERRA2_CELLAVG.uwnd

    # Get rotated winds - rotations are counter clockwise 330 N == 30
    MERRA2_rot = pd.concat((df.MERRA2.uwnd, df.MERRA2.vwnd) * Affine.rotation(30),
                           axis=1, keys=['urot','vrot'])
    MERRA2_CELLAVG_rot = pd.concat((df.MERRA2_CELLAVG.uwnd, df.MERRA2_CELLAVG.vwnd) * Affine.rotation(30),
                                   axis=1, keys=['urot','vrot'])
    NCEP_rot = pd.concat((df.NCEP.uwnd, df.NCEP.vwnd) * Affine.rotation(30),
                         axis=1, keys=['urot','vrot'])
    vrot_diff = NCEP_rot.vrot - MERRA2_rot.vrot
    vrot_avg_diff = NCEP_rot.vrot - MERRA2_CELLAVG_rot.vrot
    
    fig, ax = plt.subplots(2, 1, figsize=(10,15))

    vwnd_diff.rolling(360, win_type='boxcar', center=True).mean().plot(ax=ax[0], label='v')
    uwnd_diff.rolling(360, win_type='boxcar', center=True).mean().plot(ax=ax[0], label='u')
    vrot_diff.rolling(360, win_type='boxcar', center=True).mean().plot(ax=ax[0], label='330 deg.')
    ax[0].axhline(0.,ls='--', c='k')
    ax[0].set_ylabel('m/s')
    ax[0].set_title('NCEP - MERRA2')
    ax[0].legend()

    vwnd_avg_diff.rolling(360, win_type='boxcar', center=True).mean().plot(ax=ax[1], label='v')
    uwnd_avg_diff.rolling(360, win_type='boxcar', center=True).mean().plot(ax=ax[1], label='u')
    vrot_avg_diff.rolling(360, win_type='boxcar', center=True).mean().plot(ax=ax[1], label='330 deg.')
    ax[1].axhline(0.,ls='--', c='k')
    ax[1].set_ylabel('m/s')
    ax[1].set_title('NCEP - MERRA2_CELLAVG')

    plt.show()

if __name__ == "__main__":
    main()
    

