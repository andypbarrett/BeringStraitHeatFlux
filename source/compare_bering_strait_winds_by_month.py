import matplotlib.pyplot as plt
import calendar

import xarray as xr
import pandas as pd
import numpy as np

from affine import Affine

from compare_bering_strait_winds_scatter import get_data, get_limits, plot_panel

merra2_filepath = 'MERRA2_10m_wind_a3mooring.19900101to20190228.nc4'
merra2_cellavg_filepath = 'MERRA2_10m_wind_a3mooring.cell_average.19900101to20190228.nc4'
ncep_filepath = 'NCEP_10m_wind_a3mooring.19900101to20190228.nc4'

winddct = {'uwnd': 'U10M', 'vwnd': 'V10M', 'rotated': 'V(330) 10M'}
reandct = {'NCEP': 'NCEP',
           'MERRA2': 'MERRA2 - nearest cell',
           'MERRA2_CELLAVG': 'MERRA2 - cell average'}

def main(xreanalysis, yreanalysis, wind):

    df = get_data()

    if wind is 'rotated':
        for col in df.columns.levels[0]:
            urot, vrot = (df[col].uwnd, df[col].vwnd) * Affine.rotation(30)
            df.loc[:, (col,'rotated')] = vrot

    vmin, vmax = get_limits(df)

    fig, ax = plt.subplots(3, 4, figsize=(20,15))

    ax = ax.flatten()

    xlabel = f'{winddct[wind]} {reandct[xreanalysis]}'
    ylabel = f'{winddct[wind]} {reandct[yreanalysis]}'
    
    for i, (name, group) in enumerate(df.groupby(df.index.month)):

        plot_panel(group[xreanalysis][wind], group[yreanalysis][wind], ax[i],
                   xlabel=xlabel, ylabel=ylabel,
                   vmin=vmin, vmax=vmax, title=calendar.month_name[name])

    plt.subplots_adjust(wspace=0.3, hspace=0.35)

    #plt.tight_layout()
    plt.show()
    #plt.savefig(f'{xreanalysis}_vs_{yreanalysis}_{wind}.png')
    
if __name__ == "__main__":
    main('NCEP', 'MERRA2', 'rotated')
    
