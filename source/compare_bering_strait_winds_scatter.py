# Compare wind from MERRA2 and NCEP for the Bering Strait A3 mooring

from extract_winds import extract_merra2_winds, extract_merra2_winds_avg, extract_ncep_winds
import matplotlib.pyplot as plt

import xarray as xr
import pandas as pd
import numpy as np

#merra2_filepath = 'MERRA2_10m_wind_a3mooring.19960101to20141231.nc4'
#merra2_cellavg_filepath = 'MERRA2_10m_wind_a3mooring.cell_average.19960101to20141231.nc4'
#ncep_filepath = 'NCEP_10m_wind_a3mooring.19960101to20190228.nc4'
merra2_filepath = 'MERRA2_10m_wind_a3mooring.19900101to20190228.nc4'
merra2_cellavg_filepath = 'MERRA2_10m_wind_a3mooring.cell_average.19900101to20190228.nc4'
ncep_filepath = 'NCEP_10m_wind_a3mooring.19900101to20190228.nc4'

def load_one(f):
    """Loads a single file and returns a pandas dataframe"""
    return xr.open_dataset(f).to_dataframe()

def get_data():
    """Loads MERRA2 and NCEP data into Pandas dataframe"""
    df = pd.concat([load_one(f) for f in [merra2_filepath, merra2_cellavg_filepath, ncep_filepath]],
                   axis=1, keys=['MERRA2', 'MERRA2_CELLAVG', 'NCEP'])
    return df

def rmse(x, y):
    df = pd.concat([x, y], axis=1, keys=['x','y']).dropna()
    df['error2'] = (df.x - df.y)**2
    return np.sqrt(df.error2.mean())

def bias(x, y):
    return x.mean() - y.mean()

def correlation(x, y):
    df = pd.concat([x, y], axis=1, keys=['x','y']).dropna()
    return df.corr()['x']['y']

def plot_panel(x, y, ax=None, xlabel='', ylabel='', vmin=-20, vmax=20, title=None):
    
    if not ax:
        ax = plt.subplot(111)
        
    ax.plot(x, y, '.')
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.plot([vmin,vmax], [vmin,vmax], '-', color='0.3')
    ax.set_aspect('equal')
    ax.set_title(title)
    
    ax.text(0.02, 0.92, f'RMSE: {rmse(x, y):5.2f}', transform=ax.transAxes)
    ax.text(0.02, 0.86, f'Correlation: {correlation(x, y):5.2f}', transform=ax.transAxes)
    ax.text(0.02, 0.80, f'Bias: {bias(x, y):5.2f}', transform=ax.transAxes)
    
    return ax

def get_limits(df):
    xmin = np.round(df.min().min()/5.)*5.
    xmax = np.round(df.max().max()/5.)*5.
    return xmin, xmax

def main():

    df = get_data()

    vmin, vmax = get_limits(df)
    
    fig, ax = plt.subplots(2, 2, figsize=(10, 6))

    plot_panel(df['NCEP']['uwnd'], df['MERRA2']['uwnd'], ax[0,0],
               xlabel='U10M NCEP', ylabel='U10M MERRA2 - nearest cell',
               vmin=vmin, vmax=vmax)
    plot_panel(df['NCEP']['uwnd'], df['MERRA2_CELLAVG']['uwnd'], ax[0,1],
               xlabel='U10M NCEP', ylabel='U10M MERRA2 - NCEP cell',
               vmin=vmin, vmax=vmax)
    
    plot_panel(df['NCEP']['vwnd'], df['MERRA2']['vwnd'], ax[1,0],
               xlabel='V10M NCEP', ylabel='V10M MERRA2 - nearest cell',
               vmin=vmin, vmax=vmax)
    plot_panel(df['NCEP']['vwnd'], df['MERRA2_CELLAVG']['vwnd'], ax[1,1],
               xlabel='V10M NCEP', ylabel='V10M MERRA2 - NCEP cell',
               vmin=vmin, vmax=vmax)
    
    plt.tight_layout()

    #fig.savefig('compare_bering_strait_winds_1996_scatter.png')

    plt.show()

if __name__ == "__main__":
    main()
    
