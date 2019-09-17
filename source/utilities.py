import os
import pandas as pd
import xarray as xr

from constants import MOORING_DIRPATH, MOORING_FILEPATH, REANALYSIS_DIRPATH

def read_mooring(var='HEAT', column=None):
    """
    Reads Heat flux, transport or temperature data derived from Bering Strait
    A3 mooring data.

    Arguments: type - type of data to read
                      HEAT - heatflux
                      TRANSPORT - water tranport
                      TEMPERATURE - water temperature

    Returns: pandas dataframe
    """

    names = ['Mooring', 'Year', 'Month',
             'Mean', 'Error',
             'MeanCorr', 'CorrErr']
    
    filepath = os.path.join( MOORING_DIRPATH,
                             'Month',
                             MOORING_FILEPATH[var] )
    if var == 'WIND':
        df = pd.read_csv(filepath, header=0, index_col=0, parse_dates=True)
    else:
        df = pd.read_csv(filepath, header=None, comment='%',sep='\s+', 
                         names=names,
                         parse_dates={'date': ['Year', 'Month']},
                         index_col='date')

    if column:
        return df[column]
    else:
        return df

def get_composite_dates(s, N=5, highest=True):
    """
    Returns dates of N highest or lowest values in series s.
    """
    if highest:
        ascending=False
    else:
        ascending=True
    return s.dropna().sort_values(ascending=ascending)[:N].index

def MonClim(ds):
    """
    Calculates monthly climatolgy from monthly data

    :param ds: xarray DataSet
    :return: xarray DataSet of monthly climatology
    """
    return ds.dropna(dim='time').groupby('time.month').mean(dim='time')

def get_reanalysis(variable='SLP'):
    """
    Returns time series of monthly reanalysis grids for a given variable
    (SLP, U10M, V10M).
    """
    fileGlob = os.path.join(REANALYSIS_DIRPATH, variable, 
                            'MERRA2_???.tavgM_2d_slv_Nx.{:s}.??????.nc4'.format(variable))
    ds = xr.open_mfdataset(fileGlob)
    da = ds[variable].load()
    da = da * 1e-2 # Convert to hPa
    da.attrs['units'] = "hPa"
    return da

def lagDate(dates, lag):
    """
    Gets dates for month lag: e.g. if Month is 6 and lag=-1, returns month 5
    """
    if lag > 12:
        print ("Doesn't handle lags greater than 12 months")
        return
    
    if lag == 0: return dates

    newdates = []
    for y, m in zip(dates.year, dates.month):
        if (m == 1) & (lag < 0):
            m = 13+lag
            y = y-1
        else:
            m = m - 1
        newdates.append(dt.datetime(y,m,1))
    return newdates
    
def add_construction_grid(fig, hori=None, vert=None, novert=False, nohori=False):

    if not hori: hori=np.linspace(0,1,11)
    if not vert: vert=np.linspace(0,1,11)

    print (vert)
    
    ll = []
    if not nohori:
        for h in hori:
            ll.append( mpl.lines.Line2D([0,1],[h,h], transform=fig.transFigure, figure=fig, ls='--', lw=1) )

    if not novert:
        for v in vert:
            ll.append( mpl.lines.Line2D([v,v],[0,1], transform=fig.transFigure, figure=fig, ls='--', lw=1) )

    fig.lines.extend(ll)
    
    return

def monAnom(x, year_begin='1981', year_end='2010'):
    """
    Calculates anomalies of x relative to the 1981-2010 period.

    N.B. Ideally, I would like to use an arbitrary mean field
    """
                                     
    def _anom(x):
        return x - x.sel(time=slice(year_begin,year_end)).mean(dim='time')

    return x.groupby('time.month').apply(_anom)
                                     
