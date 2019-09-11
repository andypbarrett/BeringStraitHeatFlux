# Extracts U and V winds for Bering Strait from NCEP and MERRA2
import pandas as pd
import xarray as xr
import numpy as np

import datetime as dt

import os
import re
import glob

a3_lat = 66.338
a3_lon = -168.968

ncep_diri='/projects/arctic_scientist_data/Reanalysis/NCEP/hourly'
merra2_diri='/projects/arctic_scientist_data/Reanalysis/MERRA2/hourly'

def ncep_filepath(varname, year):
    """Return NCEP filepath"""
    vardict = {'U10M': 'uwnd', 'V10M': 'vwnd'}
    return os.path.join(ncep_diri, varname, f'{vardict[varname]}.10m.gauss.{year}.nc')


def ncep_filelist(varname, ybeg=1990, yend=2019):
    """Returns a list of ncep filepaths"""
    return sorted([ncep_filepath(varname, y) for y in range(ybeg, yend+1)])


def extract_one_ncep(f, var):
    ds = xr.open_dataset(f)
    return ds[var].sel(lat=a3_lat, lon=360+a3_lon, method='nearest')
    

def extract_ncep_winds(dbeg='19900101', dend='20190228', verbose=False):
    """Extract winds for Bering Strait from NCEP

    Returns: a tuple of xarray DataArrays for uwnd and vwnd
    """

    ybeg=int(dbeg[0:4])
    yend=int(dend[0:4])
    
    if verbose: print ('Extracting U10M for NCEP')
    uwnd = xr.concat([extract_one_ncep(f, 'uwnd') for f in ncep_filelist('U10M', ybeg=ybeg, yend=yend)], dim='time')
    uwnd = uwnd.drop(['lat','lon'])
    uwnd = uwnd.sel(time=slice(dbeg,dend))
    
    if verbose: print ('Extracting V10M for NCEP')
    vwnd = xr.concat([extract_one_ncep(f, 'vwnd') for f in ncep_filelist('V10M', ybeg=ybeg, yend=yend)], dim='time')
    vwnd = vwnd.drop(['lat','lon'])
    vwnd = vwnd.sel(time=slice(dbeg,dend))
    
    return uwnd, vwnd
        

def time_from_filename(f):
    return dt.datetime.strptime(re.search('\d{8}', f).group(), '%Y%m%d')
    

def merra2_filelist(varname, dbeg='19900101', dend='20190228'):
    """Returns a filelist of MERRA2 files

    NB Can't use a generator because filename format is not consistent MERRA2_??? changes
    depending on MERRA2 production run
    """

    dtbeg = dt.datetime.strptime(dbeg, '%Y%m%d')
    dtend = dt.datetime.strptime(dend, '%Y%m%d')
    
    globpath = os.path.join(merra2_diri, varname, '????', '??',
                            f'MERRA2_???.tavg1_2d_slv_Nx.{varname}.????????.nc4')
    files = sorted(glob.glob(globpath))
    return [f for f in files if (time_from_filename(f) >= dtbeg) & (time_from_filename(f) <= dtend)]


def make_time_coords(day, time):
    return day + dt.timedelta(minutes=time.astype(float))


def extract_one_merra2(f, var):
    """Extracts 00, 06, 12, 18 10m wind from one NCEP file"""
    date = dt.datetime.strptime(re.search('\d{8}', f).group(), '%Y%m%d')
    with xr.open_dataset(f) as ds:
        sub = ds[var].sel(time=[0,360,720,1080], lat=a3_lat, lon=a3_lon, method='nearest')
        sub['time'] = [make_time_coords(date, t) for t in sub.time.values]
    return sub


def extract_one_merra2_avg(f, var, verbose=False):
    """Extracts 00, 06, 12, 18 10m wind from one NCEP file"""
    date = dt.datetime.strptime(re.search('\d{8}', f).group(), '%Y%m%d')
    with xr.open_dataset(f) as ds:
        sub = ds[var].sel(time=[0,360,720,1080],
                          lat=slice(64.7602,66.6648),
                          lon=slice(-169.6875, -167.8125)).mean(dim=['lat','lon'])
        sub['time'] = [make_time_coords(date, t) for t in sub.time.values]
    return sub


def extract_merra2_winds(dbeg='19900101', dend='20190228', verbose=False):
    """Extract winds for Bering Strait from NCEP

    Returns: a tuple of xarray DataArrays for uwnd and vwnd
    """
    if verbose: print ('Extracting U10M for MERRA2 for closest cell to A3')
    uwnd = xr.concat([extract_one_merra2(f, 'U10M') for f in merra2_filelist('U10M', dbeg=dbeg, dend=dend)], dim='time')
    uwnd = uwnd.drop(['lat','lon'])
    
    if verbose: print ('Extracting V10M for MERRA2 for closest cell to A3')
    vwnd = xr.concat([extract_one_merra2(f, 'V10M') for f in merra2_filelist('V10M', dbeg=dbeg, dend=dend)], dim='time')
    vwnd = vwnd.drop(['lat','lon'])

    return uwnd, vwnd


def extract_merra2_winds_avg(dbeg='19900101', dend='20190228', verbose=False):
    """Extract winds for Bering Strait from NCEP

    Returns: a tuple of xarray DataArrays for uwnd and vwnd
    """
    if verbose: print ('Extracting U10M for MERRA2 for cells within NCEP cell close to A3')
    uwnd = xr.concat([extract_one_merra2_avg(f, 'U10M') for f in merra2_filelist('U10M', dbeg=dbeg, dend=dend)], dim='time')
    
    if verbose: print ('Extracting V10M for MERRA2 for cells within NCEP cell close to A3')
    vwnd = xr.concat([extract_one_merra2_avg(f, 'V10M') for f in merra2_filelist('V10M', dbeg=dbeg, dend=dend)], dim='time')

    return uwnd, vwnd


def extract_winds(reanalysis, dbeg='19900101', dend='20190228', verbose=False, cell_average=False):

    if reanalysis == 'NCEP':
        uwnd, vwnd = extract_ncep_winds(dbeg=dbeg, dend=dend, verbose=verbose)
        #ds = xr.Dataset({'uwnd': uwnd, 'vwnd': vwnd})

        #filo = f'ncep_10m_wind_a3mooring.{dbeg}to{dend}.nc4'
        #if verbose: print ('Writing to '+filo)
        #ds.to_netcdf(filo, encoding={'uwnd': {'zlib': True, 'complevel': 9},
        #                             'vwnd': {'zlib': True, 'complevel': 9},})
        
    elif reanalysis == 'MERRA2':
        if cell_average:
            uwnd, vwnd = extract_merra2_winds_avg(dbeg=dbeg, dend=dend, verbose=verbose)
        else:
            uwnd, vwnd = extract_merra2_winds(dbeg=dbeg, dend=dend, verbose=verbose)
        
    else:
        print ('Reanalysis not available')
        return
    
    ds = xr.Dataset({'uwnd': uwnd, 'vwnd': vwnd})

    filo = f'{reanalysis}_10m_wind_a3mooring.{dbeg}to{dend}.nc4'
    if cell_average: filo = filo.replace('a3mooring', 'a3mooring.cell_average')
    if verbose: print ('Writing to '+filo)
    ds.to_netcdf(filo, encoding={'uwnd': {'zlib': True, 'complevel': 9},
                                 'vwnd': {'zlib': True, 'complevel': 9},})

    return


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Extracts Bering Strait winds from NCEP and MERRA')
    parser.add_argument('reanalysis', type=str, help='Name of reanalysis - NCEP, MERRA2')
    parser.add_argument('--dbeg', type=str, default='19900101',
                        help='Fist date to extract')
    parser.add_argument('--dend', type=str, default='20190228',
                        help='Last date to extract')
    parser.add_argument('--cell_average', action='store_true',
                        help='Average MERRA2 cells that fall in NCEP cell')
    parser.add_argument('--verbose', '-v', action='store_true')
    
    args = parser.parse_args()
    
    extract_winds(args.reanalysis, dbeg=args.dbeg, dend=args.dend,
                  cell_average=args.cell_average, verbose=True)
