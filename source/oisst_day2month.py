#----------------------------------------------------------------------
# Calculates monthly mean SST from the NOAA OISST v2 data set
#----------------------------------------------------------------------
import xarray as xr
import os
import glob

from utilities import monAnom

DIRI = '/projects/arctic_scientist_data/NOAA_OISSTv2'

def main():

    fileList = glob.glob( os.path.join(DIRI, 'sst', 'daily', 'sst.day.mean.*.v2.nc') )

    for f in fileList:

        print ( 'Processing {:s}'.format(f) )
        
        ds = xr.open_dataset(f)
        dsMon = ds.resample(time='M', label='left').mean(dim='time')

        fout = f.replace('day','month').replace('daily','monthly')
        print ( 'Writing monthly data to {:s}'.format(fout) )
        dsMon.to_netcdf(fout)

    return

if __name__ == "__main__":
    main()

