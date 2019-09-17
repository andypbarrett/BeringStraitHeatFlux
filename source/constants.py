import socket
import cartopy.crs as ccrs

mooring_dirpath = {
    'nsidc-abarrett-442': '/home/apbarret/Data/BeringStraitHeatFlux/Mooring',
    'nsidc-abarrett': '/home/apbarret/data/BeringStraitHeatFlux/Mooring',
}
MOORING_DIRPATH=mooring_dirpath[socket.gethostname()]

reanalysis_dirpath = {
    'nsidc-abarrett-442': '/home/apbarret/Data/BeringStraitHeatFlux/MERRA2',
    'nsidc-abarrett': '/disks/arctic5_raid/abarrett/MERRA2/monthly', 
}
REANALYSIS_DIRPATH = reanalysis_dirpath[socket.gethostname()]

MOORING_FILEPATH = {
    'HEAT': 'BeringStrait_Monthlymeans_HEAT_Oct2017.txt',
    'PH': 'BeringStrait_Monthlymeans_PH_April2019.txt',
    'TEMPERATURE': 'BeringStrait_Monthlymeans_TEMPERATURE_Oct2017.txt',
    'TRANSPORT': 'BeringStrait_Monthlymeans_TRANSPORT_Oct2017.txt',
    'WIND': 'BeringStrait_MonthlyMeans_WIND_May2019.txt',
    }

REANALYSIS_FILEPATH = {
    'SLP': '####PLACE_HOLDER####',
    'U10M': '####ANOTHER_PLACE_HOLDER####'
}

MAP_PROJ = ccrs.NorthPolarStereo()  # Map projection for plots
MAP_EXTENT = [-4651194.319071749, 2835721.446215284, -2835721.446215284, 4651194.319071749]  # Map extent in NorthPolarStereo



