import socket

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



