# Find the ranks of the transport components relatig to the 10 highest and
# ten lowest transports

import pandas as pd

from utilities import read_mooring

def main():

    ph = read_mooring('PH', 'MeanCorr').loc['1991':'2016']
    transport = read_mooring('TRANSPORT', 'MeanCorr').loc['1991':'2016']
    wind = transport - ph

    # Merge values to allow sorting
    merge = pd.concat([transport, ph, wind], axis=1,
                      keys=['TRANSPORT', 'PH', 'WIND']).dropna()
    merge = merge.join(merge.rank(pct=True), rsuffix='_Rank')
    merge = merge.round(2)
    
#    highest = merge.sort_values('TRANSPORT', ascending=False)[0:10]
#    lowest = merge.sort_values('TRANSPORT')[0:10]

    merge = merge.sort_values('TRANSPORT', ascending=False)
    
    print ('Top 5th percentile Transports')
    print ('---------------------')
    print (merge[merge.TRANSPORT_Rank > 0.95])
    print ('')
    print ('')
    print ('Bottom 5th percentile Transports')
    print ('---------------------')
    print (merge[merge.TRANSPORT_Rank < 0.05])


if __name__ == "__main__":
    main()
