# Order parameter calculator ordermap

A program for calculating order parameter heatmaps

## Basic usage

Either add the root folder to your PATH and run

'''
ordermap -h
'''

to get more info, or just run

'''
/path/to/folder/ordermap/ordermap -h
'''

Example run:
'''
ordermap -f trajectories/md_fitted.xtc -s trajectories/start_nopbc.tpr -o order/order \
         -center protein -leafdiv -sel1 "POPC C 29-42" -timeseries -sepcarbs
'''
