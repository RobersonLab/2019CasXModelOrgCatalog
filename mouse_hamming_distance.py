#!/usr/bin/env python

from scipy.spatial import distance
import pandas as pd

#########
# input #
#########
input = pd.read_csv( "R_analysis/mouse_target_upset.csv.gz" )
input = input.drop( columns = [ 'Target' ] )

mouse_names = list( input.columns )
mouse_number = input.shape[1]

###############################
# setup russell-rao dataframe #
###############################
out = pd.DataFrame( columns = mouse_names, index = mouse_names )

######################################
# loop through and calculate indexes #
######################################
for row_idx in range( mouse_number ):
	for col_idx in range( mouse_number ):
		mouse_rowname = mouse_names[ row_idx ]
		mouse_colname = mouse_names[ col_idx ]

		dist = distance.hamming( input[ mouse_rowname ], input[ mouse_colname ] )
		out[ mouse_colname ][ mouse_rowname ] = dist
out.to_csv( "R_analysis/mouse_hamming_distance.csv", header=True, index=True )

