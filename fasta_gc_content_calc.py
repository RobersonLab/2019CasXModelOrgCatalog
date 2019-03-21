#!/usr/bin/env python

##########
# Import #
##########
import argparse
import pyfaidx
import sys

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split( '/' )[-1].split( '\\' )[-1]
VERSION = '1.0.1'

if __name__ == '__main__':
	############
	# argparse #
	############
	parser = argparse.ArgumentParser( prog=SCRIPT_NAME, epilog="%s v%s" % ( SCRIPT_NAME, VERSION ) )
	
	parser.add_argument( 'fasta_file', help="Path to indexed FASTA file to be scanned for gRNA sites" )
	parser.add_argument( '--output_file', help="Defaults to gc_percent.csv", default="gc_percent.csv" )
	parser.add_argument( "--quiet", default=True, action='store_false', dest='verbose' )

	args = parser.parse_args()

	if args.verbose:
		print "%s v%s\n\nOptions\n=======" % ( SCRIPT_NAME, VERSION )
		print "FASTA: %s" % ( args.fasta_file )
		print "Output file: %s" % ( args.output_file )
		print "Verbose: %s" % ( str( args.verbose ) )
	
	#######################
	# Read the FASTA file #
	#######################
	gc_bases = 0
	total_good_bases = 0
	
	with pyfaidx.Fasta( args.fasta_file, as_raw=True ) as FASTA:
		for sequence in FASTA:
			sequence = str( sequence ).upper()
			count_it = sequence.count
			
			gc = sum( map( count_it, ('G', 'C') ) )
			at = sum( map( count_it, ('A', 'T') ) )
			
			total_good_bases += gc + at
			gc_bases += gc
	
	with open( args.output_file, 'w' ) as OUTFH:
		OUTFH.write( "%.4f" % ( float( gc_bases ) / float( total_good_bases ) ) )
