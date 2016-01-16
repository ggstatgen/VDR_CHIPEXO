import sys
import HTSeq
import collections
import argparse
__author__ = 'Giuseppe Gallone'

assert sys.version_info[:2] >= ( 2, 4 )

def main():
	parser = argparse.ArgumentParser(description='Counts raw reads in bed intervals using the HTSeq interface.')
	parser.add_argument('--bed',help='Bed input file name',required=True)
	parser.add_argument('--bam',help='Bam input file name',required=True)
	args = parser.parse_args()

	input_bam = args.bam
	input_bed = args.bed
	#input_bam = sys.argv[1]
    	#input_bed = sys.argv[2]

	features =  HTSeq.GenomicArrayOfSets( "auto", stranded=False )
	counts = collections.Counter( )
	#almnt_file = HTSeq.SAM_Reader( input_sam )
	almnt_file = HTSeq.BAM_Reader( input_bam )

	for line in open( input_bed ):
		line = line.rstrip("\n")
		fields = line.split( "\t" )
		iv = HTSeq.GenomicInterval( fields[0], int(fields[1]), int(fields[2]) )
		features[ iv ] += fields[3]
	
	for almnt in almnt_file:
		if not almnt.aligned:
			count[ "_unmapped" ] += 1
			continue
		gene_ids = set()
		for iv, val in features[ almnt.iv ].steps():
			gene_ids |= val
		if len(gene_ids) == 1:
			gene_id = list(gene_ids)[0]
			counts[ gene_id ] += 1
			counts[ "_mapped" ] += 1
		elif len(gene_ids) == 0:
			counts[ "_no_feature" ] += 1
		else:
			counts[ "_ambiguous" ] += 1

	for gene_id in counts:
   		print gene_id, counts[ gene_id ]
	
if __name__ == "__main__":
    main()
