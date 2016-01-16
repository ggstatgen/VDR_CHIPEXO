#!/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin/python

#Scripts to randomly subsample single end reads from fastq files

#usage
#python subsample.py <fraction> <input file> <output file>

import sys, random, itertools
import HTSeq

fraction = float( sys.argv[1] )
in_fasta  = iter( HTSeq.FastqReader( sys.argv[2] ) )
out_fasta = open( sys.argv[3], "w" )

for read, in itertools.izip( in_fasta ):
   if random.random() < fraction:
      read.write_to_fastq_file( out_fasta )
      
out_fasta.close()
