####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: gff2bed.py 2861 2010-02-23 17:36:32Z andreas $
##
##
####
####
'''
gff2bed.py - convert from gff to bed
====================================

Purpose
-------

convert gff to bed formatted files.

Usage
-----

Type::

   python gff2bed.py --help

for command line usage.

Code
----
'''

import sys, re, string, optparse, time, os, itertools
import Experiment as E
import GFF, GTF, Bed

def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: gff2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage=globals()["__doc__"])
    
    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                       help="input file is gtf [default=%default] ")

    # parser.add_option( "--format", dest="format", type="choice",
    #                    help = "bed output format [%default]",
    #                    choices = ("bed3", "bed4", "bed6" ) )

    parser.add_option( "--name", dest="name", type="choice",
                       help = "field to use as the name field [%default]",
                       choices = ("gene_id", "transcript_id", "class", "family", "feature", "source" ) )
    
    parser.add_option( "--track", dest="track", type="choice",
                       choices=("feature", "source", None),
                       help="use feature/source field to define tracks [default=%default] ")

    parser.set_defaults(
        track = None,
        # format = "bed4",
        name = "gene_id",
        is_gtf = False )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    ninput, noutput = 0, 0
    
    is_gtf, name = options.is_gtf, options.name

    if is_gtf:
        iterator = GTF.iterator( options.stdin )
    else:
        iterator = GFF.iterator( options.stdin )

    if options.track:
        all_input = list( iterator )

        if options.track == "feature":
            grouper = lambda x: x.feature
        elif options.track == "source":
            grouper = lambda x: x.source

        all_input.sort( key = grouper )

        bed = Bed.Bed()
        for key, vals in itertools.groupby( all_input, grouper ):
            options.stdout.write( "track name=%s\n" % key )
            for gff in vals:
                ninput += 1
                bed.fromGFF( gff, is_gtf = is_gtf, name = name )
                options.stdout.write( str(bed) + "\n" )
                noutput += 1

    else:
        bed = Bed.Bed()
        for gff in iterator:
            ninput += 1
            bed.fromGFF( gff, is_gtf = is_gtf, name = name )
            options.stdout.write( str(bed) + "\n" )

            noutput += 1

    E.info( "ninput=%i, noutput=%i" % (ninput, noutput) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
