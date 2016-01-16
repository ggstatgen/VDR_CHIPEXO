#!/bin/bash
#script to do liftover from hg18 to hg19 for multiple bed files
#eg /net/isi-backup/giuseppe/scripts/liftOver GM_consensus_slop_1000.bed /net/isi-scratch/giuseppe/GATK_RESOURCES/chain/hg19tob37.chain GM_consensus_slop_1000_liftover19_37.bed unlifted_1000.bed

#liftOver - Move annotations from one assembly to another
#usage:
#   liftOver oldFile map.chain newFile unMapped
#oldFile and newFile are in bed format by default, but can be in GFF and
#maybe eventually others with the appropriate flags below.
#The map.chain file has the old genome as the target and the new genome
#as the query.
#
#***********************************************************************
#WARNING: liftOver was only designed to work between different
#         assemblies of the same organism, it may not do what you want
#         if you are lifting between different organisms.
#***********************************************************************
#
#options:
#   -minMatch=0.N Minimum ratio of bases that must remap. Default 0.95
#   -gff  File is in gff/gtf format.  Note that the gff lines are converted
#         separately.  It would be good to have a separate check after this
#         that the lines that make up a gene model still make a plausible gene
#         after liftOver
#   -genePred - File is in genePred format
#   -sample - File is in sample format
#   -bedPlus=N - File is bed N+ format
#   -positions - File is in browser "position" format
#   -hasBin - File has bin value (used only with -bedPlus)
#   -tab - Separate by tabs rather than space (used only with -bedPlus)
#   -pslT - File is in psl format, map target side only
#   -minBlocks=0.N Minimum ratio of alignment blocks or exons that must map
#                  (default 1.00)
#   -fudgeThick    (bed 12 or 12+ only) If thickStart/thickEnd is not mapped,
#                  use the closest mapped base.  Recommended if using 
#                  -minBlocks.
#   -multiple               Allow multiple output regions
#   -minChainT, -minChainQ  Minimum chain size in target/query, when mapping
#                           to multiple output regions (default 0, 0)
#   -minSizeT               deprecated synonym for -minChainT (ENCODE compat.)
#   -minSizeQ               Min matching region size in query with -multiple.
#   -chainTable             Used with -multiple, format is db.tablename,
#                               to extend chains from net (preserves dups)
#   -errorHelp              Explain error messages


if [ ! $# == 4  ]; then
	echo "Usage: `basename $0` <PATH> <EXT> <CHAIN> <STRING>"
        echo "<PATH> Directory containing the data files (e.g. /home/me/files)"
	echo "<EXT> file extension (e.g. bed, bedgraph, bdg)"
	echo "<CHAIN> full path to chain file to use"
	echo "<STRING> mnemonic string to use for lifted over files (eg hg19)"
	exit
fi

PDATA=$1;
EXT=$2;
PCHAIN=$3;
PSTRING=$4;

PCODE="/net/isi-backup/giuseppe/scripts";

for FILE in ${PDATA}/*.${EXT};
        do 
        ID=`basename ${FILE} ".${EXT}"`;
	
	SCRIPT=liftover_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	echo "${PCODE}/liftOver ${FILE} ${PCHAIN} ${PDATA}/${ID}_${PSTRING}.${EXT} ${PDATA}/${ID}_unmapped.${EXT}" >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/liftover_${ID}.err -o ${PDATA}/liftover_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
