#!/bin/bash

#script to run XXmotif on the cluster
#http://xxmotif.genzentrum.lmu.de/

#Usage: XXmotif OUTDIR SEQFILE [options] 
#
#	OUTDIR:  output directory for all results
#	SEQFILE: file name with sequences from positive set in FASTA format
#
#Options:
#	--negSet <FILE>				sequence set which has to be used as a reference set
#	--zoops					use zero-or-one occurrence per sequence model (DEFAULT)
#	--mops					use multiple occurrence per sequence model
#	--oops					use one occurrence per sequence model
#	--revcomp				search in reverse complement of sequences as well (DEFAULT: NO)
#	--background-model-order <NUMBER>	order of background distribution (DEFAULT: 2, 8(--negset) )
#	--pseudo <NUMBER>			percentage of pseudocounts used (DEFAULT: 10)
#	-g|--gaps <NUMBER>			maximum number of gaps used for start seeds [0-3] (DEFAULT: 0)
#	--type <TYPE>				defines what kind of start seeds are used (DEFAULT: ALL)
#						 - possible types: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM
#	--merge-motif-threshold <MODE>			defines the similarity threshold for merging motifs (DEFAULT: HIGH)
#						 - possible modes: LOW, MEDIUM, HIGH
#	--no-pwm-length-optimization		do not optimize length during iterations (runtime advantages)
#	--max-match-positions <INT>		max number of positions per motif (DEFAULT: 17, higher values will lead to very long runtimes)
#
#	--batch					suppress progress bars (reduce output size for batch jobs)
#	--maxPosSetSize <NUMBER>		maximum number of sequences from the positive set used [DEFAULT: all]
#	--help					print this help page
#	--trackedMotif <SEED>			inspect extensions and refinement of a given seed (DEFAULT: not used)
#
#Using conservation information
#	--format FASTA|MFASTA			defines what kind of format the input sequences have (DEFAULT: FASTA)
#	--maxMultipleSequences <NUMBER>		maximum number of sequences used in an alignment [DEFAULT: all]
#
#Using localization information
#	--localization				use localization information to calculate combined P-values 
#						(sequences should have all the same length)
#	--downstream <NUMBER>			number of residues in positive set downstream of anchor point (DEFAULT: 0)
#
#Start with self defined motif:
#	-m|--startMotif <MOTIF>			Start motif (IUPAC characters)
#	-p|--profileFile <FILE>			profile file
#	--startRegion <NUMBER>			expected start position for motif occurrences relative to anchor point (--localization)
#	--endRegion <NUMBER>			expected end position for motif occurrences relative to anchor point (--localization)



if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <M_THRESH>"
        echo "<PATH> data path for the .fa files"
	echo "<M_THRES> motif similarity threshold [LOW|MEDIUM|HIGH]"
        exit
fi

PDATA=$1;
SIM_TH=$2;
PCODE="/net/isi-scratch/giuseppe/tools/XXmotif";
POUT=${PDATA}/d_XXmotif_th_${SIM_TH};
mkdir ${POUT};

for FILE in ${PDATA}/*.fa;
	do ID=`basename ${FILE} ".fa"`;

        SCRIPT=XXmotif_${ID}_th_${SIM_TH}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	#echo 'source activate' >>${PDATA}/${SCRIPT};

	# --merge-motif-threshold ${SIM_TH}
        echo "${PCODE}/XXmotif ${POUT} ${FILE} --zoops --revcomp --localization  --batch --merge-motif-threshold ${SIM_TH}" >>${PDATA}/${SCRIPT};

        nice -5  qsub -v "BASH_ENV=~/.bashrc" -e ${POUT}/XXmotif_${ID}_th_${SIM_TH}.err -o  ${POUT}/XXmotif_${ID}_th_${SIM_TH}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        #rm ${PDATA}/${SCRIPT};  
done
