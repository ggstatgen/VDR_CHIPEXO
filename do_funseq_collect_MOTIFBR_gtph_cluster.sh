#!/bin/bash

#wrapper to run do_funseq_collect_MOTIFBR_gtph_REVIEW.pl on the cluster

#Input: a directory full of files created by do_pscanchip_out_intersect_vdrbvs.pl, in PscanCHIP format, each for a different TF whose PWM intervals around your sequences you've chosen to download

#This script will convert those files into a format suitable for Funseq. Collate the outputs of this and rename "ENCODE.tf.bound.union.bed", then save in the /data dir of Funseq2 (backing up the previous version)

#synopsis
#USAGE: do_funseq_collect_MOTIFBR_gtph.pl -i=<INFILE> -m=<MOTIF_NAME_ID> -subset=<SUBSET_BED> -sym -enh
#<INFILE> vcf output of funseq2
#<MOTIF_NAME_ID> Jaspar name and Jaspar id for the pwm to consider. Must be in the form <name_id>, e.g.: NFYA_MA0060.1
#(optional)<SUBSET_BED> only collect data for VDR-BVs breaking motifs intersecting with this bed file(default=no)
#(optional)<sym> flag; set it if you are dealing with SYM variants from alleleseq (default=no)
#(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)


if [ ! $# == 4 ]; then
	echo "Usage: `basename $0` <FUNSEQ_OUT_VCF_PATH> <PWM_PATH> <JASPAR_ENCODE_FILE> <ASYM|SYM>"
	echo "<FUNSEQ_OUT_VCF> - path to the vcf file with the annotated variants"
	echo "<PWM_PATH> - absolute path for the pscanchip files (.ris extension), one per TF"
	echo "<JASPAR_ENCODE_FILE> file containing the pwm desired in the ENCODE format, converted from the Jaspar format with RSA-tools (see gdrive for details)"        
	echo "<ASYM|SYM> one of [asym|sym]"
	echo "NOTE1: the vcf file is assumed to be called Output_noDBRECUR.vcf. Modify if not";
	echo "NOTE2: this will only work with filenames with the structure Pscanchip_allsites_BVs_<TF>_<ID>_sites.<EXT>";
	echo "so modify the filename accordingly";
        exit
fi

PVCF=$1;
PDATA=$2;
PJASPARPWMS=$3;
PSYMASYM=$4;

FUNSEQ_FILE="Output_noDBRECUR.vcf";
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_funseq_collect_MOTIFBR_gtph_REVIEW.pl";

for FILE in ${PDATA}/Pscanchip*.ris;
        do
	SEED=`echo ${FILE} | grep -Po "BVs_(.*)_sites"`;
	TF=`echo ${SEED} | cut -d _ -f 2`;
	ID=`echo ${SEED} | cut -d _ -f 3`;

	SCRIPT=script_collectMOTIFBR_${TF}_${ID}_${PSYMASYM}.sh;

        echo '#!/bin/bash' >>${PVCF}/${SCRIPT};
        echo '' >>${PVCF}/${SCRIPT};
	#echo 'source activate'  >>${PVCF}/${SCRIPT};

	if [ "${PSYMASYM}" = "asym"  ]
	then
		echo "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl ${PCODE} -i=${PVCF}/${FUNSEQ_FILE} -m=${TF}_${ID} -j=${PJASPARPWMS}" >>${PVCF}/${SCRIPT};
	elif [ "${PSYMASYM}" = "sym"  ]
	then
		echo "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl ${PCODE} -i=${PVCF}/${FUNSEQ_FILE} -m=${TF}_${ID} -j=${PJASPARPWMS} -sym" >>${PVCF}/${SCRIPT};
	else
		echo "field <ASYM|SYM> not recognised. Aborting."
		exit 1
	fi
	
	nice -5 qsub -e ${PVCF}/script_collectMOTIFBR_${TF}_${ID}_${PSYMASYM}.err -o ${PVCF}/script_collectMOTIFBR_${TF}_${ID}_${PSYMASYM}.out -q medium_jobs.q ${PVCF}/${SCRIPT};
        rm ${PVCF}/${SCRIPT};  
done
