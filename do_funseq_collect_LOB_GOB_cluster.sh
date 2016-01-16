#!/bin/bash

#wrapper to run do_funseq_collect_MOTIFBR_gtph_REVIEW.pl on the cluster

#synopsis
#<INFILE> vcf output of funseq2
#<MOTIF_NAME_ID> Jaspar name and Jaspar id for the pwm to consider. Must be in the form <name_id> e.g.: NFYA_MA0060.1
#<JASPAR_ENCODE_PWM> file containing the pwm desired in the ENCODE format, converted from the Jaspar format with RSA-tools (see gdrive for details)
#(optional)<SUBSET_BED> only collect data for VDR-BVs breaking motifs intersecting with this bed file(default=no)
#(optional)<sym> flag; set it if you are dealing with SYM variants from alleleseq (default=no)(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)


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
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_funseq_collect_LOB_GOB_REVIEW.pl";

for FILE in ${PDATA}/Pscanchip*.ris;
        do
	SEED=`echo ${FILE} | grep -Po "BVs_(.*)_sites"`;
	TF=`echo ${SEED} | cut -d _ -f 2`;
	ID=`echo ${SEED} | cut -d _ -f 3`;

	SCRIPT=script_collect_LOB_GOB_${TF}_${ID}_${PSYMASYM}.sh;

        echo '#!/bin/bash' >>${PVCF}/${SCRIPT};
        echo '' >>${PVCF}/${SCRIPT};

	#for SYM variants
	if [ "${PSYMASYM}" = "asym" ]
	then
		echo "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl ${PCODE} -i=${PVCF}/${FUNSEQ_FILE} -m=${TF}_${ID} -j=${PJASPARPWMS} -sym" >>${PVCF}/${SCRIPT};
	elif [  "${PSYMASYM}" = "sym" ]
	then
	#for ASYM variants
		echo "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl ${PCODE} -i=${PVCF}/${FUNSEQ_FILE} -m=${TF}_${ID} -j=${PJASPARPWMS}" >>${PVCF}/${SCRIPT};
	else
		echo "field <ASYM|SYM> not recognised. Aborting."
		exit 1
	fi	

	nice -5 qsub -e ${PVCF}/script_collect_LOB_GOB_${TF}_${ID}_${PSYMASYM}.err -o ${PVCF}/script_collect_LOB_GOB_${TF}_${ID}_${PSYMASYM}.out -q medium_jobs.q ${PVCF}/${SCRIPT};
        rm ${PVCF}/${SCRIPT};  
done
