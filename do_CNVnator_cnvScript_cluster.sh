#!/bin/bash

#download  the scripts from here
#http://alleleseq.gersteinlab.org/alleleSeq_cnvScript.zip
#then, create ONE SUCH DIRECTORY per sample (root sample)

#This script will create a customised .sh file in each of these directories and run it on the cluster

#STRUCTURE OF THE FILE
#SAMPLEID=NA12878                        ## individual's ID
#NAME=test                               ## e.g. hiseq_pcrfree_hc_130506
#FASTA=/path/to/fasta/                   ## path to FASTAs
#BINSIZE=100                             ## binsize e.g. 100 for high coverage (trio data), 1000 for low coverage (1KG data)
#ROOT=/path/to/NA12878.root              ## path of (tree) ROOT file from CNVnator
#SNP=/path/to/NA12878.snv                ## path to SNV file
#
### after running CNVnator on the BAM and produces a ROOT file, this generates the histogram from the ROOT file using CNVnator and the FASTAs
#ln -s $ROOT
#cnvnator -root $ROOT -outroot his.$SAMPLEID.$NAME.root -his $BINSIZE -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -d $FASTA
#cnvnator -root his.$SAMPLEID.$NAME.root -stat $BINSIZE
#cnvnator -root his.$SAMPLEID.$NAME.root -eval $BINSIZE > binSize$BINSIZE.log
#
### prepare addRD file
#rd=$(grep "Average RD per bin (1-22) is" binSize$BINSIZE.log | sed 's/Average RD per bin (1-22) is //g'  | awk '{printf("%d\n"),$SAMPLEID + 0.5}')
#./print_addRDcpp.sh $ROOT $rd./$BINSIZE
#make addRD
#
### run addRD
#ln -s $SNP
#./addRD $SNP >& rd.$SAMPLEID.$NAME.cnvnator.log


if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH_ROOT> <PATH_FA> <PATH_SNP> <BINSIZE>"
        echo "<PATH_ROOT> - absolute path for the .root files and CNVscripts directories (named by sample ID)"
	echo "<PATH_FA> absolute path for the .fa files by chromosome"
	echo "<PATH_SNP> - path to snp files in alleleseq format"
	echo "<BINSIZE> - eg 1000 for low depth 1kg sequence"
	echo "NOTE: the script will attempt to create ONE DIRECTORY per sample by COPYING THE alleleSeq_cnvScript directory"
	echo "SO MAKE SURE IT IS THERE"
        exit
fi

PDATA=$1;
PFASTA=$2;
PSNP=$3;
BIN=$4;
NAME="VDR_chip_exo";

PCODE="/net/isi-scratch/giuseppe/tools/CNVnator_v0.3/src";

for FILE in ${PDATA}/NA*.root;
        do
	ID=`echo ${FILE} | grep -Po "NA\d{5}"`; 

	#delete previously existing cnvscript directory and create new one by copying and renaming an existing alleleSeq_cnvScript directory
	THISDIR="${PDATA}/${ID}_cnvScript";
	if [ -d "$THISDIR" ]; then
		echo "Deleting $THISDIR ..";
	        rm -rf $THISDIR;
	fi	
	echo "Creating clean $THISDIR .."
	cp -R ${PDATA}/alleleSeq_cnvScript $THISDIR;
	#enter its cnvscript directory and prepare the file
	cd ${THISDIR};
	#mae script and launch it

        SCRIPT=CNVscript_${ID}.sh;
        echo '#!/bin/bash' >>${THISDIR}/${SCRIPT};
        echo '' >>${THISDIR}/${SCRIPT};
        echo 'source activate'  >>${THISDIR}/${SCRIPT};
        #export paths, just to be sure
        echo 'export ROOTSYS=/net/isi-scratch/giuseppe/tools/root' >>${THISDIR}/${SCRIPT};
        echo 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/net/isi-scratch/giuseppe/tools/root/lib' >>${THISDIR}/${SCRIPT};

	echo "SAMPLEID=${ID}"  >>${THISDIR}/${SCRIPT};
	echo "NAME=${NAME}" >>${THISDIR}/${SCRIPT};
	echo "FASTA=${PFASTA}" >>${THISDIR}/${SCRIPT};
	echo "BINSIZE=${BIN}" >>${THISDIR}/${SCRIPT};
	echo "ROOT=${FILE}" >>${THISDIR}/${SCRIPT};
	echo "THISDIR=${THISDIR}"  >>${THISDIR}/${SCRIPT};
	echo "PCODE=${PCODE}"  >>${THISDIR}/${SCRIPT};
	echo "SNP=${PSNP}/vcf2snp_${ID}.snv" >>${THISDIR}/${SCRIPT};

	## after running CNVnator on the BAM and produces a ROOT file, this generates the histogram from the ROOT file using CNVnator and the FASTAs
	echo "ln -s \$ROOT ${THISDIR}"  >>${THISDIR}/${SCRIPT};
	echo "\${PCODE}/cnvnator -root \$ROOT -outroot \${THISDIR}/his.\$SAMPLEID.\$NAME.root -his \$BINSIZE -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -d \$FASTA"  >>${THISDIR}/${SCRIPT};
	echo "\${PCODE}/cnvnator -root \${THISDIR}/his.\$SAMPLEID.\$NAME.root -stat \$BINSIZE"  >>${THISDIR}/${SCRIPT};
	echo "\${PCODE}/cnvnator -root \${THISDIR}/his.\$SAMPLEID.\$NAME.root -eval \$BINSIZE > \${THISDIR}/binSize\$BINSIZE.log"  >>${THISDIR}/${SCRIPT};

	## prepare addRD file
	echo 'rd=$(grep "Average RD per bin (1-22) is" ${THISDIR}/binSize1000.log | sed '\''s/Average RD per bin (1-22) is //g'\''  | awk '\''{printf("%d\n"),$SAMPLEID + 0.5}'\'')'  >>${THISDIR}/${SCRIPT};
	echo "\${THISDIR}/print_addRDcpp.sh \$ROOT \$rd./\$BINSIZE \${THISDIR}"  >>${THISDIR}/${SCRIPT};
	echo "make --directory \${THISDIR} -f \${THISDIR}/Makefile addRD"  >>${THISDIR}/${SCRIPT};

	## run addRD
	echo "ln -s \$SNP ${THISDIR}"  >>${THISDIR}/${SCRIPT};
	echo "\${THISDIR}/addRD \$SNP >& \${THISDIR}/rd.\$SAMPLEID.\$NAME.cnvnator.log"  >>${THISDIR}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/CNVscript_${ID}.err -o ${PDATA}/CNVscript_${ID}.out -q medium_jobs.q ${THISDIR}/${SCRIPT};
        #rm ${THISDIR}/${SCRIPT};  
done
