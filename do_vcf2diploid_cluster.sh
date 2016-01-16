#!/bin/sh

#java -jar vcf2diploid.jar -id sample_id -chr file.fa ... [-vcf file.vcf ...]

#where sample_id is the ID of individual whose genome is being constructed
#(e.g., NA12878), file.fa is FASTA file(s) with reference sequence(s), and
#file.vcf is VCF4.0 file(s) with variants. One can specify multiple FASTA and
#VCF files at a time. Splitting the whole genome in multiple files (e.g., with
#one FASTA file per chromosome) reduces memory usage.
#Amount of memory used by Java can be increased as follows

#java -Xmx4000m -jar vcf2diploid.jar -id sample_id -chr file.fa ... [-vcf file.vcf ...]

#Important notes
#===============
#
#All characters between '>' and first white space in FASTA header are used
#internally as chromosome/sequence names. For instance, for the header
#
#>chr1 human
#
#vcf2diploid will upload the corresponding sequence into the memory under the
#name 'chr1'.
#Chromosome/sequence names should be consistent between FASTA and VCF files but
#omission of 'chr' at the beginning is allows, i.e. 'chr1' and '1' are treated as
#the same name.
#
#The output contains (file formats are described below):
#1) FASTA files with sequences for each haplotype.
#2) CHAIN files relating paternal/maternal haplotype to the reference genome.
#3) MAP files with base correspondence between paternal-maternal-reference
#sequences.
#File formats:
#* FASTA -- see http://www.ncbi.nlm.nih.gov/blast/fasta.shtml
#* CHAIN -- http://genome.ucsc.edu/goldenPath/help/chain.html
#* MAP file represents block with equivalent bases in all three haplotypes
#(paternal, maternal and reference) by one record with indices of the first
#bases in each haplotype. Non-equivalent bases are represented as separate
#records with '0' for haplotypes having non-equivalent base (see
#clarification below).
#
#Pat Mat Ref           MAP format
#X   X   X    ____
#X   X   X        \
#X   X   X         --> P1 M1 R1
#X   X   -    -------> P4 M4  0
#X   X   -       ,--->  0 M6 R4
#-   X   X    --'  ,-> P6 M7 R5
#X   X   X    ----'
#X   X   X
#X   X   X

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <VCF_PATH> <REF> <DESC>"
        echo "<VCF_PATH> path to dir containing .vcf.gz files, one per sample, or one per trio"
        echo "<REF> FASTA file with reference sequence (can be full/chr masked/nonmasked)"
	echo "<DESC> short 1-word description of the ref/vcf (eg hg19)"
	echo "THIS DOESNT WORK FOR TRIO VCF FILES. ID MUST BE THE CHILD'S ID ONLY USE TO CREATE SCRIPT THEN SEND MANUALLY TO GSUB"
        exit
fi

PDATA=$1;
PREF=$2;
PCODE="/net/isi-scratch/giuseppe/tools/vcf2diploid";
PDESC=$3;

POUT1=${PDATA}/VCF2DIPLOID_${PDESC};
mkdir ${POUT1};

for FILE in ${PDATA}/*.vcf.gz;
	do BASENAME=`basename ${FILE} ".vcf.gz"`;
	ID=`echo ${FILE} | egrep -o "NA[0-9]*"`; #single
	#ID=`echo ${FILE} | egrep -o "YRI_Y[0-9]*"`; #trios

	POUT=${POUT1}/${ID};
	mkdir ${POUT};

	SCRIPT=vcf2diploid_${BASENAME}.sh;

	echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};

	#unzip vcf
	echo "gunzip -f ${FILE}" >>${POUT}/${SCRIPT};
	echo "cd ${POUT}" >> ${POUT}/${SCRIPT};
        echo "java -jar ${PCODE}/vcf2diploid.jar \\
              -id  ${ID}    \\
              -chr ${PREF}  \\
              -vcf ${PDATA}/${BASENAME}.vcf" >>${POUT}/${SCRIPT};
	#zip vcf
	echo "gzip ${PDATA}/${BASENAME}.vcf" >> ${POUT}/${SCRIPT};

        nice -5 qsub -e ${POUT1}/vcf2diploid_${ID}.err -o ${POUT1}/vcf2diploid_${ID}.out -q medium_jobs.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT};
	unset POUT;
done
