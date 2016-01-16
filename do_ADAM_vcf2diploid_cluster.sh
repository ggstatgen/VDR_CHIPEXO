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

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <VCF_PATH>"
        echo "<VCF_PATH> path to dir containing .vcf files";
	echo "NOTE: using /net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/by_chr/ vcf for the reference. NOTE THESE MUST BE UNZIPPED, unzip if necessary";
        exit
fi

PDATA=$1;

#1000 genomes reference by chromosome. You will need to append the actual file, which is in the format chr11.fa, chr12.fa...
PREFDIR="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/by_chr";
PCODE="/net/isi-scratch/giuseppe/tools/vcf2diploid";


for FILE in ${PDATA}/*.vcf;
	do
	ID=`echo ${FILE} | grep -Po "(chr\d+)"`;
	BASENAME=`basename ${FILE} ".vcf"`;
	#create reference file
	PREF=${PREFDIR}/${ID}.fa;

        POUT=${PDATA}/${ID};
        mkdir ${POUT};

	SCRIPT=vcf2dip_${ID}.sh;

	echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};

	#unzip vcf
	echo "cd ${POUT}" >> ${POUT}/${SCRIPT};
        echo "java -jar ${PCODE}/vcf2diploid.jar \\
              -id  SAMPLE_c \\
              -chr ${PREF}  \\
              -vcf ${PDATA}/${BASENAME}.vcf" >>${POUT}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/vcf2diploid_${ID}.err -o ${PDATA}/vcf2diploid_${ID}.out -q fgu217.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT};
	unset POUT;
done
