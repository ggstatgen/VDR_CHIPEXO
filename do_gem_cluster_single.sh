#!/bin/bash

#script to run GEM on the clusters

#hg19 mappability sizes
#from https://groups.google.com/forum/?fromgroups=#!topic/macs-announcement/9DeYzN3vFWA
# 40	2540757438


#typical command line 
#java -Xmx10G -jar gem.jar --g /net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes --d Read_Distribution_ChIP-exo.txt --s 2540757438 --expt /net/isi-scratch/giuseppe/VDR/MY_FASTQ_CUTADAPT/BAM_BOWTIE/VDR_GM06986_Peconic20301_trimmed.sorted.bam --f SAM --outBED --out /net/isi-scratch/giuseppe/VDR/MY_FASTQ_CUTADAPT/BAM_BOWTIE/GM06986 --genome /net/isi-scratch/giuseppe/indexes/Hsap/hg19_genome --k 

#other settings worth playing with
#--k_neg_dinu_shuffle -  Use di-nucleotide shuffled sequences as negative sequences for motif finding
# instead of using a piece of sequence 300 bp from binding, it uses a shuffled set

#--seed - the seed k-mer to jump start k-mer class discovery. The width of the seed k-mer will be used to set k 

# --fold [value]	Fold cutoff to filter predicted events (default=3)

# --q [value]	significance level for q-value, specified as -log10(q-value). For example, to enforce a q-value threshold of 0.001, set this value to 3. (default=2, i.e. q-value=0.01)

# --sd [value]	Shape deviation cutoff to filter predicted events (default=-0.40).

# --w2 [n]	Size of sliding window to estimate lambda parameter for Possion distribution when there is no control data (default=5,000, must be larger than 1000).

# --w3 [n]	Size of sliding window to esitmate lambda parameter for Possion distribution when there is no control data (default=10,000, must be larger than w2).

#For ChIP-exo data, add one more option --smooth 3 to estimate the read distribution more accurately.
# --sl Sort GEM output by location (default is sorted by P-value)

#high memory requirements, so use fgu217


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <MIN_KMER> <MAX_KMER>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files) and the minimum and maximum k-mer length"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools";
PCHROMSIZE="/net/isi-scratch/giuseppe/indexes/chrominfo";
#PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_genome"; 
PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/chromosomes_hg19"; #fasta files by chromosome


for FILE in ${PDATA}/*.bam;
	do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        #do ID=`echo ${FILE} | egrep -o "38380_[0-9]*"`;        

        SCRIPT=gem_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        #--q 3 \\
        echo "java -Xmx10G -jar ${PCODE}/gem/gem.jar \\
              --g ${PCHROMSIZE}/hg19.chrom.sizes \\
              --d ${PCODE}/gem/Read_Distribution_ChIP-exo.txt \\
              --s 2540757438 \\
              --expt ${FILE} \\
              --f SAM \\
              --outBED \\
              --smooth 3 \\
              --sl \\
              --out ${PDATA}/${ID} \\
              --genome ${PGENOME} \\
              --k_min $2 \\
              --k_max $3" >>${PDATA}/${SCRIPT};
        #nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q fgu217.q ${PDATA}/${SCRIPT};
        nice  -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q fgu205.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
