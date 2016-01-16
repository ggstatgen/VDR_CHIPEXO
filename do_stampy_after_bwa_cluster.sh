#!/bin/bash
#script to run stampy on the output of bwa, on the cluster
#input: bam files obtained from bwa.
#output: .bam files, sorted bam files, bam indexes
#the stampy doc says stampy won't need sorted bams from bwa


#command
#stampy.py -g <PATH>/hg19 -h <PATH>hg19 -t8 --bamkeepgoodreads -M bwa.bam
#./stampy.py --sensitive --bwa=/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa -g mm9 -h mm9 -t10 --bamkeepgoodreads -M bwa.bam

PCODE="/net/isi-scratch/giuseppe/tools";
PINDEX="/net/isi-scratch/giuseppe/indexes/stampy/hg19_stampy";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/hg19_stampy";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Mmus/mm10/mm10_stampy";
PPYTHON="~/src/pypa-virtualenv-3a798b3/ve/bin/python";
PBWA="/home/giuseppe/local/bin" #this refers to a symbolic link to the latest bwa. If you don't use this, it will try to use the cgat bwa
PSAMTOOLS="/home/giuseppe/local/bin/";

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the bam files from bwa"
        exit
fi

PDATA=$1;
POUT=${PDATA}/d_STAMPY;
mkdir ${POUT};

for FILE in ${PDATA}/*.bam;
	do
 
        BASEFILE=`basename ${FILE} ".bam"`; 
	SCRIPT=stampy_${BASEFILE}.sh;
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
	#change required options here:
        #echo "${PCODE}/stampy-1.0.21/stampy.py --sensitive --solexa --bwa=/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa -g ${PINDEX} -h ${PINDEX} -t10 --bamkeepgoodreads -M ${FILE}" >>${PDATA}/${SCRIPT};
        
	#echo "${PPYTHON} ${PCODE}/stampy-1.0.21/stampy.py --sensitive --bwa=${PBWA}/bwa -g ${PINDEX} -h ${PINDEX} -t10 --bamkeepgoodreads -M ${FILE}" >>${PDATA}/${SCRIPT};        
#	echo "${PPYTHON} ${PCODE}/stampy-1.0.21/stampy.py --bwa=${PBWA}/bwa -g ${PINDEX} -h ${PINDEX} -t10 --bamkeepgoodreads -M ${FILE}" >>${PDATA}/${SCRIPT};

        echo "${PPYTHON} ${PCODE}/stampy-1.0.28/stampy.py --bwa=${PBWA}/bwa -g ${PINDEX} -h ${PINDEX} -t10 --bamkeepgoodreads -M ${FILE}" >>${PDATA}/${SCRIPT};


        echo "${PSAMTOOLS}/samtools view -bS ${FILE} > ${PDATA}/${BASEFILE}_stampy.bam;" >> ${PDATA}/${SCRIPT};
        #echo "${PSAMTOOLS}/samtools sort ${PDATA}/${BASEFILE}.bam ${PDATA}/${BASEFILE}.sorted;" >>${PDATA}/${SCRIPT};
        #echo "${PSAMTOOLS}/samtools index ${PDATA}/${BASEFILE}.sorted.bam;" >>${PDATA}/${SCRIPT};

	#notice: if you specify h_vmem, Stampy won't work
	#nice -5 qsub -l h_vmem=5G -pe dedicated 8 -v "BASH_ENV=~/.bashrc" -e ${POUT}/stampy_${BASEFILE}.err -o ${POUT}/${BASEFILE}.sam -q medium_jobs.q ${PDATA}/${SCRIPT};
	nice -5 qsub -pe dedicated 10 -v "BASH_ENV=~/.bashrc" -e ${POUT}/stampy_${BASEFILE}.err -o ${POUT}/${BASEFILE}.sam -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
done
