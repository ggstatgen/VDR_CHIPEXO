#!/bin/bash

#Create idx and hash for Stampy aligner
#/net/isi-scratch/giuseppe/tools/stampy-1.0.20/stampy.py -G hg19_stampy hg19.fa
#/net/isi-scratch/giuseppe/tools/stampy-1.0.20/stampy.py -g hg19_stampy -H hg19_stampy


# -G PREFIX, --build-genome=PREFIX   Build genome index PREFIX.stidx from fasta file(s) on command line
# -H PREFIX, --build-hash=PREFIX     Build hash PREFIX.sthash
# -g PREFIX, --genome=PREFIX         Use genome index file PREFIX.stidx
# -h PREFIX, --hash=PREFIX           Use hash file PREFIX.sthash


#Build a genome (.stidx) file:
#     ./stampy.py --species=human --assembly=hg18_ncbi36 \
#                 -G hg18 /data/genomes/hg18/*.fa.gz
#
#Build a hash (.sthash) file:
#
#     ./stampy.py -g hg18 -H hg18


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <DATA_DIR> <PREFIX>"
        echo "DATA_DIR - full path for the reference genome(s) to index with Stampy"
        echo "PREFIX -  prefix for the outputfiles - eg hg19_stampy"
        exit
fi

PCODE="/net/isi-scratch/giuseppe/tools/stampy-1.0.21";
PPYTHON="~/src/pypa-virtualenv-3a798b3/ve/bin/python";
PDATA=$1;
PPREFIX=$2;

for FILE in ${PDATA}/*.fa;
        do
        BASEFILE=`basename ${FILE} ".fa"`;
        SCRIPT=stampy_${BASEFILE}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
        
        #build idx
        echo "${PPYTHON} ${PCODE}/stampy.py -G ${PDATA}/${PPREFIX} ${FILE}" >>${PDATA}/${SCRIPT};
        #build hash
        echo "${PPYTHON} ${PCODE}/stampy.py -g ${PDATA}/${PPREFIX} -H ${PDATA}/${PPREFIX}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/stampy_${BASEFILE}.err -o ${PDATA}/stampy_${BASEFILE}.out  -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
