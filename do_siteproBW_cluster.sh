#!/bin/bash

#this is to run siteproBW on the cluster
#Usage: siteproBW <-w bigwig -b bed> [options]
#
#sitepro -- Average profile around given genomic sites
#
#Options:
#  --version             show program's version number and exit
#  -h, --help            Show this help message and exit.
#  -w WIG, --bw=WIG      input bigWIG file. Multiple bigWIG files can be given
#                        via -w (--bw) individually (eg -w WIG1.bw, -w
#                        WIG2.bw). WARNING! multiple bigwig and bed files are
#                        not allowed.
#  -b BED, --bed=BED     BED file of regions of interest. (eg, binding sites or
#                        motif locations) Multiple BED files can be given via
#                        -b (--bed) individually (eg -b BED1.bed -b BED2.bed).
#                        WARNING! multiple wig and bed files are not allowed.
#  --span=SPAN           Span from the center of each BED region in both
#                        directions(+/-) (eg, [c - span, c + span], where c is
#                        the center of a region), default:1000 bp
#  --pf-res=PF_RES       Profiling resolution, default: 50 bp
#  --dir                 If set, the direction (+/-) is considered in
#                        profiling. If no strand info given in the BED, this
#                        option is ignored.
#  --dump                If set, profiles are dumped as a TXT file
#  --confid              If set, it will draw 95% confidence interval for each
#                        step.
#  --name=NAME           Name of this run. If not given, the body of the bed
#                        file name will be used,
#  -l LABEL, --label=LABEL
#                        Labels of the wig files. If given, they are used as
#                        the legends of the plot and in naming the TXT files of
#                        profile dumps; otherwise, the bigWIG file names will
#                        be used as the labels. Multiple labels can be given
#                        via -l (--label) individually (eg, -l LABEL1 -l
#                        LABEL2). WARNING! The number and order of the labels
#                        must be the same as the bigWIG files.

#one idea is to use a consensus .bed (for example the 1000 consensus I've been using with GATK) and then one bw for each file


if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <BW_PATH> <INTERVAL_FILE> <RES> <SPAN>"
        echo "BW_PATH - the directory containing the .bw files"
        echo "INTERVAL_FILE - the .bed used to build the model"
	echo "RES - Profiling resolution, (eg 50 bp)"
	echo "SPAN - Span from the center of each BED region, (eg 1000 bp)"
        exit
fi

PDATA=$1;
INTERVAL=$2;
RES=$3;
SPAN=$4;
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
#Got this from the CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html

for FILE in ${PDATA}/*.bw;
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        do ID=`basename ${FILE} ".bw"`;

        SCRIPT=sitepro_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
	echo "${PCODE}/siteproBW -w ${FILE} -b ${INTERVAL} --pf-res=${RES} --span=${SPAN} --name=${PDATA}/${ID} --confid  --dir" >>${PDATA}/${SCRIPT};        
        nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
