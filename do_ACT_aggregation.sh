#!/bin/bash

#10/10/2013
#compute aggregated signal for peaks around some genomic feature
#uses ACT from the gerstein lab
#http://info.gersteinlab.org/ACT_Tool
#The aggregation script agg-py takes values from multiple points on a single genomic signal track and creates an average signal profile around a set of anchor points, such as Transcription Start Sites (TSS's). 

#parameters
#ACT.py [--nbins=#] [--mbins=#] [--radius=#] [--median] [--mean] [--region] [--point] [--annotationparser=<name>] [--signalparser=<name>] [--mingenelen=#] [--output=<file>] annotationFile signalFile(s)
#eg python ../ACT.py --mean --radius=1000 --nbins=40 --output=testmean.out ch22_genome.pos PolII.chr22.sgr

#where
#annotationFile - tab-delimited file containing several start and end positions, chromosome locations, and strand annotations, such as a bed file
#signal file - tab-delimited file containing signals and positions (sgr, converted from wig for instance)

#--nbins. Specifies how many bins will be in each flanking region (for a total of 2xn bins). For jobs that do not use scaled bins in the region between start and stop sites, there will be 2xn bins: n bins on each side of the start site.

#--region. Optional tag. If included, the program will aggregate "radius" base pairs upstream of the start site, downstream of the stop site, AND the "region" between the start and stop site. The way it will do this is to create "mbins" number of bins between the start and stop site for each segment in the annotations file, and scale the bins based on the size of the segment. Note that, as a result, mbins must be smaller than the smallest segment in the annotation file.

#--mbins. Optional field to be defined if stop and start sites (regions) are to be considered in the aggregation process. Defines the number of bins assigned to the region between start and stop sites (for a total number of m+2n bins).

#--radius. Specifies number of base pairs the program should analyze in both directions from either the start site or the start/stop pair, depending on the option specified. (in the output, the program will present radius base pairs divided into nbins upstream and downstream of start sites or upstream of start sites and downstream of stop sites, depending on if --regions is used or not)

#--mingenelen: Tells the program not to consider annotations for start/stop pairs that are shorter than a certain distance.

#--mean. The program will automatically ensure that only the median signal of an individual gene (or number of probes, if the density option is selected) contributes to the final "averaging" calculation in each bin to avoid bias against shorter genes. If use mean is selected, the program will take the mean signal across all genes per bin and report the final result. Otherwise, it will use the median signal (of the median signal of each gene). 

#annotation file looks like this
#chr7    115444712       115491487       +
#chr7    115468538       115491487       +
#chr7    115733786       115740164       +
#chr7    115733786       115740126       +
#chr7    115759067       115793292       +
#chr7    115933089       116030129       +
#chr7    116008875       116030129       +

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <DATA> <NBINS> <MBINS> <RADIUS>"
        echo "DATA - path to the directory containing the signal file(s) to annotate";
	echo "NBINS - how many bins on each side of the start site";
	echo "MBINS - to be defined if stop and start sites (regions) are to be considered in the aggregation process";
	echo "RADIUS - number of base pairs the program should analyze in both directions from either the start site or the start/stop pair";
        exit;
fi

PCODE="/net/isi-scratch/giuseppe/tools/Agg";
PDATA=$1;
NBINS=$2;
MBINS=$3;
RADIUS=$4;
PANNO="/net/isi-scratch/giuseppe/indexes/Hsap/Ensembl73/ensembl73_known_protein_coding_genes_UCSC_sorted.bed";


for FILE in ${PDATA}/*.srg;
        do
        BASEFILE=`basename ${FILE} ".srg"`;
        SCRIPT=ACT_${BASEFILE}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        #echo "python ${PCODE}/ACT.py --nbins=${NBINS} --mbins=${MBINS} --radius=${RADIUS} --region --output=${PDATA}/ACT_${BASEFILE}.out ${PANNO} ${FILE}" >>${PDATA}/${SCRIPT};
	echo "python ${PCODE}/ACT.py --nbins=${NBINS}  --radius=${RADIUS} --output=${PDATA}/ACT_${BASEFILE}.out ${PANNO} ${FILE}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/ACT_${BASEFILE}.err -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
