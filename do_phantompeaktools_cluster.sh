#!/bin/bash

#script to run modified spp for distance estimation
#see https://code.google.com/p/phantompeakqualtools/

#typical run
#Rscript run_spp.R [options]

#(1) Determine strand cross-correlation peak / predominant fragment length OR print out quality measures
#Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
#You can run the program on multiple datasets in parallel and append all the quality information to the same <outFile> for a summary analysis.

#this is how you'd use it to call actual peaks
#Rscript run_spp.R -c data.bam -i input.bam -p 15 -fdr=0.05 -savn -savr -savd -savp
#or maybe
#Rscript run_spp.R -c=chipSampleRep1.tagAlign.gz -i=controlSampleRep0.tagAlign.gz -npeak=300000 -odir=/peaks/reps -savr -savp -rf -out=/stats/phantomPeakStatsReps.tab

#Usage: Rscript run_spp_nodups.R <options>
#MANDATORY ARGUMENTS
#-c=<ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz) 
#MANDATORY ARGUMENTS FOR PEAK CALLING
#-i=<Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
#OPTIONAL ARGUMENTS
#-s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-500:5:1500
#-speak=<strPeak>, user-defined cross-correlation peak strandshift
#-x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10)
#-p=<nodes> , number of parallel processing nodes, default=0
#-fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling
#-npeak=<numPeaks>, threshold on number of peaks to call
#-tmpdir=<tempdir> , Temporary directory (if not specified R function tempdir() is used)
#-filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name
#OUTPUT ARGUMENTS
#-odir=<outputDirectory> name of output directory (If not set same as ChIP file directory is used)
#-savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)
#-savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width peaks with regions of enrichment)
#-savd=<rdatafile> OR -savd, save Rdata file
#-savp=<plotdatafile> OR -savp, save cross-correlation plot
#-out=<resultfile>, append peakshift/phantomPeak results to a file
#     format:Filename<tab>numReads<tab>estFragLen<tab>corr_estFragLen<tab>PhantomPeak<tab>corr_phantomPeak<tab>argmin_corr<tab>min_corr<tab>Normalized SCC (NSC)<tab>Relative SCC (RSC)<tab>QualityTag)
#-rf, if plot or rdata or narrowPeak file exists replace it. If not used then the run is aborted if the plot or Rdata or narrowPeak file exists
#-clean, if present will remove the original chip and control files after reading them in. CAUTION: Use only if the script calling run_spp.R is creating temporary files


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
#location of R script, change if needed
#PCODE="/net/isi-scratch/giuseppe/tools/spp_package";
PCODE="/net/isi-scratch/giuseppe/tools/phantompeakqualtools";
#PRSCRIPT="/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin"
#PRSCRIPT="/net/isi-scratch/giuseppe/tools/R-3.0.1/bin"
PRSCRIPT="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin/"


for FILE in ${PDATA}/*.bam;
	#change this depending on the data
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        do 
        #ID=`echo ${FILE} | egrep -o "38380_[0-9]*"`;        
        
        ID=`basename ${FILE} ".bam"`;
        
        SCRIPT=ppqtools_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};


	#change required options here:
	#-s=-500:1:500 -x=35:45
        echo "${PRSCRIPT}/Rscript ${PCODE}/run_spp.R -c=${FILE} -savp -s=-1500:1:1500 -x=35:45 -out=${PDATA}/ppqtools_${ID}_StatsReps.tab" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/ppqtools_${ID}.err -o ${PDATA}/ppqtools_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        #nice -5 qsub -v "BASH_ENV=~/.bashrc" -e ${PDATA}/ppqtools_${ID}.err -o ${PDATA}/ppqtools_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
