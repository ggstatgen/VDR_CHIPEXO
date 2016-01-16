#!/usr/bin/Rscript


#######################################################
## Script for performing plots on fastqscreen outputs
#######################################################


# script assumes fastq-screen output will be in the same format as the following example:


#   #Fastq_screen version: 0.4.1
#   Library	%Unmapped	%One_hit_one_library	%Multiple_hits_one_library	%One_hit_multiple_libraries	%Multiple_hits_multiple_libraries
#   Human	99.54	0.32	0.06	0.03	0.05
#   Mouse	18.51	66.31	14.68	0.42	0.08
#   Rat	99.56	0.00	0.00	0.42	0.02
#   Lizard	99.96	0.00	0.00	0.01	0.03
#   Ecoli	99.99	0.01	0.00	0.00	0.00
#   Yeast	100.00	0.00	0.00	0.00	0.00
#   
#   %Hit_no_libraries: 18.12




# User-provided information
#############################

## Input Directory
# this should be a folder containing all of your fastqscreen results
# this script will plot for every file ending "_screen.txt" in the folder

indir <- "/net/isi-scratch/harryc/fastqscreen"

## Output File
# output graphs will be placed in the input directory, in a single pdf, with the filename specified by the user
# (the below filename should end with ".pdf")

outfile <- "fastqscreen_plots.pdf"






#############################################################
# script
##########


# changes working directory and sets infiles
setwd(indir)
infiles <- Sys.glob('*_screen.txt')


# loops through infiles, plotting to pdf
pdf(outfile,width=14)
for(infile in infiles){
    
    # reads in data
    data <- read.delim(infile,stringsAsFactors=F,header=F)
    
    # removes first and last row 
    # (these should only contain a fastqscreen version number and a summary of hit libraries respectively)
    table <- data[ 3:(length(data[,1])-1) , -1]
    row.names(table) <- data[ 3:(length(data[,1])-1) , 1 ]
    colnames(table) <- data[ 2 , -1 ]
    
    # takes hits summary
    nohits <- as.numeric( strsplit( as.character(data[length(data[,1]),])[1] , ":" )[[1]][2] )
    
    # removes "%unmapped" from data and adds "No Hits"
    tmp_data <- t(as.matrix( table[,which(colnames(table)!="%Unmapped")] ))
    tmp_data2 <- cbind( tmp_data , c(nohits,rep("0.00",length(tmp_data[,1])-1)) )
    colnames(tmp_data2)[length(colnames(tmp_data2))] <- "No Hits"
    
    # plot
    cols <- c("blue","darkblue","red","darkred")
    title <- gsub("_screen.txt","",infile)
    par(mar=c(3.1, 5.1, 4.1, 2.1))
    layout(rbind(1,2), heights=c(9,1))
    
    barplot(tmp_data2,col=cols,ylim=c(0,100),ylab="% Mapped",main=title)
    
    legend_names <- c("One hit / one library", "Multiple hits / one library", "One hit / multiple libraries", "Multiple hits / multiple libraries")
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend('center',legend_names,col=cols,ncol=4,bty ="n", cex=0.8, pch=15)
    
    
}
dev.off()




