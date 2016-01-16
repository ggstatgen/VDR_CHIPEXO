#!/bin/bash
# Converts pairwise IDR peak overlap output to narrowPeak
#EXMAPLE (Giuseppe)
#"chr1" "start1" "stop1" "sig.value1" "chr2" "start2" "stop2" "sig.value2" "idr.local" "IDR"
#"1" "chr1" 0 565729 2311.39898 "chr1" 0 565729 2083.05609 3.74912644350855e-05 3.74912644350855e-05
#"2" "chr12" 0 92464 8.23536 "chr12" 0 95470 12.67976 1 0.252690284525126
#compares column 3 and 7, gets the longer value (wider peak)
#compares column 4 and 8, gets the bigger value (taller peak)


if [[ "$#" -lt 1 ]]
    then
    echo 'Converts pairwise IDR peak overlap output to narrowPeak' 1>&2
    echo "USAGE: $(basename $0) [idrOverlapFile] [oDir]" 1>&2
    echo '[idrOverlapFile]: overlap output file from pairwise IDR analysis' 1>&2
    echo '[oDir]: output directory' 1>&2
    exit 1
fi

# overlap file
ovFile=$1
if [[ ! -e ${ovFile} ]]
    then
    echo "ERROR:${ovFile} does not exist" 1>&2
    exit 1
fi

# Output directory
oDir=$(dirname ${ovFile})
[[ $# -gt 1 ]] && oDir=$2
if [[ ! -d ${oDir} ]]
    then
    mkdir ${oDir}
fi
oDir=$(echo ${oDir} | sed -r 's:/$::g')

# Create output file
oFile="${oDir}/$(basename ${ovFile} .gz).npk"
if grep -q -E '\.gz$' ${ovFile}
    then
    #zcat ${ovFile} | sed 1d | sed -r 's/"//g' | sort -k11g,11g | awk '{if ($3 <=$7) st=$3 ; else st=$7 ; if ($4 >= $8) sto=$4 ; else sto=$8 ; printf "%s\t%d\t%d\t%d\t%s\t.\t%s\t%f\t%f\n",$2,st,sto,NR,$5,$9,-log($10)/log(10),-log($11)/log(10)}' | gzip -c > ${oFile}
    zcat ${ovFile} | sed 1d | sed -r 's/"//g' | sort -k11g,11g | awk '{if ($3 <=$7) st=$3 ; else st=$7 ; if ($4 >= $8) sto=$4 ; else sto=$8 ; printf "%s\t%d\t%d\t%d\t%s\t.\t%s\t%f\t%f\n",$2,st,sto,NR,$5,$9,-log($10)/log(10),-log($11)/log(10)}' > ${oFile}
else
    #sed 1d ${ovFile} | sed -r 's/"//g' | sort -k11g,11g | awk '{if ($3 <=$7) st=$3 ; else st=$7 ; if ($4 >= $8) sto=$4 ; else sto=$8 ; printf "%s\t%d\t%d\t%d\t%s\t.\t%s\t%f\t%f\n",$2,st,sto,NR,$5,$9,-log($10)/log(10),-log($11)/log(10)}' | gzip -c > ${oFile}
    sed 1d ${ovFile} | sed -r 's/"//g' | sort -k11g,11g | awk '{if ($3 <=$7) st=$3 ; else st=$7 ; if ($4 >= $8) sto=$4 ; else sto=$8 ; printf "%s\t%d\t%d\t%d\t%s\t.\t%s\t%f\t%f\n",$2,st,sto,NR,$5,$9,-log($10)/log(10),-log($11)/log(10)}' > ${oFile}
fi
