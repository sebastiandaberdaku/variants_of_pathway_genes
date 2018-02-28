#!/bin/bash

# This script extracts the snp records in VCF format from the complete dbSNP Human database.
# The script uses multiprocessing to speed-up the process.
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    "$@" 
    printf '%.3d' $? >&3
    )&
}
task(){
	vcftools --gzvcf "00-All.vcf.gz" --snps "$1" --recode --recode-INFO-all --stdout | gzip -c > "$2/${2}_dbSNP/$(basename $ids _rsIDs.txt).vcf.gz"
}

# Uncomment the following line to download the current dbSNP Human database in VCF format.
#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz

N=12 # number of processes to be run in parallel
open_sem $N
for pathway in "insulin_signaling" "tryptophan_metabolism"
do
	for ids in "$pathway"/"$pathway"_dbSNP/*_rsIDs.txt
	do
    	run_with_lock task $ids $pathway 
	done
done