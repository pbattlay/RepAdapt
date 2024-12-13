#!/bin/bash
# 1 CPU
# 10 Go

#SBATCH -J 06c.Local_Depth
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 10G
#SBATCH --time=0-6:00:00

# Global variables
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*{fasta,fa,fasta.gz,fa.gz} | xargs -n 1 basename)
ANNOTATION="clean.gff"
ALIGNEDFOLDER="06_bam_files"
METRICSFOLDER="99_metrics_merged"
SVFOLDER="97_Local_Depth"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
SCRIPTNAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp $SCRIPT $LOG_FOLDER/${TIMESTAMP}_${SCRIPTNAME}

mkdir $SVFOLDER

# Load needed modules
module load samtools
module load bedtools

# if this is the first array
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then

# make a BED file of 5000 bp windows from FASTA index
cat $GENOMEFOLDER/${GENOME}.fai | awk -v w=5000 '{chr = $1; chr_len = $2;
    for (start = 0; start < chr_len; start += w) {
        end = ((start + w) < chr_len ? (start + w) : chr_len);
        print chr "\t" start "\t" end;
    }
}' | sort -k1,1 -k2,2n > 02_info_files/windows.bed

# and a bed file of each gene
cat $GENOMEFOLDER/$ANNOTATION | awk '$3 == "gene" {print $1"\t"$4"\t"$5}' \
| sort -k1,1 -k2,2n > 02_info_files/genes.bed

fi

# Fetch filename from the array
file=$(cut -f2 02_info_files/datatable.txt | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d")
bamfile=${file}.realigned.bam

samtools depth -aa $ALIGNEDFOLDER/$bamfile > $SVFOLDER/${file}.depth

# gene depth analysis
echo \n">>> Computing depth of each gene for $file <<<"\n
cat $SVFOLDER/${file}.depth | awk '{print $1"\t"$2"\t"$2"\t"$3}' | bedtools map -a 02_info_files/genes.bed -b stdin -c 4 -o mean > $SVFOLDER/${file}-genes.bed

# window depth analysis
echo \n">>> Computing depth of each window for $file <<<"\n
cat $SVFOLDER/${file}.depth | awk '{print $1"\t"$2"\t"$2"\t"$3}' | bedtools map -a 02_info_files/windows.bed -b stdin -c 4 -o mean > $SVFOLDER/${file}-windows.bed

# overall genome depth
echo \n">>> Computing depth of whole genome for $file <<<"\n
cat $SVFOLDER/${file}.depth | awk '{sum += $3; count++} END {if (count > 0) print sum/count; else print "No data"}' > $SVFOLDER/${file}-wg.txt

echo " >>> Cleaning a bit...
"
rm -rf $SVFOLDER/${file}.depth
echo "
DONE! Check your files"



