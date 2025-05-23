#!/bin/bash
# 1 CPU
# 30 Go

#SBATCH -J "07.mpileup"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00-72:00:00

cd $SLURM_SUBMIT_DIR

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp $SCRIPT $LOG_FOLDER/${TIMESTAMP}_${NAME}

begin=`date +%s`

# Load needed modules
module load bcftools/1.11

# Global variables
INFO="02_info_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*{fasta,fa,fasta.gz,fa.gz} | xargs -n 1 basename)
INDGENOME=$GENOMEFOLDER/${GENOME}.fai
VCF="./07_raw_VCFs"
# POP="02_info_files/popmap.txt"
BAM="02_info_files/bammap.txt"
PLD="02_info_files/ploidymap.txt"
REGION_FILE=$(ls 02_info_files/all_scafs*pos | sed "${SLURM_ARRAY_TASK_ID}q;d")

    for scaf in $(cut -f1 $REGION_FILE)
    do
    echo ">>> Genotyping scaffold $scaf"
    done

    # Action!
    parallel -j8 "bcftools mpileup -Ou -f $GENOMEFOLDER/$GENOME --bam-list $BAM -q 5 -r {} -I -a FMT/AD,FMT/DP | bcftools call -S $PLD -G - -f GQ -mv -Ov > $VCF/${DATASET}_{}.vcf" :::: $REGION_FILE

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
