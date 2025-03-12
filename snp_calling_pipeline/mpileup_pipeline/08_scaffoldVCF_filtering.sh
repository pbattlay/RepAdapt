#!/bin/bash

#SBATCH --job-name="08.FiltVCF"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00-12:00:00

module load vcftools bcftools
module load StdEnv/2020 intel/2020.1.217 tabix/0.2.6

cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
echo "STARTING AT $TIMESTAMP"
echo $SCRIPT
cp $SCRIPT $LOG_FOLDER/${TIMESTAMP}_${NAME}
begin=`date +%s`

# Variables
VCF="07_raw_VCFs"
FILTVCF="08_filtered_VCFs"

# Pull from the array...
REGION_FILE=$(ls 02_info_files/all_scafs*pos | sed "${SLURM_ARRAY_TASK_ID}q;d")

    echo "
    >>> Filtering through BCFtools first!
    "
    parallel -j8 "bcftools filter -e 'AC=AN || MQ < 30' $VCF/${DATASET}_{}.vcf -Oz > $FILTVCF/${DATASET}_{}_filtered.vcf.gz" :::: $REGION_FILE

    >>> Preparation for concatenation of VCF files
    "
    parallel -j8 "tabix -p vcf $FILTVCF/${DATASET}_{}_filtered.vcf.gz" :::: $REGION_FILE

echo "
>>> Cleaning a bit...
"
for scaf in $(cat $REGION_FILE)
do
  rm $FILTVCF/${DATASET}_${scaf}_filteredTmp.vcf.gz
done

echo "
DONE! Check you files"

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
