#!/bin/bash

#SBATCH --job-name="09.concatVCF"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00-12:00:00

module load vcftools bcftools samtools

cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
echo $SCRIPT
cp $SCRIPT ${LOG_FOLDER}/${TIMESTAMP}_${NAME}

# Variables
FILTVCF="08_filtered_VCFs"
CONCATVCF="09_final_vcf"

begin=`date +%s`

# Concatenate all the scaffold-VCF files into one global VCF file
vcf-concat $(ls -1 $FILTVCF/*_filtered.vcf | perl -pe 's/\n/ /g') > ${CONCATVCF}/${DATASET}_full_concatened.vcf && bgzip ${CONCATVCF}/${DATASET}_full_concatened.vcf
#vcf-concat $(ls -1 $FILTVCF/*_filtered.imputed.vcf.gz | perl -pe 's/\n/ /g') > ${CONCATVCF}/${DATASET}_full_concatened_imputed.vcf && bgzip ${CONCATVCF}/${DATASET}_full_concatened_imputed.vcf

# Add final maf filtering here...
bcftools view --min-af 0.01:minor ${CONCATVCF}/${DATASET}_full_concatened.vcf.gz -Oz -o ${CONCATVCF}/${DATASET}_full_concatened_maf01.vcf.gz
bcftools view --min-af 0.05:minor ${CONCATVCF}/${DATASET}_full_concatened.vcf.gz -Oz -o ${CONCATVCF}/${DATASET}_full_concatened_maf05.vcf.gz
#bcftools view --min-af 0.05:minor ${CONCATVCF}/${DATASET}_full_concatened_imputed.vcf.gz -Oz -o ${CONCATVCF}/${DATASET}_full_concatened_imputed_maf05.vcf.gz

###### combine windowed analysis merging

# make a header row
echo -e "location\t$(cut -f2 02_info_files/datatable.txt | sort | uniq | paste -s -d '\t')" > 97_Local_Depth/depthheader.txt

# get a list of sample names
cut -f2 02_info_files/datatable.txt | sort | uniq > 02_info_files/samples.txt

### combine windowed depth analysis results
# get just the second (depth) column of each output file
while read samp; do cut -f2 97_Local_Depth/${samp}-windows.sorted.tsv > 97_Local_Depth/${samp}-windows.sorted.depthcol ; done < 02_info_files/samples.txt

# combine both columns of the first output file and add the second column of all other files
paste $(sed 's/^/97_Local_Depth\//' 02_info_files/samples.txt | sed 's/$/-windows.sorted.tsv/' | head -n 1) $(sed 's/^/97_Local_Depth\//' 02_info_files/samples.txt | sed 's/$/-windows.sorted.depthcol/' | tail -n +2) > 97_Local_Depth/combined-windows.temp

# add header
cat 97_Local_Depth/depthheader.txt 97_Local_Depth/combined-windows.temp > 09_final_vcf/${DATASET}_combined_windows.tsv

### combine gene depth analysis results
# get just the second (depth) column of each output file
while read samp; do cut -f2 97_Local_Depth/${samp}-genes.sorted.tsv > 97_Local_Depth/${samp}-genes.sorted.depthcol ; done < 02_info_files/samples.txt

# combine both columns of the first output file and add the second column of all other files
paste $(sed 's/^/97_Local_Depth\//' 02_info_files/samples.txt | sed 's/$/-genes.sorted.tsv/' | head -n 1) $(sed 's/^/97_Local_Depth\//' 02_info_files/samples.txt | sed 's/$/-genes.sorted.depthcol/' | tail -n +2) > 97_Local_Depth/combined-genes.temp

# add header
cat 97_Local_Depth/depthheader.txt 97_Local_Depth/combined-genes.temp > 09_final_vcf/${DATASET}_combined_genes.tsv

### make a table of whole-genome depths from samtools
while read samp; do echo -e $samp"\t"$(cat 97_Local_Depth/$samp-wg.txt); done < 02_info_files/samples.txt > 09_final_vcf/${DATASET}_combined_wg.tsv

echo "
DONE! Check your files"

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed s
