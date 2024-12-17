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

### combine windowed depth analysis results
# get just the second (depth) column of each output file
cut -f2 02_info_files/datatable.txt | sort | uniq | sed 's/^/97_Local_Depth\//' | sed 's/$/-windows.sorted.tsv/' | while read tsv; do cut -f2 $tsv > ${tsv/-windows.sorted/-windows.sorted-depthcol} ; done

# combine both columns of the first output file and add the second column of all other files
paste $(cut -f2 02_info_files/datatable.txt | sort | uniq | sed 's/^/97_Local_Depth\//' | sed 's/$/-windows.sorted.tsv/' | head -n 1) $(cut -f2 02_info_files/datatable.txt | sort | uniq | sed 's/^/97_Local_Depth\//' | sed 's/$/-windows.sorted-depthcol.tsv/' | tail -n +2) > 97_Local_Depth/combined-windows.temp

# add header
cat 97_Local_Depth/depthheader.txt 97_Local_Depth/combined-windows.temp > 09_final_vcf/${DATASET}_combined_windows.tsv

### combine gene depth analysis results
# get just the second (depth) column of each output file
cut -f2 02_info_files/datatable.txt | sort | uniq | sed 's/^/97_Local_Depth\//' | sed 's/$/-genes.sorted.tsv/' | while read tsv; do cut -f2 $tsv > ${tsv/-genes/-genes-depthcol} ; done

# combine both columns of the first output file and add the second column of all other files
paste $(cut -f2 02_info_files/datatable.txt | sort | uniq | sed 's/^/97_Local_Depth\//' | sed 's/$/-genes.sorted.tsv/' | head -n 1) $(cut -f2 02_info_files/datatable.txt | sort | uniq | sed 's/^/97_Local_Depth\//' | sed 's/$/-genes.sorted-depthcol.tsv/' | tail -n +2) > 97_Local_Depth/combined-genes.temp

# add header
cat 97_Local_Depth/depthheader.txt 97_Local_Depth/combined-genes.temp > 09_final_vcf/${DATASET}_combined_genes.tsv

### make a table of whole-genome depths from samtools
cut -f2 02_info_files/datatable.txt | sort | uniq | while read file; do echo -e $file"\t"$(cat 97_Local_Depth/$file-wg.txt); done > 09_final_vcf/${DATASET}_combined_wg.tsv

echo "
DONE! Check your files"

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed s
