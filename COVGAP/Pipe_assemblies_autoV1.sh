#!/bin/bash
#SBATCH --qos=30min
#SBATCH --mem=50G
#SBATCH --cpus-per-task=12
date
ml purge
ml Miniconda2/4.3.30

#to be put into another script
readarray -t array < samplelist.l # from COVGAP a samplelist is made, here we are reading it for the array job
export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}
echo "Processing file "$sample_id".fa, fasta taken from vcf/"$sample_id".vcf, clean vcf drafted in procvcf/"$sample_id"_proc.vcf,  annotated vcf drafted at anno_vcfs/"$sample_id"_anno.vcf"
####add the fasta clearing, division and vcf calling###
ml MAFFT/7.467-GCCcore-7.3.0-with-extensions
source activate /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/COVGAP3/Downstream/Env/Rload

#to be put into another script
echo "Processing file "$sample_id".fa, fasta taken from genomes/"$sample_id".fa, alignment drafted in alignments/"$sample_id".alignment, vcf drafted at vcfs/"$sample_id".vcf.
Starting alignment with mafft.."
#mkdir alignments/
mafft --thread 12 --reorder --keeplength --mapout --kimura 1 --addfragments ../genomes/"$sample_id".fa --auto /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/referenceNC.fasta > alignments/"$sample_id".alignment
echo "..alignment done. Now drafting the vcf file.."
#mkdir vcfs/
snp-sites -v -o vcfs/"$sample_id".vcf alignments/"$sample_id".alignment

####now preprocessing of vcf files#####
source deactivate
source activate /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Env/SnpEff4.5/ 
echo "preprocessing the vcf file.. " 
newhead=$(head -n 3 vcfs/"$sample_id".vcf | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}' OFS="\t")
newbody=$(tail -n +4 vcfs/"$sample_id".vcf | sed 's/1/NC_045512.2/1' | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}' OFS="\t")
#mkdir processed_vcfs/
echo "$newhead" > processed_vcfs/"$sample_id"_proc.vcf
echo "$newbody" >> processed_vcfs/"$sample_id"_proc.vcf
echo "..preprocess done, Now annotating the  vcf file.."
#mkdir annotated_vcfs/
snpEff -v NC_045512.2 processed_vcfs/"$sample_id"_proc.vcf > annotated_vcfs/"$sample_id"_anno.vcf
ml purge
ml R
echo "Now reporting the findings from the annotated file into a better shape, available at Reports/"$sample_id".variant.regions.annotated.primary.alleles.report.tab "
#mkdir Reports/
Rscript /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Genomevisual_assembliesV1.R "$sample_id" annotated_vcfs/ Reports/
echo "Sample $sample_id processed."
touch "$sample_id".flag
date
