#!/bin/bash
#This pipeline was constructed on an orignal backbone of a variant calling for Influenza pipeline by Daniel Wuettrich, adapted and for COVID-19 by Tim Roloff and extended with new features by Alfredo Mari, 
#in case of troubles please contact Alfredo Mari at: alfredo.mari@unibas.ch
#Rational: this pipeline assumes the user is in possess of already demultiplexed Illumina reads obtained throught Tilling amplification as described in the artic protocol for nanopore
#from this, it calls the variants compared to th Wuhan reference, outputting Statistics , consensus sequences and plots describing the variant position, and incidence.
#Overall steps summary:
#-trims the reads from the adaptors by trimmomatic
#-aligns them to the reference with bwa
#-outputs the alignment in different formats for later pprocesses and cleans it up from duplicates
#-calls the variants with pilon
#-annotates them with SnpEff no filter to the variants is applied
#-calls the consensus fasta sequence by applying the variants and a coverage threshold of 750 reads per region
#-outputs a record of variants, their functions, whether they fall into the commonly used region for positive-negative test via qPCR, binding sites of Lopinavir or Chloroquin 
#-outputs coverage plot
#-outputs a trackplot including coverage, main COVId relevant pharmaceuticals binding sites (Lopinavir and Chloroquine, as well as COVID-19 antigen binding sites)
#for details on how the binding sites for Lopinavir and Chloroquine were called, please see the dir /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/build* and read the readme file

#this script is intened to run as a subscript of COV_GAP#.sh, please visit COV_GAP#.sh for more info, 
#to run it singularly -these passages are all embedded in COV_GAP#.sh:
#prepare a file called samplelist.l containing "mock" at the beginning, followed by the name of samples 
#sbatch --array=1-[total number of samples] --job-name= [JOB name] Pipe_Illumina_auto_V11.sh


#SBATCH --time=03:00:00
#SBATCH --qos=6hours
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#time recording
date
time1=$(date | awk '{print $4}')
sec1=`date +%s -d ${time1}`
#loading modules
echo "###### Starting preparatory operations.."
echo "###### Loading modules.."
ml purge
ml Miniconda2/4.3.30
#sourcing the environments
echo "###### Sourcing environments.."
source activate /scicore/home/egliadr/GROUP/Software/conda_environments/Variant_calling
export PATH="/scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality:$PATH"
export softwarepath="/scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality"
export _JAVA_OPTIONS="-Xmx10g"

echo "###### Reading samplename array: INPUT: samplelist.l"
readarray -t array < samplelist.l # from COVGAP a samplelist is made, here we are reading it for the array job

echo "###### Defining the reference as /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/referenceNC.fasta"
export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}
reference=/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/referenceNC.fasta

echo "###### Defining the reads. input: ./reads/"$sample_id"_R*.fastq.gz"
raw_r1=../reads/"$sample_id"_R1.fastq.gz
raw_r2=../reads/"$sample_id"_R2.fastq.gz

echo "###### Creating the trimmomatic dir: result/"$sample_id"/trimmomatic"
mkdir -p result/"$sample_id"/trimmomatic

################# Part1: QUALITY CONTROL ##################
#Software needed: Trimmomatic
#Dependencies needed: Primers.fa and adapters.fa, where the adapters and primers sequences are stored
#Summary:
#-length based filter, phred score filter >20, 
#-adapter trim, 
#-primer trim, 
#-a new length trim to trash the reds too much affected by the primer cut.
#-preparing the reference to be used: copied locally 
echo "###### Starting PART1: Quality Control: ###line67###"
#Step1:first trim
echo "###### Step1: Performing first trim. ###line69###
USED: trimmomatic, 
INPUT: 
-../reads/"$sample_id"_R*.fastq.gz, 
OUTPUT: 
-result/"$sample_id"/trimmomatic/r*.firsttrimmed.fastq.gz, 
-result/"$sample_id"/trimmomatic/read_trimm_info"
trimmomatic PE -threads "$SLURM_CPUS_PER_TASK" -phred33 "$raw_r1" "$raw_r2" result/"$sample_id"/trimmomatic/r1.firsttrimmed.fastq.gz result/"$sample_id"/trimmomatic/r1.firsttrimmed.not-paired.fastq.gz result/"$sample_id"/trimmomatic/r2.firsttrimmed.fastq.gz result/"$sample_id"/trimmomatic/r2.firsttrimmed.not-paired.fastq.gz SLIDINGWINDOW:4:20 MINLEN:125 2> result/"$sample_id"/trimmomatic/read_trimm_info
#Step2:adaptor trim
echo "###### Step2: Performing adaptor trim. ###line78###
USED: trimmomatic, 
INPUT: 
-result/"$sample_id"/trimmomatic/r*.firsttrimmed.fastq.gz, 
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/adapters.fa, 
OUTPUT: 
-result/"$sample_id"/trimmomatic/r*.no_adaptors.fastq.gz, 
-result/"$sample_id"/trimmomatic/adaptor_read_trimm_info"
trimmomatic PE -threads "$SLURM_CPUS_PER_TASK" -phred33 result/"$sample_id"/trimmomatic/r1.firsttrimmed.fastq.gz result/"$sample_id"/trimmomatic/r2.firsttrimmed.fastq.gz result/"$sample_id"/trimmomatic/r1.no_adaptors.fastq.gz result/"$sample_id"/trimmomatic/r1.no_adaptors.not-paired.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors.not-paired.fastq.gz ILLUMINACLIP:/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/adapters.fa:2:30:10 2> result/"$sample_id"/trimmomatic/adaptor_read_trimm_info
#Step3:primer trim
echo "###### Step3: Performing primer trim. ###line88###
USED; trimmomatic, 
INPUT: 
-result/"$sample_id"/trimmomatic/r*.no_adaptors.fastq.gz, 
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/primers.fa, 
OUTPUT: 
-result/"$sample_id"/trimmomatic/r*.no_adaptors_no_primers.fastq.gz, 
-result/"$sample_id"/trimmomatic/primer_read_trimm_info"
trimmomatic PE -threads "$SLURM_CPUS_PER_TASK" -phred33 result/"$sample_id"/trimmomatic/r1.no_adaptors.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors.fastq.gz result/"$sample_id"/trimmomatic/r1.no_adaptors_no_primers.fastq.gz result/"$sample_id"/trimmomatic/r1.no_adaptors_no_primers.not-paired.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors_no_primers.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors_no_primers.not-paired.fastq.gz ILLUMINACLIP:/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/primers.fa:2:30:10 2> result/"$sample_id"/trimmomatic/primer_read_trimm_info
#Step4:final length trim to eliminate shortened reads
echo "###### Step4: Performing last length trim. ###line98###
USED: trimmomatic, 
INPUT: 
-result/"$sample_id"/trimmomatic/r*.no_adaptors_no_primers.fastq.gz, 
OUTPUT: 
-result/"$sample_id"/trimmomatic/r*.no_adaptors_no_primers_trimmed.fastq.gz, 
-result/"$sample_id"/trimmomatic/no_adaptor_no_primer_quality_read_trimm_info"
trimmomatic PE -threads "$SLURM_CPUS_PER_TASK" -phred33 result/"$sample_id"/trimmomatic/r1.no_adaptors_no_primers.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors_no_primers.fastq.gz result/"$sample_id"/trimmomatic/r1.no_adaptors_no_primers_trimmed.fastq.gz result/"$sample_id"/trimmomatic/r1.no_adaptors_no_primers_trimmed_not-paired.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors_no_primers_trimmed.fastq.gz result/"$sample_id"/trimmomatic/r2.no_adaptors_no_primers_not-paired.fastq.gz SLIDINGWINDOW:4:20 MINLEN:125 2> result/"$sample_id"/trimmomatic/no_adaptor_no_primer_quality_read_trimm_info

echo "###### Defining working reads as result/"$sample_id"/trimmomatic/r*.no_adaptors_no_primers_trimmed.fastq.gz.. ###line107###"
r1=result/"$sample_id"/trimmomatic/r1.no_adaptors_no_primers_trimmed.fastq.gz
r2=result/"$sample_id"/trimmomatic/r2.no_adaptors_no_primers_trimmed.fastq.gz

#Step5:copying the reference into local directory
echo "###### Step5: Creating Mapping directory and copying the reference inside ###line112###"
mkdir -p result/"$sample_id"/Mapping/
cp "$reference" result/"$sample_id"/Mapping/reference.fa

##################### Part2: ALIGNMENT and READ MAPPING #######################
#Software needed: BWA, SAMTools, VariantBam
#dependencies needed: Json files containing mapping and unmapping rules 
#Summary: 
#-The reads are mapped against the chosen reference (NC_045512.2 -Wuhan1), 
#-further converted into bam format for later handling, and sorted. 
#-The alignments are then screened for mapping and unmapping reads using variantbam
#-Statistics are outputted concerning the ratio of mapped and unmapped reads in a dedicated file 
#-Resorting of the filtered bam
#-Selective downsampling of high coverage regions (>1000 reads per section) for graphical purpose only later on, to better identify droputs
#-Coverge stats are outputted as well
echo "###### Starting PART2: ALIGNMENT and READ MAPPING: ###line127###"
#Step1: Aligning the reads to the reference
echo "###### Step1: Performing alignment to the reference ###line129###
USED: bwa, 
INPUT: 
-result/"$sample_id"/Mapping/reference.fa, 
-result/"$sample_id"/trimmomatic/r*.no_adaptors_no_primers_trimmed.fastq.gz, 
OUTPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.sam"
bwa index result/"$sample_id"/Mapping/reference.fa
bwa mem -t "$SLURM_CPUS_PER_TASK" result/"$sample_id"/Mapping/reference.fa "$r1" "$r2" > result/"$sample_id"/Mapping/"$sample_id".alignment.sam
#Step2: Conversion from Sam into bam, indexing and sorting
echo "###### Step2: Converting SAM files into BAM, indexing and sorting, ###line139### 
USED: samtootls, 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.sam, 
OUTPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.bam, 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.bam"
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T result/"$sample_id"/Mapping/temp_sort -o result/"$sample_id"/Mapping/"$sample_id".alignment.bam result/"$sample_id"/Mapping/"$sample_id".alignment.sam
samtools rmdup result/"$sample_id"/Mapping/"$sample_id".alignment.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.bam
samtools index result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.bam
#Step3: Map filter, discriminating the mapping from the unmapping reads.
echo "###### Step3: Filtering just the mapped reads, ###line150### 
USED: variantbam, 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.bam, 
OUTPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sam, 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sam"
ml VariantBam #sourcing the module
variant result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.bam -r /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/map_rules.json > result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sam -v
variant result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.bam -r /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/unmap_rules.json >  result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sam -v
#Step4: Mapping stats.
echo "###### Step4: Creating the mapping stats. ###line161###
USED: samtools, bash, bc, 
INPUT:
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sam, 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sam, 
OUTPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.stats.tab "
countmap="$(samtools view -c result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sam)"
countunmap="$(samtools view -c result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sam)"
tot=$((countmap+countunmap))
percmap=$(bc <<< "scale=4;$countmap/$tot*100")
percunmap=$(bc <<< "scale=4;$countunmap/$tot*100")
echo "$sample_id" > result/"$sample_id"/Mapping/"$sample_id".alignment.stats.tab
echo "Total read count= $tot , 100 %" >> result/"$sample_id"/Mapping/"$sample_id".alignment.stats.tab
echo "Mapped read count= $countmap , $percmap %" >> result/"$sample_id"/Mapping/"$sample_id".alignment.stats.tab
echo "Unmapped read count= $countunmap , $percunmap %" >> result/"$sample_id"/Mapping/"$sample_id".alignment.stats.tab
#Step5: Resorting and indexing.
echo "###### Step5: Resorting and indexing of the mapping bams. ###line178### 
USED: samtools, 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.*mapped.reads.only.sam, 
OUTPUT:
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam"
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T result/"$sample_id"/Mapping/temp_sort result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sam -o result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T result/"$sample_id"/Mapping/temp_sort result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sam -o result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sorted.bam
samtools index result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam
samtools index result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sorted.bam
#Step6: Selective downsampling anf the alignments only for the regions displaying coverage > 1000 bp per section, sorting of it and indexing for later usage.
echo "###### Step6: Selective downsampling to 1000 reads max, will be used only for cosmetic reaons in the plots. ###line189###
USED: variantbam, samtools sort, samtools index, 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam, 
OUTPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted.bam"
variant result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam -m 1000 > result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sam -v
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T result/"$sample_id"/Mapping/temp_sort result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sam -o result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted.bam
samtools index result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted.bam
#Step7: Coverage stats are being outputted.
echo "###### Step7: Brewing coverage stats ###line199###
USED: samtools depth,awk 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam, 
OUTPUT: 
-result/"$sample_id"/Mapping/"$sample_id".average.coverage.tab "
echo "$sample_id" > result/"$sample_id"/Mapping/"$sample_id".average.coverage.tab
samtools depth -a result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> result/"$sample_id"/Mapping/"$sample_id".average.coverage.tab


#################### Part3: VARIANT CALL ####################
#Software needed: Pilon, vcflib, Java, snpEff
#dependencies needed: secondary_revealer1.sh, SnpEFF SARS genome repository 
#Summary: 
#-Pilon calls the variants applying no particular quality or information filter
#-For the primary variants: vcffilter individuates SNPs applying a cutoff of AF>0.7, DP>50, parallely it individuates the deletions and insertion with the DEL/INS filter, grep eleiminates the duplicates
#-For the seconday alleles: vcffilter creates a file with only AF>0 are collected -this step is done to speed up the process-
#-the bash script secondary_revealer1.sh loops through the AF>0 vcf and identifies the secondary alleles at a given position with an 0<AF<0.99, outputs in the "sample_id".secondary_alleles.vcf and report.tech 
#-finally snpeff annotates the primary and the secondary alleles.
echo "###### Starting Part3: VARIANT CALL ###line218###"

#Step1: variant call
echo "###### Step1 (variantcall): Create variantcall/ directory.."
mkdir -p result/"$sample_id"/variantcall/
echo "###### Calling variants, no filters at this stage. ###line223### 
USED: pilon, 
INPUT: 
-result/"$sample_id"/Mapping/reference.fa, 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".vcf"
pilon --threads "$SLURM_CPUS_PER_TASK" --genome result/"$sample_id"/Mapping/reference.fa --frags result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam --changes --variant --outdir result/"$sample_id"/variantcall/ --output "$sample_id"

#Step2: variant filtering to identify primary alleles
echo "###### Step2: variant filtering: purging modules and Loading vcflib ###line233###"
ml purge
ml vcflib
echo "###### Step2.1: Filtering the VCF, applied cutoff for AF= 0.7, for DP= 50. ###line236###
USED: vcffilter(vcflib), 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".vcf, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".depth.filtered.vcf"
vcffilter -f "SVTYPE = DEL & ! IMPRECISE" --or -f "SVTYPE = INS & ! IMPRECISE" --or -f "AF > 0.7 & DP > 50" result/"$sample_id"/variantcall/"$sample_id".vcf > result/"$sample_id"/variantcall/"$sample_id".depth.filtered.vcf
echo "###### Step2.2: Removing duplicates. ###line243###
USED: grep, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".depth.filtered.vcf, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".final.vcf"
grep -v "<DUP>" result/"$sample_id"/variantcall/"$sample_id".depth.filtered.vcf  > result/"$sample_id"/variantcall/"$sample_id".final.vcf
#Step3: variant filtering and identification of secondary alleles
echo "###### Step3: identification of secondary alleles: prefilter of all the AF > 0 -the following identifire uses a loop and given the vast majority of 0, otherwise would take too long-, ###line251###
USED: vcffilter, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".vcf, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".afM0.vcf"
vcffilter -f "AF > 0" result/"$sample_id"/variantcall/"$sample_id".vcf > result/"$sample_id"/variantcall/"$sample_id".afM0.vcf
echo "###### Step3.1: Picking the minority alleles, ###line258###
USED: /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/minority_revealer.sh, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".afM0.vcf, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".minority_alleles.vcf, 
-result/"$sample_id"/variantcall/"$sample_id".minority_alleles_report.tech"
/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/minority_revealer.sh result/"$sample_id"/variantcall/"$sample_id".afM0.vcf result/"$sample_id"/variantcall/"$sample_id"
#Step4: primary and secondary variant annotation, updated the snpEff version to the 4.5-COVID19
echo "###### Step4: variant annotation: purging modules and environments, ###line267###
USED: SnpEff version 4.5 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".final.vcf, 
-result/"$sample_id"/variantcall/"$sample_id".minority_alleles.vcf, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".annotated.variants.vcf, 
-result/"$sample_id"/variantcall/"$sample_id".annotated.minorityvariants.vcf"
ml purge 
ml Miniconda2/4.3.30
source deactivate
source activate /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Env/SnpEff4.5/
snpEff -v NC_045512.2 result/"$sample_id"/variantcall/"$sample_id".final.vcf > result/"$sample_id"/variantcall/"$sample_id".annotated.variants.vcf
snpEff -v NC_045512.2 result/"$sample_id"/variantcall/"$sample_id".minority_alleles.vcf > result/"$sample_id"/variantcall/"$sample_id".annotated.minorityvariants.vcf
source deactivate 

#################### Part4: CONSENSUS CALL ####################
#Software needed: bedtools,bcftools, seqtk,
#dependencies needed: none 
#Summary: 
#-First create the bedgraph files of the converage of both the subsampled and the native bam files, it will be used later for plotting
#-Creation of a bedfiles featuring the low coverage areas set as < 50
#-Subtraction from this bedfile of the position of the variants (they have anyway more than 50 reads to be called), this is done to include deletions, otherwiseregarded as regions with 0 coverage. -finalmask
#-Compression and indexing of the variant file
#-Creation of the consensus from the variant file, the reference, and the final mask featuring low coverage regions. Those regions will be called N
#-Calculation of the N percentage per sequence, this is later useful to tag the samples as High coverage (%N<10) or Low coverage (%N>10)
#-File tagging
#-Renaming of the header of the consensus sequence with the samplename
echo "###### Starting PART4: Consensus Call ###line295###"
ml purge
ml BEDTools 
ml BCFtools
#Step1: plotting the coverage:
echo "###### Step1: plotting the coverage, producing bedgraphs files first, ###line300### 
USED: bedtools, 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.*trimmed.sorted.bam, 
-OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".rd.*uptrim.bedgraph "
bedtools genomecov -ibam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam -bga | awk '{print "\t"$2"\t"$3"\t"$4}' |sed 's/^/NC_045512.2/g' > result/"$sample_id"/variantcall/"$sample_id".rd.nouptrim.bedgraph
bedtools genomecov -ibam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted.bam -bga | awk '{print "\t"$2"\t"$3"\t"$4}' |sed 's/^/NC_045512.2/g' > result/"$sample_id"/variantcall/"$sample_id".rd.1000uptrim.bedgraph
#Step2: identifying a first N mask for low coverage and then correcting it for the deletions
echo "###### Step2: identifying the regions to be masked as Ns, all the positions showing less than 50 reads, then overriding the variants on top of the mask (deletions are considered no coverage zones, this way we retain them) ###line309###
USED: bedtools, 
INPUT: 
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".lessthan50.bed, 
-result/"$sample_id"/variantcall/"$sample_id".finalmask.bed"
bedtools genomecov -bga -ibam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam | awk '$4<50' > result/"$sample_id"/variantcall/"$sample_id".lessthan50.bed
bedtools subtract -a result/"$sample_id"/variantcall/"$sample_id".lessthan50.bed -b result/"$sample_id"/variantcall/"$sample_id".final.vcf > result/"$sample_id"/variantcall/"$sample_id".finalmask.bed
#Step3: compression of variantfile -it is here and not before becuase of the incompatibility of BCF modules and vcflib modules
echo "###### Step3: Compression of files for further processing, ###line319###
USED: bgzip, tabix, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".final.vcf, 
-result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz, 
-result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz.tbi"
bgzip -c -i result/"$sample_id"/variantcall/"$sample_id".final.vcf > result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz
tabix -p vcf result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz -f > result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz.tbi
#Step4: consensus calling
echo "###### Step4: First consensus calling, ###line330###
USED: bcftools, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz, 
-result/"$sample_id"/Mapping/reference.fa, 
-result/"$sample_id"/variantcall/"$sample_id".finalmask.bed, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta"
bcftools consensus result/"$sample_id"/variantcall/"$sample_id".consvariants.vcf.gz -f result/"$sample_id"/Mapping/reference.fa -m result/"$sample_id"/variantcall/"$sample_id".finalmask.bed -o result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta
#Step5: Elaborating the stats on the Ns 
echo "#####  Step5: elaborating Nstats, crucial for tagging the samples, ###line340###
USED: seqtk, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".Nstats.tab"
ml seqtk/1.2-foss-2018b
echo "$sample_id" > result/"$sample_id"/variantcall/"$sample_id".Nstats.tab
seqtk comp result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta | awk '{Tot+=$2; N+=$9} END { print "Total length = ",Tot; print "Ns = ",N; print "Ns Perc = ", N/Tot*100,"%"}' >> result/"$sample_id"/variantcall/"$sample_id".Nstats.tab
#Step6: Sample tagging: classify the sequence as bad or good based on the percentage of Ns, <10% will be considered good, >10% will be considered trash, however, no sequence will be dropped
echo "###### Step6: sample tagging, labelling of HighCov and LowCov according to the 10% Ns threshold ###line350### 
USED: seqtk, awk,
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta,  
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".classification.tab,
-result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.*.sorted."$sampletag".bam"
nperc="$(seqtk comp result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta | awk '{Tot+=$2; N+=$9} END {print N/Tot*100}')"
echo "$nperc"
if (( $(echo "$nperc < 10" |bc -l) ));
then
echo "Nperc= $nperc, sample $sample_id will be considered as a Good sample, and labelled accordingly: HighCov"
sampletag="HighCov"
else
echo "Nperc= $nperc, sample $sample_id will be considered as a Bad sample, and labelled accordingly: LowCov"
sampletag="LowCov"
fi
echo "$sample_id = $sampletag with Ns: $nperc % of the genome" > result/"$sample_id"/variantcall/"$sample_id".classification.tab
##some renaming according to the labels
mv result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted."$sampletag".bam
mv result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sorted.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sorted."$sampletag".bam
mv result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted."$sampletag".bam

#Step7:change the name of the header of the consensus to the sample name and not the reference
echo "###### Step7: renaming the header of the reference with the sample name, handy for later analysis, ###line374###
USED: sed, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".consensus."$sampletag".fasta 
NOTE: this is the final consensus. the initial one is removed in this step"
sed "s/NC_045512.2/"$sample_id"_"$sampletag"/g" result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta > result/"$sample_id"/variantcall/"$sample_id".consensus."$sampletag".fasta
rm result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta


#################### Part5: VISUALISATION AND REPORT DRAFTING ####################
#Software needed: R
#dependencies needed: 
#-Environments and scripts: 
#---Seqinr and Rvisual (Envs), 
#---Ner.R, Covplotter.R, Genomevisual6.1.R
#-Files:
#---SARS-COV2 genome in gff format (GCA_009858895.3_ASM985889v3_genomic.gff in dependencies)
#---Genomic ranges for Lopinavir, Remdesivir, HC, and Bcells epitopes (*.ranges.txt in dependencies)
#Summary: 
#-The script Ner.R is creating a bed-like file listing the position of the Ns in the genome, later useful for the final Genomevisual script.
#-The script Covplotter is creating two coverage plots based on the bedgraphs of part4.
#The script Genomevisual is creating two reports: one for primary alleles and one for secondary alleles, basically adding to the annotated vcf file the information about potential farmaco-resistance.
#on top of this, it is creating a genome track plot -> a graphical representation of the report. Note that this graphics features ONLY the dominant alleles. no secondary alleles involved.
echo "###### Starting PART5: Visualisation and report drafting"
#Step1: bed-like file featuring the Npositions
echo "###### Step1: purging modules and loading Seqinr env, assessing the Npositions across the genomes, ###line401###
USED: /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Ner.R, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".consensus."$sampletag".fasta, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".NPosition.tab"
source deactivate
ml purge
ml Miniconda2/4.3.30
source activate /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Env/Seqinr/
Rscript /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Ner.R result/"$sample_id"/variantcall/"$sample_id".consensus."$sampletag".fasta result/"$sample_id"/variantcall/"$sample_id".NPosition.tab

#Step2: Coverage plots, report drafting and plots 
echo "###### Step2: Coverage plot, drafting trackplots and reports of variants, secondary variants. ###line414###
USED: /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Covplotter.R, 
INPUT: 
-result/"$sample_id"/variantcall/"$sample_id".rd.nouptrim.bedgraph, 
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".*uptrim.coverage.pdf"
source deactivate 
source activate /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Env/RVisual
Rscript /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Covplotter.R result/"$sample_id"/variantcall/"$sample_id".rd.nouptrim.bedgraph result/"$sample_id"/variantcall/"$sample_id".nouptrim.coverage.pdf "$sample_id"
Rscript /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Covplotter.R result/"$sample_id"/variantcall/"$sample_id".rd.1000uptrim.bedgraph result/"$sample_id"/variantcall/"$sample_id".1000uptrim.coverage.pdf "$sample_id"
echo "###### Report drafting and trackplot drafting, ###line348###
USED: /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Genomevisual6.4.1, 
INPUT:
-result/"$sample_id"/variantcall/"$sample_id".NPosition.tab
-result/"$sample_id"/variantcall/"$sample_id".annotated.variants.vcf
-result/"$sample_id"/variantcall/"$sample_id".annotated.minorityvariants.vcf 
-result/"$sample_id"/variantcall/"$sample_id".minority_alleles_report.tech
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/ranges_epitopes.txt
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Mproranges.txt
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/ranges_remdesivir1.txt
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/ranges_chloroquine.txt
-/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/GCA_009858895.3_ASM985889v3_genomic1.gff
OUTPUT: 
-result/"$sample_id"/variantcall/"$sample_id".variant.regions.annotated.primary.alleles.report.*tag*.tab
-result/"$sample_id"/variantcall/"$sample_id".variant.regions.annotated.minority.alleles.report.*tag*.tab
-result/"$sample_id"/variantcall/"$sample_id".nocoverage.trackplot.*tag*.pdf
-result/"$sample_id"/variantcall/"$sample_id".complete.trackplot.*tag*.pdf"
Rscript /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Genomevisual6.4.1.R "$sample_id" result/"$sample_id"/variantcall/ "$sampletag"


#### FINAL  TOUCH ####
#these last commands serve the overscript, which counts every three minutes if the number of flagfiles is the same as the sample scheduled to be run. when yes, it starts summarising plots and reports.
#it also includes the date information, for time run assessment. 
echo "Final closing features: date and time, and creation of flag files. ###line447###"
touch final_"$sample_id".flag
date
time2=$(date | awk '{print $4}')
sec2=`date +%s -d ${time2}`
diffsec=`expr ${sec2} - ${sec1}`
echo Start ${time1}
echo Finish ${time2}
echo Took ${diffsec} seconds
echo Total time elapsed `date +%H:%M:%S -ud @${diffsec}`
