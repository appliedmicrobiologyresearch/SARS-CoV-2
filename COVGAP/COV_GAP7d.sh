#!/bin/bash
#SBATCH --time=04:30:00
#SBATCH --qos=6hours
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#time tracking
date
time1=$(date | awk '{print $4}')
sec1=`date +%s -d ${time1}`
#choose whether we are dealing with assemblies or reads,
run=$1
letter=$2
frags=$3
dir=analysis_"$run""$letter"
if [ "$frags" == "GENOMES"  ];
	then
	echo "##################################################### COVGAP: COVID-19 Genome Analysis Pipeline #########################################################
Alfredo Mari -University of Basel- alfredo.mari@unibas.ch
The mode GENOMES was chosen, fasta files will be treated as completed assempled genomes. NOTE: with the GENOMES mode certain features cannot be established:
-coverage, minority alleles for example. If you wish to retrieve this information you should access the raw reads from the original source.  
   Running command was:
sbatch --output=Autorun_RUNNAME_REPNUMB.out --job-name=COV-RUNNAME /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/COV_GAP7.sh STUDYNUMBER REPNUMBER
   Essential steps:
   aligns the sequences provided and produces an alignmet with mafft
-finds the variants with snpsites
-annotates the variants --> snpEff
-produces a report whether these variant affect drug related pockets: Remdesivir, Lopinavir, Chloroquine
   Output list:
-alignments --> alignments/
-variant files --> vcfs/, processed_vcfs/
-annotated variants --> annotated_variants/
-report of the variants position, annotation, potential affection of drug binding site: available Lopinavir, Chloroquine, Remdesivir -> Reports/"
	mkdir "$dir"
	cd "$dir"
	mkdir alignments/
	mkdir vcfs/
	mkdir processed_vcfs/
	mkdir annotated_vcfs/
	mkdir Reports/	
	echo "Creating directories and samplelists.."
        echo "mock" > samplelist.l
        rm ../genomes/*map       #previous alignments may have created a map file that messes up the array number count, we remove them preventively
	ls ../genomes/ | sed 's/.fa//g' >> samplelist.l
        echo "..done."
        echo "Defining the number of subjobs to run.."
        arraynumber="$(wc -l samplelist.l | awk '{print $1-1}')"
        echo "Done: $arraynumber"
        echo "Copying the latest version of the script from Illumina/Softwares/Pipe_assemblies_auto_V1.sh  in the working dir.."
        cp /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Pipe_assemblies_autoV1.sh .
        echo "..done."
        echo "The Following command has been launched:"
        echo "sbatch --array=1-"$arraynumber" --job-name= "$run" Pipe_assemblies_autoV1.sh "
        sbatch --array=1-"$arraynumber" --job-name="$run" Pipe_assemblies_autoV1.sh
        echo "Wait for the child jobs to be ready: it may take up to 30 min"
	
	 while true; do
          files=("./"*.flag)
          length=${#files[@]}
          rlength=$(($length-1))

          if [ "$length" != "$arraynumber" ];
          then
                echo "Number of jobs to be completed: $arraynumber .."
                echo "Number of jobs completed: $length"
                echo "..still waiting for completion.."
                sleep 180 ## checks the directory every 3 min to see whether the jobs are done
          else
		cat Reports/* > All."$run".primary.alleles.report.tab
		echo "..done! for doublecheck go to  $dir/Outdir/slurm_*.out"
                break
	  fi
	  done
	  rm *flag
	  mkdir Outdir/
	  mv slurm* Outdir/
	else
	#modules
	ml purge
	ml poppler
	ml Miniconda2/4.3.30
	ml MUSCLE
	source activate /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Env/dt
	mkdir "$dir"
	cd "$dir"
	echo "##################################################### COVGAP: COVID-19 Genome Analysis Pipeline #########################################################
Alfredo Mari -University of Basel- alfredo.mari@unibas.ch
The mode READS was chosen, fragments will be treated as raw Illumina reads.
  Running command was:
sbatch --output=Autorun_RUNNAME_REPNUMB.out --job-name=COV-RUNNAME /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/COV_GAP7.sh RUNNAME REPNUMBER
  Essential steps:
-trims the illumina adaptors and the V3 primers out, trimming reads by a sliding window of 4 bp in which the phred score is less than 20--> trimmomatic
-aligns the reads to the reference MN908947.3 --> bwa
-polishes the alignment (filters only mapped reads, removes duplicates, downsamples to 1000 reads the regions with too high coverage for portability) --> Variantbam, samtools
-calls the variants applying a depth support of at least 70% of the reads to call the SNP (calculated on the average depth)--> pilon
-annotates the variants --> snpEff
-calls the consensus masking as Ns any region covered with less than 50 reads --> bcftools
-produces coverage and trackplots, as well as reports indicating where the snps fall na dif they are possible candidate to antiviral resistance --> R
  Output list:
-consensus sequences including the variants, Ns are called for every region below 50 reads in coverage -> Sequences/
-report of the variants position, annotation, potential affection of drug binding site: available Lopinavir, Chloroquine, Remdesivir -> Reports/
-trackplots -with and without coverage displayed- of the report details in a compact pdf -> Reports/
-original coverage plot, prior to the downsampling at 1000 reads per high coverage region -> Reports/
-bam and sam files containing the alignments for maped, unmapped and trimmed reads -> Alignments/ "
	echo "Creating directories and samplelists.."
	echo "mock" > samplelist.l
	ls ../reads/ | grep 'R1' | sed 's/_R1.fastq.gz//g' >> samplelist.l
	echo "..done."
	echo "Defining the number of subjobs to run.."
	arraynumber="$(wc -l samplelist.l | awk '{print $1-1}')"
	echo "Done: $arraynumber"
	echo "Copying the latest version of the script from Illumina/Softwares/Pipe_Illumina_auto_V10.2.sh  in the working dir.."
	cp /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/Pipe_Illumina_auto_V10.3.sh .
	echo "..done."
	echo "The Following command has been launched:"
	echo "sbatch --array=1-"$arraynumber" --job-name= "$run" Pipe_Illumina_auto_V10.3.sh "
	sbatch --array=1-"$arraynumber" --job-name="$run" Pipe_Illumina_auto_V10.3.sh

	echo "Wait for the child jobs to be ready: it may take up to 30 min"

	while true; do
	  files=("./"*.flag) 
	  length=${#files[@]}
	  rlength=$(($length-1))  

	  if [ "$length" != "$arraynumber" ];
	  then
	       	echo "Number of jobs to be completed: $arraynumber .."
		echo "Number of jobs completed: $length"
		echo "..still waiting for completion.."
		sleep 180 ## checks the directory every 3 min to see whether the jobs are done
	  else
		echo "Number of jobs to be completed: $arraynumber .."
	  	echo "Number of jobs completed: $length"
		echo "All the children jobs completed, now wrapping up the final files.."
	        echo "Report drafted at: $dir/Reports/All."$run".report.tab"
		echo "N statistics drafted at: $dir/Reports/All."$run".Nstats.tab"
		echo "Alignment statistics -mapped/unmapped- drafted at: $dir/Reports/All."$run".alignmentStats.tab"
		echo "Trackplot summary available at: $dir/Reports/All."$run".complete.trackplots.pdf"
		echo "Trackplot summary -without coverage- at: $dir/Reports/All."$run".nocov.trackplots.pdf"
		echo "Original coverage plot -before regional downsampling  drafted at: $dir/Reports/All."$run".pre-trim.coverage.pdf"
		mkdir -p Reports/
		cat result/*/variantcall/*.variant.regions.annotated.primary.alleles.report.HighCov.tab > Reports/All."$run".HighCov.primary.alleles.report.tab
		cat result/*/variantcall/*.variant.regions.annotated.primary.alleles.report.LowCov.tab > Reports/All."$run".LowCov.primary.alleles.report.tab
	    	cat result/*/variantcall/*.variant.regions.annotated.minority.alleles.report.HighCov.tab > Reports/All."$run".HighCov.minority.alleles.report.tab
		cat result/*/variantcall/*.variant.regions.annotated.minority.alleles.report.LowCov.tab > Reports/All."$run".LowCov.minority.alleles.report.tab
		cat result/*/variantcall/*Nstats* > Reports/All."$run".Nstats.tab
	        cat result/*/Mapping/*alignment.stats.tab > Reports/All."$run".alignmentStats.tab
	        cat result/*/Mapping/*average.coverage.tab > Reports/All."$run".coverageStats.tab
		cat result/*/variantcall/*classification.tab > Reports/All."$run".NClassificationChart.tab
		pdfunite result/*/variantcall/*.complete.trackplot.HighCov.pdf Reports/All."$run".complete.trackplots.HighCov.pdf
		pdfunite result/*/variantcall/*.complete.trackplot.LowCov.pdf Reports/All."$run".complete.trackplots.LowCov.pdf
	        pdfunite result/*/variantcall/*.nocoverage.trackplot.HighCov.pdf Reports/All."$run".nocov.trackplots.HighCov.pdf
		pdfunite result/*/variantcall/*.nocoverage.trackplot.LowCov.pdf Reports/All."$run".nocov.trackplots.LowCov.pdf
	        pdfunite result/*/variantcall/*nouptrim.coverage.pdf Reports/All."$run".pre-trim.coverage.pdf
	        echo "Outputting the consensus sequences divided in LowCov (Ns > 10%) and HighCov (Ns < 10%) at $dir/Sequences/.."
		mkdir -p Sequences/HighCov
		mkdir -p Sequences/LowCov
		mv result/*/variantcall/*.consensus.HighCov.fasta Sequences/HighCov
		mv result/*/variantcall/*.consensus.LowCov.fasta Sequences/LowCov 
	#	echo "..done, Now attaching the reference sequence and aligning all the consensus to the reference with muscle.."
	#	cat /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/reference.fasta Sequences/HighCov/* > Sequences/HighCov/All.HighCovSeq.fasta
	#	muscle -in Sequences/HighCov/All.HighCovSeq.fasta -diags -quiet -out Sequences/HighCov/All.HighCovSeq.alignment
	#	cat /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/reference.fasta Sequences/LowCov/* > Sequences/LowCov/All.LowCovSeq.fasta
	#        muscle -in Sequences/LowCov/All.LowCovSeq.fasta -diags -quiet -out Sequences/LowCov/All.LowCovSeq.alignment
	#	echo "..done." 
		echo "Outputting the bam files for mapped, unmapped and trimmed reads at $dir/Alignments/-mapped,unmapped,trimmed-.."
		mkdir -p Alignments/mapped/HighCov
		mkdir -p Alignments/mapped/LowCov
		mkdir -p Alignments/unmapped/HighCov
		mkdir -p Alignments/unmapped/LowCov
		mkdir -p Alignments/coverage_trimmed/HighCov
		mkdir -p Alignments/coverage_trimmed/LowCov
	
		mv result/*/Mapping/*.alignment.removed.duplicates.mapped.reads.only.sorted.HighCov* Alignments/mapped/HighCov
		mv result/*/Mapping/*.alignment.removed.duplicates.mapped.reads.only.sorted.LowCov* Alignments/mapped/LowCov
		mv result/*/Mapping/*.alignment.removed.duplicates.unmapped.reads.only.sorted.HighCov* Alignments/unmapped/HighCov
		mv result/*/Mapping/*.alignment.removed.duplicates.unmapped.reads.only.sorted.LowCov* Alignments/unmapped/LowCov
		mv result/*/Mapping/*.alignment.removed.duplicates.mapped.1000trimmed.sorted.HighCov* Alignments/coverage_trimmed/HighCov
		mv result/*/Mapping/*.alignment.removed.duplicates.mapped.1000trimmed.sorted.LowCov* Alignments/coverage_trimmed/LowCov
		
		mkdir Single_output
		mv result/ Single_output
		mv slurm* Single_output
		mv snp* Single_output
	
		cd Reports/
		cp /scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/QCsheet_Illumina_v01_DT.Rmd .
		echo "..done, outputting the QC plots and reports.."
		Rscript -e "rmarkdown::render(input ='QCsheet_Illumina_v01_DT.Rmd')"
		#echo "Now making the directories and files accessible and modifiable by all of us.."
		mv QCsheet_Illumina_v01_DT.html "$run".html
		echo "..done! for doublecheck go to  $dir/slurm_*.out"
        	cd ..
		break
  	fi
	done

	rm *flag
	cd ../../
	echo "Now making the directory $run accessible by all of us (permission 775 will apply).."
	chmod 775 -R $run/
	date
	time2=$(date | awk '{print $4}')
	sec2=`date +%s -d ${time2}`
	diffsec=`expr ${sec2} - ${sec1}`
	echo Start ${time1}
	echo Finish ${time2}
	echo Took ${diffsec} seconds
	echo Total time elapsed `date +%H:%M:%S -ud @${diffsec}`
	echo "In case of troubles please contact me at alfredo.mari@unibas.ch"
	echo "Cheers!"
fi
