#!/usr/bin/bash

#Genome-wide variant calling specially for re-sequencing data
#Dependencies: GATK, BWA, Picard, Samtools, Bedtools, SnpEff
#You need to set the PATH of particular tools when it is written in the code as PATH_to_Picard, PATH_to_GATK

#Step1: Mapping fastq files using bwa

genome_dir=/shared_resources/genomes/soybean

for fn in 'ls /shared_resources/shamim_resource/RNASeq_files/*.fastq'; do
    input_dir="/shared_resources/shamim_resource/RNASeq_files"
    base=$(basename $fn "_R1.fastq")
    bwa mem -M -t 16 ${genome_dir}/GMax.fa ${input_dir}/${base}_R1.fastq ${input_dir}/${base}_R2.fastq > Samfiles/${base}.sam
    echo "bwa mem -M -t 16 ${genome_index_dir}/GMax.fa ${input_dir}/${base}_R1.fastq ${input_dir}/${base}_R2.fastq > Samfiles/${base}.sam"
done


#step2:Sort SAM file by coordinate and convert to BAM using Picard tools

java -jar PATH_to_Picard/picard.jar SortSam INPUT=bwa_aligned_reads.sam OUTPUT=sorted_bwa_aligned_reads.bam SORT_ORDER=coordinate


#Step3: Collect Alignment & Insert Size Metrics 
# Dependencies: Picard and Samtools

java -jar PATH_to_Picard/picard.jar CollectAlignmentSummaryMetrics R=${genome_dir}/GMax.fa I=sorted_bwa_aligned_reads.bam O=output_alignment_metrics.txt
java -jar PATH_to_Picard/picard.jar CollectInsertSizeMetrics INPUT=sorted_bwa_aligned_reads.bam OUTPUT=output_insert_metrics.txt HISTOGRAM_FILE=output_insert_size_histogram.pdf
samtools depth -a sorted_bwa_aligned_reads.bam > depth_out.txt


#Step4:Marking Duplicates Using Picard tools
java -jar PATH_to_Picard/picard.jar MarkDuplicates INPUT=sorted_bwa_aligned_reads.bam OUTPUT=out_dedup_reads.bam METRICS_FILE=output_metrics.txt


#Step5: Build BAM Index using Picard tools

java -jar PATH_to_Picard/picard.jar BuildBamIndex INPUT=out_dedup_reads.bam


#Step6: Creating Realignment Targets using GATK
java -jar Path_to_GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${genome_dir}/GMax.fa -I out_dedup_reads.bam -o out_realignment_targets.list

#Step7: 	Realign Indels using GATK
java -jar PATH_to_GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${genome_dir}/GMax.fa -I out_dedup_reads.bam -targetIntervals out_realignment_targets.list -o realigned_reads.bam


#Step8: Variant Calling using GATK
java -jar PATH_to_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${genome_dir}/GMax.fa -I realigned_reads.bam -o output_raw_variants.vcf


#Step9: Extract SNPs & Indels using GATK
java -jar PATH_to_GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${genome_dir}/GMax.fa -V output_raw_variants.vcf -selectType SNP -o unfiltered_snps.vcf
java -jar PATH_to_GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${genome_dir}/GMax.fa -V output_raw_variants.vcf -selectType INDEL -o unfiltered_indels.vcf


#Step10: Filter SNPs and INDELs using GATK
java -jar PATH_to_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${genome_dir}/GMax.fa -V unfiltered_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps.vcf
java -jar PATH_to_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${genome_dir}/GMax.fa -V unfiltered_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels.vcf


# step 11:  Annotate SNPs and Predict Effects using SnpEffs
java -jar snpEff.jar -v snpeff_db filtered_snps.vcf > filtered_snps_final.ann.vcf


###There are other additional steps such as Base Quality Score Recalibration for improving the results
















