#! /bin/bash

# 2020-09-02

####### Input files
expression_matrix=$1
genome_sequence=$2
genome_annotation=$3
rnaseq=$4

####### Software
softwareDir="/home/galaxy/tools/easyMF/"

####### Parameters
thread=$5
##### fastp
min_read_length=$6
quality_value=$7
##### HISAT2
min_intron_length=$8
max_intron_length=$9

####### Create work direction
curDir=$(date +%s%N | cut -c 1-12)
resDir="/home/galaxy/tools/easyMF/01_MATRIX_PREPARATION/$curDir/"

mkdir -p ${resDir}01_genome ${resDir}02_raw_data ${resDir}03_clean_data ${resDir}04_hisat_result ${resDir}05_uniq_map ${resDir}06_expression/gene_exp ${resDir}06_expression/trans_exp


####### Commands
cat ${rnaseq} > ${resDir}Raw_fastq.tar
tar -xvf ${resDir}Raw_fastq.tar -C ${resDir}02_raw_data

${softwareDir}hisat2-2.2.0/hisat2-build -p ${thread} ${genome_sequence} ${resDir}01_genome/genome_hisat_index
hisat_genome_index="${resDir}01_genome/genome_hisat_index"

${softwareDir}samtools-1.9/samtools faidx ${genome_sequence}


dataID=$( ls ${resDir}02_raw_data | sed -r 's/_.+//g' | sort -u )

for i in $dataID
do

	## For single-end reads
	if [ `ls ${resDir}02_raw_data | grep "${i}" | wc -l` = 1 ]; then
		
		#### Check the sequence quality
		phredval=$( cat ${resDir}02_raw_data/${i}.fastq | perl ${softwareDir}/fastq_phred.pl  - | grep -o "Phred.." | sed "s/Phred//" )

		if [ ${phredval} = 33 ]; then
		
			#### Trim low-quality reads
			${softwareDir}fastp/fastp -i ${resDir}02_raw_data/${i}.fastq -o ${resDir}03_clean_data/${i}_trim.fastq -l ${min_read_length} --cut_mean_quality ${quality_value} -j ${resDir}03_clean_data/${i}_fastp.json -h ${resDir}03_clean_data/${i}_fastp.html -w ${thread}

			#### Alignment
			${softwareDir}hisat2-2.2.0/hisat2 --min-intronlen ${min_intron_length} --max-intronlen ${max_intron_length} -k 20 -p ${thread} --phred33 ${hisat_genome_index} -U ${resDir}03_clean_data/${i}_trim.fastq -S ${resDir}04_hisat_result/${i}_hisatResult.sam --summary-file ${resDir}04_hisat_result/${i}_hisatInfo

		else

			#### Trim low-quality reads
			${softwareDir}fastp/fastp --phred64 -i ${resDir}02_raw_data/${i}.fastq -o ${resDir}03_clean_data/${i}_trim.fastq -l ${min_read_length} --cut_mean_quality ${quality_value} -j ${resDir}03_clean_data/${i}_fastp.json -h ${resDir}03_clean_data/${i}_fastp.html -w ${thread}

			#### Alignment
			${softwareDir}hisat2-2.2.0/hisat2 --min-intronlen ${min_intron_length} --max-intronlen ${max_intron_length} -k 20 -p ${thread} --phred64 ${hisat_genome_index} -U ${resDir}03_clean_data/${i}_trim.fastq -S ${resDir}04_hisat_result/${i}_hisatResult.sam --summary-file ${resDir}04_hisat_result/${i}_hisatInfo

		fi

	## For paired-end reads
	else 

		#### Check the sequence quality
		phredval=$( cat ${resDir}02_raw_data/${i}_1.fastq | perl ${softwareDir}/fastq_phred.pl  - | grep -o "Phred.." | sed "s/Phred//" )

		if [ ${phredval} = 33 ]; then
		
			#### Trim low-quality reads
			${softwareDir}fastp/fastp -i ${resDir}02_raw_data/${i}_1.fastq -I ${resDir}02_raw_data/${i}_2.fastq -o ${resDir}03_clean_data/${i}_1_paired.fastq -O ${resDir}03_clean_data/${i}_2_paired.fastq -l ${min_read_length} --cut_mean_quality ${quality_value} -j ${resDir}03_clean_data/${i}_fastp.json -h ${resDir}03_clean_data/${i}_fastp.html -w ${thread}

			#### Alignment
			${softwareDir}hisat2-2.2.0/hisat2 --min-intronlen ${min_intron_length} --max-intronlen ${max_intron_length} -k 20 -p ${thread} --phred33 ${hisat_genome_index} -1 ${resDir}03_clean_data/${i}_1_paired.fastq -2 ${resDir}03_clean_data/${i}_2_paired.fastq -S ${resDir}04_hisat_result/${i}_hisatResult.sam --summary-file ${resDir}04_hisat_result/${i}_hisatInfo

		else

			#### Trim low-quality reads
			${softwareDir}fastp/fastp --phred64 -i ${resDir}02_raw_data/${i}_1.fastq -I ${resDir}02_raw_data/${i}_2.fastq -o ${resDir}03_clean_data/${i}_1_paired.fastq -O ${resDir}03_clean_data/${i}_2_paired.fastq -l ${min_read_length} --cut_mean_quality ${quality_value} -j ${resDir}03_clean_data/${i}_fastp.json -h ${resDir}03_clean_data/${i}_fastp.html -w ${thread}

			#### Alignment
			${softwareDir}hisat2-2.2.0/hisat2 --min-intronlen ${min_intron_length} --max-intronlen ${max_intron_length} -k 20 -p ${thread} --phred64 ${hisat_genome_index} -1 ${resDir}03_clean_data/${i}_1_paired.fastq -2 ${resDir}03_clean_data/${i}_2_paired.fastq -S ${resDir}04_hisat_result/${i}_hisatResult.sam --summary-file ${resDir}04_hisat_result/${i}_hisatInfo

		fi

	fi
	#### Extract Unique-mapped reads
	grep -w "NH:i:1" ${resDir}04_hisat_result/${i}_hisatResult.sam > ${resDir}05_uniq_map/${i}_hisatResult_uniq.sam
	${softwareDir}samtools-1.9/samtools view -@ ${thread} -bt ${genome_sequence}.fai ${resDir}05_uniq_map/${i}_hisatResult_uniq.sam > ${resDir}05_uniq_map/${i}_hisatResult_uniq.bam
	${softwareDir}samtools-1.9/samtools sort -@ ${thread} ${resDir}05_uniq_map/${i}_hisatResult_uniq.bam -o ${resDir}05_uniq_map/${i}_hisatResult_uniq_sort.bam
	${softwareDir}samtools-1.9/samtools view -@ ${thread} -bS ${resDir}04_hisat_result/${i}_hisatResult.sam > ${resDir}04_hisat_result/${i}_hisatResult.bam

	##### Estimate experssion level
	${softwareDir}stringtie-1.3.4d.Linux_x86_64/stringtie -e -p ${thread} -A ${resDir}06_expression/gene_exp/${i}_gene.gtf -G ${genome_annotation} -o ${resDir}06_expression/trans_exp/${i}.gtf ${resDir}05_uniq_map/${i}_hisatResult_uniq_sort.bam

done 

Rscript ${softwareDir}01_MATRIX_PREPARATION/02_integrate_expression_value.R --input1 ${resDir}06_expression/gene_exp --output1 ${expression_matrix}
rm -rf ${resDir}

