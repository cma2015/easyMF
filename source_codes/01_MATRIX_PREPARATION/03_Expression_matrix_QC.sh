#! /bin/bash

# 2020-09-02
script_name="03_Expression_matrix_QC.sh"
script_ver="1.0.0"

#Help function
usage() {
	echo "-h  Help documentation for $script_name"
	echo "-R  Raw expression matrix"
	echo "-E  High-quality expression matrix"
	echo "-g  Wether remove lowly expressed genes"
	echo "-e  Expression value of expressed genes"
	echo "-N  Minimum sample number for expression supported"
	echo "-s  Wether remove outliers samples"
	echo "-r  Rethreshold"
	echo "-l  Lowthreshold"
	echo "-b  Wether remove batch effects"
	echo "-S  Sample information"

	exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

while getopts :R:E:g:e:N:s:r:l:b:S:hv opt
	do
		case $opt in
			R) raw_expression_matrix=$OPTARG;;
			E) expression_matrix_QC=$OPTARG;;
			g) operation_gene=$OPTARG;;
			e) min_exp=$OPTARG;;
			N) min_exp_num=$OPTARG;;
			s) operation_sample=$OPTARG;;
			r) rethreshold=$OPTARG;;
			l) lowthreshold=$OPTARG;;
			b) operation_batch=$OPTARG;;
			S) sample_info=$OPTARG;;
			h) usage;;
			v) version;;
		esac
	done

softwareDir="/home/galaxy/tools/easyMF/01_MATRIX_PREPARATION/"

curDir=$(date +%s%N | cut -c 1-12)
resDir="/home/galaxy/tools/easyMF/01_MATRIX_PREPARATION/$curDir/"

mkdir -p ${resDir}
cat ${raw_expression_matrix} > ${resDir}Expression_matrix

############# Genes
if [ ! $operation_gene ]; then
	echo "Removing lowly expressed genes is NULL"
else
	if [ $operation_gene == "genes" ]; then
		Rscript ${softwareDir}03_Remove_lowly_expressed_genes.R \
			--input1 ${resDir}Expression_matrix --input2 ${min_exp} --input3 ${min_exp_num}  \
			--output1 ${resDir}Expression_matrix
	fi
fi

############# Samples
if [ ! $operation_sample ]; then
	echo "Removing outliers samples is NULL"
else
	if [ $operation_sample == "samples" ]; then
		Rscript ${softwareDir}03_Remove_outliers_samples.R \
			--input1 ${resDir}Expression_matrix --input2 ${rethreshold} --input3 ${lowthreshold}  \
			--output1 ${resDir}Expression_matrix
	fi
fi

############# Batch effects
if [ ! $operation_batch ]; then
	echo "Removing batch effects is NULL"
else
	if [ $operation_batch == "batchEffects" ]; then
		Rscript ${softwareDir}03_Remove_batch_effects.R \
			--input1 ${resDir}Expression_matrix --input2 ${sample_info} \
			--output1 ${resDir}Expression_matrix
	fi
fi

Rscript ${softwareDir}03_Log_transform.R \
			--input1 ${resDir}Expression_matrix --output1 ${resDir}Expression_matrix

cat ${resDir}Expression_matrix > ${expression_matrix_QC}
rm -rf ${resDir}

