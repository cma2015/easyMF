#!/bin/bash

# 2020-10-27


software_type=$1
data_access=$3
HQ_fq_retrieve=$2

softwareDir="/home/galaxy/tools/easyMF/"

curDir=$(date +%s%N | cut -c 1-12)
resDir="/home/galaxy/tools/easyMF/01_MATRIX_PREPARATION/$curDir/"


mkdir -p ${resDir}Raw_SRA ${resDir}Raw_fastq


if [ $software_type == "sratoolkit" ]; then
	cat $data_access | while read line
	do
		${softwareDir}sratoolkit.2.9.6-1-ubuntu64/bin/prefetch -O ${resDir}Raw_SRA $line
		${softwareDir}sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --split-3 ${resDir}Raw_SRA/$line.sra -O ${resDir}Raw_fastq/
	done
	
elif [ $software_type == "wget" ]; then
	sed "s#^#wget -P ${resDir}Raw_SRA #g" $data_access > ${resDir}Download.sh
	sh ${resDir}Download.sh

	rawDataID=$(ls ${resDir}Raw_SRA )
	cat $rawDataID | while read line
	do
		newDataID=$( echo $line | sed -r 's/\..+//g' )
		echo "mv ${resDir}Raw_SRA/$line ${resDir}Raw_SRA/$newDataID" >> ${resDir}change_name.sh 

	done
	
	sh ${resDir}change_name.sh 
	${softwareDir}sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --split-3 ${resDir}Raw_SRA/* -O ${resDir}Raw_fastq/
	
fi

tar cvf ${resDir}Raw_fastq.tar.gz -C ${resDir}Raw_fastq .
mv ${resDir}Raw_fastq.tar.gz $2

rm -rf ${resDir}


