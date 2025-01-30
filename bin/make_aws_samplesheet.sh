#!/usr/bin/env bash
aws_dir=$1
run_name=$2
r1_ext=$3
r2_ext=$4
aws s3 ls ${aws_dir} | grep "fastq.gz" | rev | cut -d " " -f1 | rev > ${run_name}_aws_file_list.txt
echo "sample,fastq_1,fastq_2" > ${run_name}_samplesheet.csv
IFS=$'\n\r'
for line in $(cat ${run_name}_aws_file_list.txt | grep "${r1_ext}")
do
	base=`basename $line $r1_ext`
	echo "${base},${aws_dir}${base}${r1_ext},${aws_dir}${base}${r2_ext}" >> ${run_name}_samplesheet.csv 
done