date

types=$1
species=mm10

[ -z "$types" ] && echo "no type input" && exit

echo "running $species $types data"

if [[ $types == PE ]]; then
    for file in `ls *_1/fastq.gz`
	do 
	    name=`echo ${file%_1.fastq.gz}`
	    echo "processing $name"
	    singularity run -B /SLOWSTOR:/SLOWSTOR   -B ./:/process /home/shaopengliu/singularity/rna-seq/RNA-seq_mm10_target_181207.simg -r PE -g mm10 -o ${name}_1.fastq.gz -p ${name}_2.fastq.gz
	    unset file name
	done
elif [[ $types == SE ]]; then
    for file in `ls *.fastq.gz`
	do
	    echo "processing $file"
	    singularity run -B ./:/process -B /SLOWSTOR:/SLOWSTOR   /home/shaopengliu/singularity/rna-seq/RNA-seq_mm10_target_181207.simg -r SE -g mm10 -o $file
	    unset file
	done
else 
    echo "the input type is $types, didn't find match"
    exit
fi

echo "whole pipe done"
date
