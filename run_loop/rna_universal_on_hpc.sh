date

types=$1
species=$2

[ -z "$types" ] && echo "no type input" && exit
[ -z "$species" ] && echo "no speceis input" && exit

echo "running $species $types data"

if [[ $types == PE ]]; then
    for file in `ls *_1/fastq.gz`
        do
            name=`echo ${file%_1.fastq.gz}`
            echo "processing $name"
	    bash /home/shaopengliu/pipe_script/github/RNA-seq_QC_analysis/pipe_code/rna_v4.sh -r PE -g $species -o ${name}_1.fastq.gz -p ${name}_2.fastq.gz
	    unset file name
	done
elif [[ $types == SE ]]; then
    for file in `ls *.fastq.gz`
        do
            echo "processing $file"
	    bash /home/shaopengliu/pipe_script/github/RNA-seq_QC_analysis/pipe_code/rna_v4.sh -r SE -g $species -o $file 
	    unset file
	done
else
    echo "the input type is $types, didn't find match"
    exit
fi


echo "whole pipe done"
date


