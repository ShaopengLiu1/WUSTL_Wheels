# to run full pipeline of Drop-seq analysis
date

# read parameters
while getopts g:o:i:h opts
do case "$opts" in
g) species="$OPTARG";;    # mm10 for now
o) R1="$OPTARG";;    # PE read 1, or the SE file
i) index="$OPTARG";;    # the meta infor for reads file
h) echo "
Analyzing DROP-seq data by following the McCroll lab instructions

Usage: bash pipe.sh -g <mm10/hg38> -i <index_file(read1)> -o <read_file(read2)>
"
exit;;
esac
done

[ -z "$R1" ] && echo "please specify out the read file" && exit
[ -z "$index" ] && echo "please specfiy out the meta file for read file" && exit
echo "the input file is $index	$R1"

# load executives
drop_tools=/home/shaopengliu/pipe_script/build_a_wheel/drop_seq/Drop-seq_tools-2.0.0
star=/usr/bin/STAR
picard_tools=/home/shaopengliu/pipe_script/build_a_wheel/drop_seq/Drop-seq_tools-2.0.0/jar/lib/picard-2.18.14.jar


# load reference file
[ -z "$species" ] && echo "please specify out the species for read file" && exit

if [[ $species == "mm10" ]];then
  star_ref='/home/Resource/Genome/mm10/STAR_index_mm10_v2.5.4b.gencode.vM15'
  fasta="/home/shaopengliu/resources/mm10/drop_seq_ref/mm10.fa"
elif [[ $species == "hg38" ]];then
  star_ref=star_ref='/home/Resource/Genome/hg38/STAR_index_gencode.v27.annotation'
  fasta="/home/shaopengliu/resources/hg38/drop_seq_ref/hg38.fa"
fi
	
if [[ $R1 == *.fastq || *.fastq.gz || *.fq.gz ]]; then
    name=`echo ${R1%.fastq*}`
    name=`echo ${name%.fq.gz}`
else
    echo "please use fastq file as input"
fi

# step0, preparation and check file exist for completing mode
if [ -d Processed_${name} ]; then
    echo "The processed file is existing"
    exit
else
    mkdir Processed_${name}
    ln -s `realpath $R1`    ./Processed_${name}/${R1}
    ln -s `realpath $index` ./Processed_${name}/${index}
    cd Processed_${name}
fi 

# step1, concate the index fastq and read fastq to bam as the input for pipe
java -jar $picard_tools FastqToSam F1=$index F2=$R1 O=unaligned_start.bam SM=$name

# step2, run alignment code
mkdir output_files
mkdir temp_files
$drop_tools/Drop-seq_alignment.sh -g $star_ref  -r  $fasta -d $drop_tools -o output_files -t temp_files -k unaligned_start.bam 

# step3, extract expression table
export TMPDIR="`pwd`/temp_files" 
cd output_files
$drop_tools/DigitalExpression I=final.bam O=out_gene_exon_tagged.dge.txt.gz SUMMARY=out_gene_exon_tagged.dge.summary.txt   MIN_NUM_TRANSCRIPTS_PER_CELL=500
unset TMPDIR
cd ../..
echo "whole pipe done"
date





