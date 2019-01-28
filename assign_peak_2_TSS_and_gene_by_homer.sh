# assign peak file to nearest TSS by homer
# then match transcript ID to gene name and ID

# usage: bash <pipe> <peak_file> <species>
# output: assigned peaks with correlated gene name

input_file=$1
species=$2
name=`echo ${input_file%%.*}`
echo "processing $name"

source /home/shaopengliu/pipe_script/github/RNA-seq_QC_analysis/pipe_code/rna_seq_qc_source.sh $species
echo "the species in $species, and annotation file is $annotation_file"

annotatePeaks.pl $input_file $species -gtf $annotation_file > temp  2> homer_annotate_record.txt 

sed '1d' temp | cut -f1-4,10-11 | sort -k1,1V -k2,2n   > assigned_TSS_from_peak.bed

cut -f 6 assigned_TSS_from_peak.bed > transid_list.txt
bash /home/shaopengliu/pipe_script/build_a_wheel/get_gene_name_id_by_transcript_id.sh transid_list.txt $species \
	&& rm transid_list.txt

# merge the gene id, name back to the Homer output
# just sort by trans id and the paste to merge
sort -k6 assigned_TSS_from_peak.bed  > temp1.txt
cut -f 4-6 merged_transid_2_gtf.txt | sort -k1 > temp2.txt
paste temp1.txt <(cut -f 2-3 temp2.txt) > temp_out  \
	&& mv temp_out assigned_TSS_from_peak.bed  \
	&& rm temp1.txt temp2.txt merged_transid_2_gtf.txt


awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6"\t"$7"\t"$8}'  assigned_TSS_from_peak.bed  | sort -k1,1V -k2,2n > out_temp \
	&& mv out_temp ${name}_assign2TSS.bed \
	&& rm temp 


echo "assign pipe done"
echo "the output file is assigned_TSS_from_peak.bed"
date





