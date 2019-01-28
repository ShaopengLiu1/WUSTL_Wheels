# get gene name and gene id by input transcript id

# input: transcript id list (only 1 col), species
# output: a new file with trans id, gene name, gene id

transid=$1
species=$2

source /home/shaopengliu/pipe_script/github/RNA-seq_QC_analysis/pipe_code/rna_seq_qc_source.sh $species


awk '$3=="transcript"' $annotation_file > temp_trans.gtf

cat <<EOF > temp_find_shared.R
args=commandArgs()
id_file=args[6]
gtf_file=args[7]
read.table(id_file, header=F, sep="\t", stringsAsFactors=F) -> trans_id
read.table(gtf_file, header=F, sep="\t", stringsAsFactors=F) -> gtf_anno
gtf_anno<-gtf_anno[c("V1","V4","V5","V9")]

gtf_anno\$gene_name= sapply(strsplit( sapply(strsplit(as.character(gtf_anno\$V9), ';'), '[', 4),' '), '[',3)
gtf_anno\$trans_id= sapply(strsplit( sapply(strsplit(as.character(gtf_anno\$V9),';'), '[', 2), ' '), '[',3)
gtf_anno\$gene_id= sapply(strsplit( sapply(strsplit(as.character(gtf_anno\$V9),';'), '[', 1), ' '), '[',2)

gtf_anno<-gtf_anno[c("V1","V4","V5","gene_id","gene_name","trans_id")]
colnames(gtf_anno)[1:3]<-c("chrom","start","end")
colnames(trans_id)<-"trans_id"
merge(gtf_anno, trans_id, by="trans_id") -> out
write.table(out, file="merged_transid_2_gtf.txt", col.names=F, sep="\t", quote=F, row.name=F)
EOF

Rscript temp_find_shared.R $transid temp_trans.gtf \
	&& rm temp_trans.gtf temp_find_shared.R

awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6}' merged_transid_2_gtf.txt | sort -k1,1V -k2,2n > temp_out \
	&& mv temp_out merged_transid_2_gtf.txt

echo "find gene from transid pipe done"
date

