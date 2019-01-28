args=commandArgs()

if(grepl("-h",args[6])) {
  #############
  # help page #
  #############
  cat("\n####################################################\n")
  cat("# Rscript for performing RUV analysis and generating input table for DEG code#")
  cat("\n####################################################\n")
  cat("\nVersion: 1.0 (Last modified on 10-11-2018)\n")
  cat("\nUsage: Rscript RUV_for_deg.R exp_table design_table\n")
  cat("\ndesign_table should have 3 columns: name treatment contrast")
  cat("\ndesign_table is same to what will be used in DEG")
  cat("\n#######################\n")
  cat("# Thanks for using!!! #")
  cat("\n#######################\n")
  cat("\n")
}
  
suppressMessages(library(RUVSeq))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

input_file=args[6]
design_file=args[7]

colors <- brewer.pal(3, "Set1")
fCount <- suppressMessages(read_delim(input_file, "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE))
fCount <- as.data.frame(fCount)
GeneNames <- fCount[,1]
gene_length<- fCount[, c(1,2)]

fCount <- fCount[,-c(1:2)]
fCount <- as.data.frame(fCount)
row.names(fCount) <- GeneNames
design_table <- suppressMessages(read_delim(design_file, "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE))
design_table <- as.data.frame(design_table)
row.names(design_table) <- design_table[,1]
design_table <- design_table[,-1]
design_table <- as.data.frame(design_table)

design=model.matrix(~design_table$contrast)

y=DGEList(counts=fCount,group=design_table$contrast)
#filter low expressed genes
keep=rowSums(cpm(y)>1)>=as.integer(0.5*length(design_table$contrast))
y=y[keep,]
y=calcNormFactors(y, method="upperquartile")
y=estimateGLMCommonDisp(y, design)
y=estimateGLMTagwiseDisp(y, design)
set=newSeqExpressionSet(as.matrix(y$counts),phenoData=design_table)

fit=glmFit(y, design)
res=residuals(fit, type="deviance")
seqUQ=betweenLaneNormalization(set, which="upper")
controls=rownames(set)


rm_RUV <- function(k) {
    seqRUVr=RUVr(seqUQ, controls, k=k, res) ######## k is the number of removed variance.
    ndddd=(normCounts(seqRUVr))
    temp_out<-as.data.frame(ndddd)
    temp_col<-colnames(temp_out)
    temp_out$name=row.names(temp_out)
    temp_out<-merge(temp_out, gene_length, by="name", all.x=T, all.y=F)
    temp_out<-temp_out[ ,c("name","length",temp_col)]
    write.table(temp_out, file = paste0("out_RUV_k", k, "_out.txt"), sep="\t", quote=F, row.names = F)
    suppressMessages(pdf(paste0("out_RUV_k", k, "_out.pdf")))
    suppressMessages(plotPCA(ndddd))
    invisible(dev.off())
    suppressMessages(pdf(paste0("RLE_out_RUV_k", k, "_out.pdf")))
    suppressMessages(plotRLE(ndddd, outline=FALSE, col=colors[as.factor(design_table$contrast)], ylim=c(-1,1)))
    invisible(dev.off())
}

#plot normalized raw PCA
raw_data=(normCounts(seqUQ))
suppressMessages(pdf("out_norm_raw_PCA.pdf"))
suppressMessages(plotPCA(raw_data))
invisible(dev.off())
suppressMessages(pdf("RLE_out_norm_raw_PCA.pdf"))
suppressMessages(plotRLE(raw_data, outline=FALSE, col=colors[as.factor(design_table$contrast)], ylim=c(-1,1)))
invisible(dev.off())

for (rmk in seq(1,dim(design_table)[1]-1,1)){
  rm_RUV(rmk)
}

#DESeq2 differential expression analysis on batch corrected
# use the DEG code for output table

