# help infor 
args=commandArgs()

if(grepl("-h",args[6])) {
  #############
  # help page #
  #############
  cat("\n###########################################################\n")
  cat("# Rscript for performing GO analysis using GOstats and KEGG #")
  cat("\n###########################################################\n")
  cat("\nVersion: 1.0 (Last modified on 10-15-2018)\n")
  cat("\nUsage: Rscript GO_analysis.R <Gene_list_file> <Hs/Mm/Dr>  <GO_Pvalue_Cutoff>\n")
  cat("\nNote: the parameters are positional, which means you have to follow the above order of parameters\n")
  stop("Please use proper parameters to run the code :)")
}


#load packages and parameter
suppressMessages(library(GOstats))
suppressWarnings(suppressMessages(library(biomaRt)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
suppressWarnings(suppressMessages(library(org.Mm.eg.db)))
suppressWarnings(suppressMessages(library(org.Dr.eg.db)))
suppressWarnings(suppressMessages(library(KEGG.db)))
suppressWarnings(suppressMessages(library(gage)))
suppressMessages(library(DBI))


input_file=args[6]
species=args[7]
cutoff=args[8]
name=unlist(strsplit(input_file, split="[.]"))[1]

#prepare data
xegdb=paste0("org.",args[7],".eg.db")
xegGO=paste0("org.",args[7],".egGO")
if(args[7]=="Hs")
{
  GO3Name="hsa"
}else if (args[7]=="Mm")
{
  GO3Name="mmu"
}else if(args[7]=="Dr")
{
  GO3Name="dre"
}

read.table(input_file, header=F, sep="\t", stringsAsFactors=F) -> gene_list



#run GO
suppressMessages(E_IDs<-select(get(xegdb),gene_list$V1,c("ENTREZID","GENENAME"),"ALIAS"))
Entrez=unique(E_IDs[!is.na(E_IDs$ENTREZID),2])
universe=mappedkeys(get(xegGO))

Entrez=Entrez[Entrez%in%universe]
paramsBP=new('GOHyperGParams',geneIds=Entrez,universeGeneIds=universe,ontology='BP',pvalueCutoff=as.numeric(args[8]),conditional=F,testDirection='over',annotation=xegdb)
hgOverBP=hyperGTest(paramsBP)
resultBP <- as.data.frame(summary(hgOverBP))
resultBP$DE_genes <- NA

for( i in 1:nrow(resultBP) )
{
  offspring <- get(resultBP[i,1], GOBPOFFSPRING)
  offspring <- offspring[!is.na(offspring)]
  egids <- unique(unlist(mget(c(resultBP[i,1], offspring), revmap(get(xegGO)), ifnotfound=NA), use.names=FALSE))
  egids <- egids[!is.na(egids)]
  Entrez_dat <- Entrez[Entrez %in% egids]
  Entrez_dat <- Entrez_dat[!is.na(Entrez_dat)]
  symbols <- suppressMessages(select(get(xegdb), keys = Entrez_dat, columns = c("SYMBOL"), keytype = "ENTREZID"))
  resultBP[i,8] <- paste(symbols$SYMBOL, collapse = " ")
}

write.csv(resultBP,file=paste0("GO_BP_in_",name,".csv"))
htmlReport(hgOverBP,file=paste0("GO_BP_in_",name,".html"),digits=5)

paramsCC=new('GOHyperGParams',geneIds=Entrez,universeGeneIds=universe,ontology='CC',pvalueCutoff=as.numeric(args[8]),conditional=F,testDirection='over',annotation=xegdb)
hgOverCC=hyperGTest(paramsCC)
resultCC <- as.data.frame(summary(hgOverCC))
resultCC$DE_genes <- NA
for( i in 1:nrow(resultCC) )
{
  offspring <- get(resultCC[i,1], GOCCOFFSPRING)
  offspring <- offspring[!is.na(offspring)]
  egids <- unique(unlist(mget(c(resultCC[i,1], offspring), revmap(get(xegGO)), ifnotfound=NA), use.names=FALSE))
  egids <- egids[!is.na(egids)]
  Entrez_dat <- Entrez[Entrez %in% egids]
  Entrez_dat <- Entrez_dat[!is.na(Entrez_dat)]
  symbols <- suppressMessages(select(get(xegdb), keys = Entrez_dat, columns = c("SYMBOL"), keytype = "ENTREZID"))
  resultCC[i,8] <- paste(symbols$SYMBOL, collapse = " ")
}
write.csv(resultCC,file=paste0("GO_CC_in_",name,".csv"))
htmlReport(hgOverCC,file=paste0("GO_CC_in_",name,".html"),digits=5)

paramsMF=new('GOHyperGParams',geneIds=Entrez,universeGeneIds=universe,ontology='MF',pvalueCutoff=as.numeric(args[8]),conditional=F,testDirection='over',annotation=xegdb)
hgOverMF=hyperGTest(paramsMF)
resultMF <- as.data.frame(summary(hgOverMF))
resultMF$DE_genes <- NA
for( i in 1:nrow(resultMF) )
{
  offspring <- get(resultMF[i,1],GOMFOFFSPRING)
  offspring <- offspring[!is.na(offspring)]
  egids <- unique(unlist(mget(c(resultMF[i,1], offspring), revmap(get(xegGO)), ifnotfound=NA), use.names=FALSE))
  egids <- egids[!is.na(egids)]
  Entrez_dat <- Entrez[Entrez %in% egids]
  Entrez_dat <- Entrez_dat[!is.na(Entrez_dat)]
  symbols <- suppressMessages(select(get(xegdb), keys = Entrez_dat, columns = c("SYMBOL"), keytype = "ENTREZID"))
  resultMF[i,8] <- paste(symbols$SYMBOL, collapse = " ")
}
write.csv(resultMF,file=paste0("GO_MF_in_",name,".csv"))
htmlReport(hgOverMF,file=paste0("GO_MF_in_",name,".html"),digits=5)


paramsKEGG=new('KEGGHyperGParams',geneIds=Entrez,universeGeneIds=universe,pvalueCutoff=as.numeric(args[8]),testDirection='over',annotation=xegdb)
hgOverKEGG=hyperGTest(paramsKEGG)
resultKEGG <- as.data.frame(summary(hgOverKEGG))
resultKEGG$DE_genes <- NA
for( i in 1:nrow(resultKEGG) )
{
  offspring <- get(paste(GO3Name,resultKEGG[i,1], sep=""), KEGGPATHID2EXTID )
  offspring <- offspring[!is.na(offspring)]
  Entrez_dat <- Entrez[Entrez %in% offspring]
  Entrez_dat <- Entrez_dat[!is.na(Entrez_dat)]
  symbols <- suppressMessages(select(get(xegdb), keys = Entrez_dat, columns = c("SYMBOL"), keytype = "ENTREZID"))
  resultKEGG[i,8] <- paste(symbols$SYMBOL, collapse = " ")
}
write.csv(resultKEGG,file=paste0("KEGG_in_",name,".csv"))



