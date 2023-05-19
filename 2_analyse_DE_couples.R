## author : Marie MICHEL
## date : 01/08/2019
## aim : differential expression analysis for couples of species
## input : count data, condition table for each couples
## output : list of COGs differentially expressed and plots of counts per species
## packages : DESeq2, upsetR, FactoMineR, dplyr

rm(list=ls())

library(DESeq2)
library(FactoMineR)
library(dplyr)
library(UpSetR)
library(tidyr)


########################################################################
## Option longueur
########################################################################
length <- "length"
# length <- "no_length"

########################################################################
## Ouverture des fichiers 
########################################################################

mypath <- getwd()
dir.create("./results")
setwd("./results")
dir.create("./PCA")
dir.create("./dispersion")
dir.create("./MAplot")
dir.create("./hist_pval")
getwd()
setwd(dir = "../data")
getwd()

setwd(dir = "./conditions")
myfiles_conditions<- list.files(pattern = ".RData", full.names=TRUE)

if (length =="length"){
  setwd(dir = "../length")
  myfiles_length<- list.files(pattern = ".RData", full.names=TRUE)
}

setwd(dir = "../sum_counts")
myfiles_counts<- list.files(pattern = ".RData", full.names=TRUE)

#initialisation pour enregistrer les pvalues
pval_list = list()
tmp_names = NULL



list_for_grep <- c("Hkan","Gpru", "Pgab", "Phco","Pcos","Nmon","Scor")



name_spe="Phco"
i=1
for( name_spe in list_for_grep){  
  
  setwd(mypath)
  setwd("./data")
  
  
  # load counts 
  setwd(dir = "./sum_counts")
  load(myfiles_counts[grep(paste0("(",name_spe,")"), myfiles_counts)])
  ## remove COG with NA
  print(table(complete.cases(sum_counts_COGspartage)))
  countData <- sum_counts_COGspartage[complete.cases(sum_counts_COGspartage),]
  
  # load conditions
  setwd(dir = "../conditions")
  load(myfiles_conditions[grep(paste0("(",name_spe,")"), myfiles_conditions)])
  conditions <- data.frame(conditions)
  
  # load length
  if (length =="length" & name_spe!="Scor"){
    setwd(dir = "../length")
    load(myfiles_length[grep(paste0("(",name_spe,")"), myfiles_length)])
    dim(sum_COGslength)
    ## remove COG with NA
    table(complete.cases(sum_COGslength))
    matrice_length = sum_COGslength[complete.cases(sum_counts_COGspartage),]
    table(complete.cases(matrice_length))
    }
  
  
  
  
  #création du ColData pour l'utilisation de DESeq2
  #ici il faut faire attention à ce que les espèces soient dans le même ordre
  #si les noms sont déjà identiques entre les fichiers et colData c'est bon
  #dans le doute on ordonne par ordre alphanumérique et on renomme 
  #(il faut que ça soit dans le meme ordre)
  #head(colData)
  conditions <- conditions[order(rownames(conditions)),]
  colData <- conditions[,"Hyperaccu"]
  colData <- data.frame(colData)
  colnames(colData) <- "Hyperaccu"
  rownames(colData) <- rownames(conditions)
  
  
  
  countData<- countData[,order(colnames(countData))]
  matrice_length <- matrice_length[,order(colnames(matrice_length))]
  ##sum_counts_COGspartage <- sum_counts_COGspartage[,order(colnames(sum_counts_COGspartage))]
  #colnames(sum_counts_COGspartage) <- rownames(colData)
  
  print("check")
  if (length =="length" & name_spe!="Scor"){print(colnames(matrice_length))}
  print(colnames(countData))
  print(colData)
  
  
  ########################################################################
  ## format analyse diff
  ########################################################################
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design= ~Hyperaccu)
  ## normalisation factors for library size
  dds <-  estimateSizeFactors(dds)
  size_fac <- sizeFactors(dds)
  
  ## normalized counts for library size
  norm_lib_counts <- counts(dds,normalized=TRUE)
  ## check consistancy between
  norm_lib_counts[1,1]
  sum_counts_COGspartage[1,1]/size_fac[1]
  
  
  if (length=="length" & name_spe!="Scor"){
    ## compute normalization with length... 
    mat_size_fac <- matrix(size_fac,ncol=length(size_fac),nrow=length(matrice_length[,1]),byrow=T)
    normFactors <- (mat_size_fac*matrice_length) / exp(rowMeans(log(mat_size_fac*matrice_length)))
    
    
    normalizationFactors(dds) <- as.matrix(normFactors)
    normFactors["1006",]
    sum_counts_COGspartage["1006", ]
    sum_COGslength["1006",]
    
    ## check consistency 
    norm_lib_length_counts <- counts(dds,normalized=TRUE)
    norm_lib_counts[1,1]
    norm_lib_length_counts[1,1]
    
  }
  
  nom_espece <- unique(conditions$ID.espece)
  ###########################################################
  ## ACP
  ###########################################################
  setwd(dir = "../../results")
  rld <- vst(dds)
  
  setwd("./PCA")
  pdf(paste0("PCA_", nom_espece[1],"_",length(nom_espece),"_",length,".pdf"))
  print(plotPCA(rld, intgroup="Hyperaccu"))
  dev.off()
  
  ###########################################################
  ## DE analysis
  ###########################################################
  ## runing analysis
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  res <- results(dds, lfcThreshold=0, name="Hyperaccu_oui_vs_non")
  
  ##diagnostics
  #dispersion
  setwd("../dispersion")
  pdf(paste0("dispersion_", nom_espece[1],"_",length(nom_espece),"_",length,".pdf"))
  plotDispEsts(dds,main=unique(conditions$ID.espece))
  dev.off()
  
  #MA plot
  setwd("../MAplot")
  pdf(paste0("MAplot_", nom_espece[1],"_",length,"_",length(nom_espece),".pdf"))
  plotMA(dds,main=unique(conditions$ID.espece))
  dev.off()
  
  #Histogram of pvalue
  setwd("../hist_pval")
  pdf(paste0("hist_pval_",nom_espece[1],"_",length(nom_espece),"_",length,".pdf"))
  hist(res$pvalue,xlab="p-value",main=paste0("Histogram of p-value ", unique(conditions$ID.espece),"_", length),col="blue")
  dev.off()
  
  #enregistrement des pvalues
  tmp_names = c(tmp_names, paste0(as.character(nom_espece[1]),"_",length(nom_espece)))
  pval = data.frame(COG_ID = as.character(rownames(res)), "value" = as.numeric(res$pvalue), stringsAsFactors = F)
  
  pval_list[[i]] = pval
  i=i+1
}

tab_pval = pval_list[[1]]

#transforme la liste en une table de pvalue
for(i in 2:length(list_for_grep))
{
  tab_pval = full_join(tab_pval, pval_list[[i]], by = "COG_ID")
}
setwd("../")
colnames(tab_pval) = c("COG_ID",tmp_names)
save(tab_pval, file=paste0("tab_pval_", length,".RData"))

setwd(mypath)
