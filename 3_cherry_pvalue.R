## author : Marie MICHEL 
## date : 07/08/2019
## aim : pvalue correction 
## input : pvalue table 
## output : list of COGs differentially expressed 
## packages : cherry


rm(list=ls())
library(tidyverse)
library(cherry)
library(here)

## read data
directory <- "results"
## load p-value table (p-value for DE for each species pairs)
load(here(directory,"tab_pval_length.RData"))

## remove COG with no p-value for any family
#tab_pval_tmp = tab_pval
tab_pval = tab_pval[complete.cases(tab_pval),]

par(mfrow=c(2,2))
for(i in 2:8){
hist(tab_pval[,i])
}
#tab_pval_tmp[is.na(tab_pval_tmp)]<-1
#which(rowSums(tab_pval_tmp[,2:8])==7)
#tab_pval= tab_pval[-which(rowSums(tab_pval_tmp[,2:8])==7),]

#on met tout en lignes ==> trois colonnes : COG - couple - pvalue
gather_pvalue_tmp <- gather(tab_pval, colnames(tab_pval)[-1], key="pair", value="pvalue")

## remove NA 
gather_pvalue <- drop_na(gather_pvalue_tmp)


#correction d'hommel
hom <- hommelFast(gather_pvalue$pvalue,simes=TRUE)

#donne les vrais positifs
cogs_H1 <- sapply(tab_pval$COG_ID, FUN=function(i) pickSimes(hom, select=gather_pvalue$COG_ID == i))

# cogs_H1_list <- list()
# for(ut in tab_pval$COG_ID){
#   print(ut)
#   cogs_H1_list[ut] <- pickSimes(hom, select=gather_pvalue$COG_ID == ut)
# }
# 
# gather_pvalue[gather_pvalue$COG_ID == "3564",]

#tableau de résultat
res <- data.frame(COG_ID=tab_pval$COG_ID,  TP_pred=cogs_H1)

## join results in one table
tab_pval_cherry<- full_join(tab_pval,res,by="COG_ID")
save(tab_pval_cherry,file="tab_pval_cherry.RData")


#permet de voir combien de COGs appartiennent à n couples
nb_H0_reject_effectif=table(res$TP_pred)
save(nb_H0_reject_effectif,file="nb_H0_reject_effectif.RData")

# save lists 
COGs_DE_cherry_1 <- res$COG_ID[which(res$TP_pred>0)]
COGs_DE_cherry_2 <- res$COG_ID[which(res$TP_pred>1)]
COGs_DE_cherry_3 <- res$COG_ID[which(res$TP_pred>2)]
length(COGs_DE_cherry_3)
save(COGs_DE_cherry_1,file="COGs_DE_cherry_1.RData")
save(COGs_DE_cherry_2,file="COGs_DE_cherry_2.RData")
save(COGs_DE_cherry_3,file="COGs_DE_cherry_3.RData")

