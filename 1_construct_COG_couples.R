## author : Marie MICHEL
## date : 16/08/2019
## input 1 : for each species, a .csv file containing reads counts for each contigs for each sample
## input 2 : the file "new_groups_all.txt" that contains the association between contigs for each species and COG groups
## input 3 : the file "conditions_all.csv" that contains all metadata about samples from Garcia de la Torre et al.  and Meier et al.
## place all input files in a "data" folder
## output : for each pairs or trio or species, a count data table at the level of COGs ready for differential analysis + length table
## packages : dplyr tidyr reshape2

## This file create a count data table at the level of COGs ready for differential analysis + length table
## for each pairs of species (Garcia de la Torre et al.)
# list_names_pair <- c("Hkan|Hbet","Gpru|Grac", "Pgab|Psem", "Phco|Pluc", "Lhav|Lmic")
# list_names_pair <- c("Hkan|Hbet","Gpru|Grac", "Pgab|Psem", "Phco|Pluc") si on supprime le couple avec un seul échantillon
## for each trio of species (Garcia de la Torre et al.)
# list_names_pair <- c("Pcos|Pgra|Pcle", "Nmon|Mper|Ncfi")
## for the specie Scenecio (Meier et al.)
# Scor

rm(list=ls())

library(dplyr)
library(tidyr)
library(reshape2)


## mypath contains the working directory 
mypath <- getwd()

setwd(dir = "./data")
dir.create("./conditions")
dir.create("./length")
dir.create("./sum_counts")
dir.create("./visu_groups")

# load contigs to COG attributions contained in new_groups_all.txt
groups_cog_all <- read.table("new_groups_all.txt",sep=";")

groups_cog <- groups_cog_all[,c(3,2)]
colnames(groups_cog)<-c("Contig.ID","COG.ID")

########################################################################
########################################################################
########################################################################
## CREATE COG MATRIX FOR PAIRS OF SPECIES (Garcia de la Torre et al. )
########################################################################
########################################################################
########################################################################

# ########################################################################
# ## select pair
# ########################################################################
#
# liste des especes par paires
# list_names_pair <- c("Hkan|Hbet","Gpru|Grac", "Pgab|Psem", "Phco|Pluc", "Lhav|Lmic")
list_names_pair <- c("Hkan|Hbet","Gpru|Grac", "Pgab|Psem", "Phco|Pluc")
#name_pair="Hkan|Hbet"
for( name_pair in list_names_pair){
  
  paire1 <- strsplit(name_pair, split ="\\|")[[1]][1]
  paire2 <- strsplit(name_pair, split="\\|")[[1]][2]
  
  
  ########################################################################
  ## fichiers comptages
  ########################################################################
  
  #ouverture des fichiers de comptage de chaque espèce
  counts1 <- read.table (paste0(paire1,"_read count and RPKM.csv"), header = TRUE, sep=",")
  counts2 <- read.table (paste0(paire2,"_read count and RPKM.csv"), header = TRUE, sep=",")
  
  #ouverture fichier conditions
  conditions <- read.csv2("conditions_all.csv", header = TRUE, sep=";",stringsAsFactors = TRUE,row.names = 1)
  
  #on garde seulement les comptages bruts
  ## remove columns with RPKM data
  counts1 <- counts1[,-(grep("RPK", colnames(counts1)))]
  counts2 <- counts2[,-(grep("RPK", colnames(counts2)))]
  
  colnames(counts1)[2] <- "Size_1"
  colnames(counts2)[2] <- "Size_2"

  ########################################################################
  ## première espèce
  ########################################################################
  
  #permet de mettre les colonnes des échantillons en ligne et la valeur sur la colonne à coté
  res_gather_1 <- gather(counts1, colnames(counts1)[3:ncol(counts1)],key="counts",value="comptage")
  
  # #verifs
  # dim(Hbet_counts)[1]*3
  # dim(res_gather_Hbet)
  # head(res_gather_Hbet)
  
  #attribue COGs aux contigs
  res_left_join_1 <- res_gather_1 %>% left_join(groups_cog)
  
  # #verifs
  # dim(res_gather_Hbet)
  # dim(res_left_join_Hbet)
  # head(res_left_join_Hbet)
  
  # #nombre de contigs supprimés et gardés
  # barplot(table(complete.cases(res_left_join_Hbet$COG.ID)))
  #
  # #comptages gardés
  # sum(res_left_join_Hbet$comptage[complete.cases(res_left_join_Hbet)])
  # #proportion comptage gardés
  # sum(res_left_join_Hbet$comptage[complete.cases(res_left_join_Hbet)])/sum(res_left_join_Hbet$comptage)
  # str(res_left_join_Hbet$COG.ID[2])

  #supprime les contigs sans COGs
  res_left_join_1_complete <- res_left_join_1[complete.cases(res_left_join_1),]
  #sépare les échantillons
  res_spread_left_join_1_complete  <- res_left_join_1_complete %>%  spread(counts, comptage)

  ########################################################################
  ## deuxième espèce
  ########################################################################

  res_gather_2 <- gather(counts2, colnames(counts2)[3:ncol(counts2)],key="counts",value="comptage")

  # #verifs
  # dim(Hbet_counts)[1]*3
  # dim(res_gather_Hbet)
  # head(res_gather_Hbet)

  #attribue COGs aux contigs
  res_left_join_2 <- res_gather_2 %>% left_join(groups_cog)

  # #verifs
  # dim(res_gather_Hbet)
  # dim(res_left_join_Hbet)
  # head(res_left_join_Hbet)

  # #nombre de contigs supprimés et gardés
  # barplot(table(complete.cases(res_left_join_Hbet$COG.ID)))
  #
  # #comptages gardés
  # sum(res_left_join_Hbet$comptage[complete.cases(res_left_join_Hbet)])
  # #proportion comptage gardés
  # sum(res_left_join_Hbet$comptage[complete.cases(res_left_join_Hbet)])/sum(res_left_join_Hbet$comptage)
  # str(res_left_join_Hbet$COG.ID[2])

  #supprime les contigs sans COGs
  res_left_join_2_complete <- res_left_join_2[complete.cases(res_left_join_2),]
  #sépare les échantillons
  res_spread_left_join_2_complete  <- res_left_join_2_complete %>%  spread(counts, comptage)

  ########################################################################
  ## Regroupement des deux espèces selon les COGs VISUALISATION
  ########################################################################

  # big_2especes = bind_rows(res_left_join_1, res_left_join_2)
  # head(big_2especes)
  # mat_final = big_2especes %>%  spread(counts, comptage)

  #mat_cog = mat_final %>%  dcast(COG.ID ~ counts, sum)


  ########################################################################
  ## Regroupement des deux espèces selon les COGs VISUALISATION
  ########################################################################
  # créé une matrice avec les COGs en commun, met des NA sinon
  COGs_partage <- res_spread_left_join_1_complete %>% full_join(res_spread_left_join_2_complete, by="COG.ID")

  #pour que ce soit plus clair on met des NA sur les lignes dupliquées (quand il y a une différence de contigs attribués à un COGs dans les deux espèces = duplication)
  #le code s'adapte au nombre d'échantillons
  COGs_partage[which(duplicated(COGs_partage[, 1:(ncol(counts1)+1)])==TRUE), c(2,4:(ncol(counts1)+1))] <- 0
  COGs_partage[which(duplicated(COGs_partage[, c(3,(ncol(counts1)+2):ncol(COGs_partage))])==TRUE), (ncol(counts1)+3):ncol(COGs_partage) ] <- 0

  #enregistrement du fichier de compatage avec chaque contigs et COGs
  setwd(dir = "./visu_groups")
  write.csv(COGs_partage,file=paste0("COGs_partage_",paire1,"_",paire2,".csv"))
  save(COGs_partage,file=paste0("COGs_partage_",paire1,"_",paire2,".RData"))

  #on change les NA en 0 pour pouvoir faire la somme
  #COGs_partage[is.na(COGs_partage)] <- 0

  #on remet les échantillons en ligne et les comptages en colonne
  #pour comptages
  gather_COGs_partage <- gather(COGs_partage[,c(3,4:(ncol(counts1)+1),(ncol(counts1)+4):ncol(COGs_partage))], colnames(COGs_partage)[c(4:(ncol(counts1)+1),(ncol(counts1)+4):ncol(COGs_partage))],key="counts",value="comptage")

  #on somme les comptages par COGs pour chaque échantillon
  sum_counts_COGspartage <- dcast(gather_COGs_partage, COG.ID ~ counts, sum)
  rownames(sum_counts_COGspartage) = sum_counts_COGspartage[,1]
  sum_counts_COGspartage <-sum_counts_COGspartage[,-1]

  # #enregistrement
  setwd(dir = "../sum_counts")
  write.csv(sum_counts_COGspartage,file=paste0("sum_counts_COGspartage_",paire1,"_",paire2,".csv"))
  save(sum_counts_COGspartage,file=paste0("sum_counts_COGspartage_",paire1,"_",paire2,".RData"))

  #pour la taille des COGs
  gather_COGs_length <- gather(COGs_partage[,c(3,2,ncol(counts1)+3)],colnames(COGs_partage[,c(2,ncol(counts1)+3)]),key="size",value="taille")
  sum_COGslength <- dcast(gather_COGs_length, COG.ID ~ size, sum)
  rownames(sum_COGslength) = sum_COGslength[,1]
  sum_COGslength <-sum_COGslength[,-1]
  sum_COGslength <- cbind(sum_COGslength, replicate(ncol(counts1)-2, sum_COGslength[,1]))
  sum_COGslength <- cbind(sum_COGslength, replicate(ncol(counts2)-2, sum_COGslength[,2]))
  sum_COGslength <- sum_COGslength[,-c(1,2)]
  colnames(sum_COGslength) <- c(colnames(counts1)[3:ncol(counts1)], colnames(counts2)[3:ncol(counts2)])

  #enregistrement
  setwd(dir = "../length")
  write.csv(sum_COGslength,file=paste0("sum_COGslength_",paire1,"_",paire2,".csv"))
  save(sum_COGslength,file=paste0("sum_COGslength_",paire1,"_",paire2,".RData"))

  ########################################################################
  ## séparer tableau conditions
  ########################################################################
  conditions <- conditions[grep((paste0("^(",name_pair,")")), rownames(conditions)),]
  
  setwd(dir = "../conditions")
  write.csv(conditions,file=paste0("conditions_",paire1,"_",paire2,".csv"))
  save(conditions,file=paste0("conditions_",paire1,"_",paire2,".RData"))
  
  setwd(dir = "..")
}


########################################################################
########################################################################
########################################################################
## CREATE COG MATRIX FOR TRIOS OF SPECIES (Garcia de la Torre et al. )
########################################################################
########################################################################
########################################################################

setwd(mypath)
setwd("./data")

list_names_pair <- c("Pcos|Pgra|Pcle", "Nmon|Mper|Ncfi")
name_pair <- "Pcos|Pgra|Pcle"
for( name_pair in list_names_pair){
  paire1 <- strsplit(name_pair, split ="\\|")[[1]][1]
  paire2 <- strsplit(name_pair, split="\\|")[[1]][2]
  paire3 <- strsplit(name_pair, split="\\|")[[1]][3]
  
  ########################################################################
  ## fichiers comptages
  ########################################################################
  
  counts1 <- read.table (paste0(paire1,"_read count and RPKM.csv"), header = TRUE, sep=",")
  counts2 <- read.table (paste0(paire2,"_read count and RPKM.csv"), header = TRUE, sep=",")
  counts3 <- read.table (paste0(paire3,"_read count and RPKM.csv"), header = TRUE, sep=",")
  
  #ouverture fichier conditions
  conditions <- read.csv2("conditions_all.csv", header = TRUE, sep=";",stringsAsFactors = TRUE,row.names = 1)
  
  #si il y a 3 échantillons
  counts1 <- counts1[,-(grep("RPK", colnames(counts1)))]
  counts2 <- counts2[,-(grep("RPK", colnames(counts2)))]
  counts3 <- counts3[,-(grep("RPK", colnames(counts3)))]
  
  colnames(counts1)[2] <- "Size_1"
  colnames(counts2)[2] <- "Size_2"
  colnames(counts3)[2] <- "Size_3"
  
  
  ########################################################################
  ## première espèce
  ########################################################################
  
  res_gather_1 <- gather(counts1, colnames(counts1)[3:ncol(counts1)],key="counts",value="comptage")
  
  #attribue COGs aux contigs
  res_left_join_1 <- res_gather_1 %>% left_join(groups_cog)
  
  #supprime les contigs sans COGs
  res_left_join_1_complete <- res_left_join_1[complete.cases(res_left_join_1),]
  #sépare les échantillons
  res_spread_left_join_1_complete  <- res_left_join_1_complete %>%  spread(counts, comptage)
  
  ########################################################################
  ## deuxième espèce
  ########################################################################
  
  res_gather_2 <- gather(counts2, colnames(counts2)[3:ncol(counts2)],key="counts",value="comptage")
  
  #attribue COGs aux contigs
  res_left_join_2 <- res_gather_2 %>% left_join(groups_cog)
  
  #supprime les contigs sans COGs
  res_left_join_2_complete <- res_left_join_2[complete.cases(res_left_join_2),]
  #sépare les échantillons
  res_spread_left_join_2_complete  <- res_left_join_2_complete %>%  spread(counts, comptage)
  
  ########################################################################
  ## troisème espèce
  ########################################################################
  
  res_gather_3 <- gather(counts3, colnames(counts3)[3:ncol(counts3)],key="counts",value="comptage")
  
  #attribue COGs aux contigs
  res_left_join_3 <- res_gather_3 %>% left_join(groups_cog)
  
  #supprime les contigs sans COGs
  res_left_join_3_complete <- res_left_join_3[complete.cases(res_left_join_3),]
  #sépare les échantillons
  res_spread_left_join_3_complete  <- res_left_join_3_complete %>%  spread(counts, comptage)
  
  ########################################################################
  ## Regroupement selon les COGs
  ########################################################################
  COGs_partage <- res_spread_left_join_1_complete %>% full_join(res_spread_left_join_2_complete, by="COG.ID") 
  COGs_partage <- COGs_partage %>% full_join(res_spread_left_join_3_complete, by="COG.ID")
  
  #on met des 0 là où il a des lignes dupliquées
  COGs_partage[which(duplicated(COGs_partage[, 1:(ncol(counts1)+1)])==TRUE), c(2,4:(ncol(counts1)+1))] <- 0
  COGs_partage[which(duplicated(COGs_partage[, c(3,(ncol(counts1)+2):(ncol(counts1)+1+ncol(counts2)))])==TRUE), (ncol(counts1)+3):(ncol(counts1)+1+ncol(counts2)) ] <- 0
  COGs_partage[which(duplicated(COGs_partage[, c(3,(ncol(counts1)+2+ncol(counts2)):ncol(COGs_partage))])==TRUE), (ncol(counts1)+3+ncol(counts2)):ncol(COGs_partage) ] <- 0
  
  
  setwd(dir = "./visu_groups")
  write.csv(COGs_partage,file=paste0("COGs_partage_",paire1,"_",paire2,"_",paire3,".csv" ))
  save(COGs_partage,file=paste0("COGs_partage_",paire1,"_",paire2,"_",paire3,".RData"))
  
  #somme des contigs attribués à chaque COG
  gather_COGs_partage <- gather(COGs_partage[,c(3,4:(ncol(counts1)+1),(ncol(counts1)+4):(ncol(counts1)+1+ncol(counts2)),(ncol(counts1)+4+ncol(counts2)):ncol(COGs_partage))], colnames(COGs_partage)[c(4:(ncol(counts1)+1),(ncol(counts1)+4):(ncol(counts1)+1+ncol(counts2)),(ncol(counts1)+4+ncol(counts2)):ncol(COGs_partage))],key="counts",value="comptage")
  sum_counts_COGspartage <- dcast(gather_COGs_partage, COG.ID ~ counts, sum)
  rownames(sum_counts_COGspartage) = sum_counts_COGspartage[,1]
  sum_counts_COGspartage <-sum_counts_COGspartage[,-1]
  
  setwd(dir = "../sum_counts")
  write.csv(sum_counts_COGspartage,file=paste0("sum_counts_COGspartage_",paire1,"_",paire2,"_",paire3,".csv" ))
  save(sum_counts_COGspartage,file=paste0("sum_counts_COGspartage_",paire1,"_",paire2,"_",paire3,".RData"))
  
  #pour la taille des COGs
  gather_COGs_length <- gather(COGs_partage[,c(3,2,ncol(counts1)+3,ncol(counts1)+3+ncol(counts2))],colnames(COGs_partage[,c(2,ncol(counts1)+3,(ncol(counts1)+3+ncol(counts2)))]),key="size",value="taille")
  sum_COGslength <- dcast(gather_COGs_length, COG.ID ~ size, sum)
  rownames(sum_COGslength) = sum_COGslength[,1]
  sum_COGslength <-sum_COGslength[,-1]
  sum_COGslength <- cbind(sum_COGslength, replicate(ncol(counts1)-2, sum_COGslength[,1]))
  sum_COGslength <- cbind(sum_COGslength, replicate(ncol(counts2)-2, sum_COGslength[,2]))
  sum_COGslength <- cbind(sum_COGslength, replicate(ncol(counts3)-2, sum_COGslength[,3]))
  sum_COGslength <- sum_COGslength[,-c(1,2,3)]
  colnames(sum_COGslength) <- c(colnames(counts1)[3:ncol(counts1)], colnames(counts2)[3:ncol(counts2)], colnames(counts3)[3:ncol(counts3)])
  
  #enregistrement
  setwd(dir = "../length")
  write.csv(sum_COGslength,file=paste0("sum_COGslength_",paire1,"_",paire2,"_",paire3,".csv"))
  save(sum_COGslength,file=paste0("sum_COGslength_",paire1,"_",paire2,"_",paire3,".RData"))
  
  ########################################################################
  ## séparer tableau condition
  ########################################################################
  
  #sélectionne les lignes correspondant aux paires
  conditions <- conditions[grep((paste0("^(",name_pair,")")), rownames(conditions)),]
  
  setwd(dir = "../conditions")
  write.csv(conditions,file=paste0("conditions_",paire1,"_",paire2,"_",paire3,".csv"))
  save(conditions,file=paste0("conditions_",paire1,"_",paire2,"_",paire3,".RData"))
  
  setwd("..")
}

########################################################################
########################################################################
########################################################################
## CREATE COG MATRIX FOR SCOR (Meier et al. )
########################################################################
########################################################################
########################################################################
setwd(mypath)
setwd(dir = "./data")

########################################################################
## fichiers comptages
########################################################################
countsB2 <- read.table ("Scor_meier_read count.csv", header = TRUE, sep="t")
conditions <- read.csv2("conditions_all.csv", header = TRUE, sep=";",stringsAsFactors = TRUE,row.names = 1)

###enlever le "mot" mapping###
countsB2$contig_id_B2 = as.character(countsB2$contig_id_B2)
for (i in 1:nrow(countsB2)) {
  countsB2[i,1] = gsub("_mapping", "", as.character(countsB2[i,1]))
}
countsB2$contig_id_B2 = as.factor(countsB2$contig_id_B2)
colnames(countsB2)[1]<-"Contig.ID"


########################################################################
## deuxième espèce : B2
########################################################################

res_gather_B2 <- gather(countsB2, colnames(countsB2)[3:ncol(countsB2)],key="counts",value="comptage") #met la matrice sous forme de liste

#attribue COGs aux contigs via la fonction join
res_left_join_B2 <- res_gather_B2 %>% left_join(groups_cog)

#supprime les contigs sans COGs
res_left_join_B2_complete <- res_left_join_B2[complete.cases(res_left_join_B2),]

#sépare les échantillons
res_spread_left_join_B2_complete  <- res_left_join_B2_complete %>%  spread(counts, comptage)


gather_COGs_partage <- gather(res_spread_left_join_B2_complete[,3:ncol(res_spread_left_join_B2_complete)], colnames(res_spread_left_join_B2_complete)[4:ncol(res_spread_left_join_B2_complete)],key="counts",value="comptage")
sum_counts_COGspartage <- dcast(gather_COGs_partage, COG.ID ~ counts, sum)  

#Je suis obligé de faire ça car dcast  trie par ordre alphanum donc on une matrice A-B-C-D et on veut A-C-B-D (A et C étant NiH et B et D non)
#sum_counts_COGspartage<-cbind(sum_counts_COGspartage[,1:3],sum_counts_COGspartage[,7:9],sum_counts_COGspartage[,4:6],sum_counts_COGspartage[,10:12])

##On passe COG_id en nom de ligne
rownames(sum_counts_COGspartage) = sum_counts_COGspartage[,1]
sum_counts_COGspartage <-sum_counts_COGspartage[,-1]

########################################################################
## Sauvegarde
########################################################################
setwd("./sum_counts")
write.csv(sum_counts_COGspartage,file="sum_counts_COGpartage_Scor_B2.csv")
save(sum_counts_COGspartage,file="sum_counts_COGpartage_Scor_B2.RData")


########################################################################
## Taille cumulée des COGs
########################################################################

#Fait la somme cumulée des longueurs en fonction du nombre de contigs pas COG
gather_COGs_length <- gather(res_spread_left_join_B2_complete[,c(2,3)], colnames(res_spread_left_join_B2_complete)[2],key="size",value="taille")
sum_COGslength_temp <- dcast(gather_COGs_length, COG.ID ~ size, sum)

#On duplique en fonction du nombre d'échantillon 
sum_COGslength <- sum_COGslength_temp[,c(1,2,rep(2,dim(sum_counts_COGspartage)[2]-1))]
sum_COGslength

#on met les COG en nom de lignes
rownames(sum_COGslength) = sum_COGslength[,1]
sum_COGslength <-sum_COGslength[,-1]


#on redonne les bons noms aux colonnes
colnames(sum_COGslength) <- colnames(countsB2)[3:ncol(countsB2)]
head(sum_COGslength)

##Sauvegarde matrice taille
setwd("../length")
write.csv(sum_COGslength,file="sum_COGslength_Scor_B2.csv")
save(sum_COGslength,file="sum_COGslength_Scor_B2.RData")


setwd("../conditions")
conditions <- conditions[grep("^Scor", rownames(conditions)),]
write.csv(conditions,file="conditions_Scor_B2.csv")
save(conditions,file="conditions_Scor_B2.RData")


## end scripts.
mypath
setwd(mypath)
