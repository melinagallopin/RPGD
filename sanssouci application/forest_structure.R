load("results/hierarchy.RData")
load("results/tab_pval_length.RData")
## p-values within a unique vector.
vect_couple<-colnames(tab_pval)
number_cog <- tab_pval[,1]
data <-c()
for (i in 2:length(vect_couple)){
  ID <- paste(number_cog,rep(c(vect_couple[i]),each = length(number_cog)))
  p_values <-tab_pval[i,!is.na(tab_pval[i,])][2:length(tab_pval[i,!is.na(tab_pval[i,])])]
  data <- rbind(data,cbind(ID,p_values))
}

basic_tree = c(c())


## Création de tab_pval sans NA (correspond à une paire d'espèce ne possèdant pas le cog associé).----

nona_tab_pval <- tab_pval[!rowSums(is.na(tab_pval))>0,]

data <-c()
vect_couple<-colnames(nona_tab_pval)
number_cog <- nona_tab_pval[,1]

for (i in 1:length(nona_tab_pval[,1])){
  for (j in 2:length(vect_couple)){
  ID <- paste(number_cog[i],c(vect_couple[j]))
  p_value <- nona_tab_pval[i,j]
  data <- rbind(data,c(ID,p_value))
  }}
p_values <-c()
vect_couple<-colnames(nona_tab_pval)
number_cog <- nona_tab_pval[,1]

for (i in 1:length(nona_tab_pval[,1])){
  for (j in 2:length(vect_couple)){
    p_value <- nona_tab_pval[i,j]
    p_values <- rbind(p_values,p_value)
  }}

## Structure d'arbre----

## Rassemblement par 
forest <- c()
data <- c()
leaves <- c()
for (i in length(tab_pval[,1])) {
  index <- length(data[,1])
  data_index <- !is.na(tab_pval[i,])
  data_index[1] <- FALSE
  ID <- paste(rep(c(tab_pval[i,1]),each = length(tab_pval[i,data_index])),vect_couple[data_index])
  p_values <-tab_pval[i,data_index]
  data <- rbind(data,cbind(ID,p_values))
  if (data_index[1,"Pgab"]){
  }
    
  }
}
rm(list = ls())