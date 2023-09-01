load("results/hierarchy.RData")
load("results/tab_pval_length.RData")

##Forest Structure ----

## Rassemblement par 
forest <- c()
data <- c()
leaves <- c()
for (i in 1:length(tab_pval[,1])) {
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

## Création de tab_pval sans NA (correspond à une paire d'espèce ne possèdant pas le cog associé).----

nona_tab_pval <- tab_pval[!rowSums(is.na(tab_pval))>0,]

## Mini-Structure d'arbre -----

get_mole_forest <- function(HierarTree) {
  mole_forest <- list(list(c(1,length(HierarTree))))
  for (i in 2:length(HierarTree)){
    for (j in 1:i){
      all_p_val <- which(unname(HierarTree[[i]])==j)
      region <- c(all_p_val[1],all_p_val[length(all_p_val)])
      H <- length(mole_forest)
      h <- H
      continue <- TRUE
      while (continue && h !=0){
        foresth <- mole_forest[[h]]
        for (k in 1:length(foresth)) {
          x <- foresth[[k]]
          if ((x[1] <= region[1]) && (x[2] >= region[2])){
            if (x[1] != region[1] || x[2] != region[2]) {
              if (h == H){
                mole_forest <- append(mole_forest, list(list(region)))
              }
                else {
                  mole_forest[[h+1]] <- append(mole_forest[[h+1]], list(region))
                }
              }
              continue <- FALSE
          }
          if (continue == FALSE){
            break
          }
        }
        h <- h-1
      }
    }
  }
}

next_tree <- function(tree) {
  nb_leaves_tree <- tree[[1]][[1]][2] - tree[[1]][[1]][1] +1
  new_tree <- tree
  for (h in 1:length(tree)) {
    treeh <- tree[[h]]
    for (i in 1:length(treeh)){
      treehi <- treeh[[i]]
      new_tree[[h]][[i]][1] <- treehi[1] + nb_leaves_tree
      new_tree[[h]][[i]][2] <- treehi[2] + nb_leaves_tree
    }
  }
}

get_forest <- function(HierarTree,data){
  mole_forest <- get_mole_forest(HierarTree)
  new_tree <- mole_forest
  H <- length(mole_forest)
  for (i in 1:(length(data)-1)){
    new_tree <- next_tree(new_tree)
    for (h in 1:H) {
      mole_forest[[h]] <- append(mole_forest[[h]],new_tree[[h]]) 
    }
  }
}

get_p_values <- function(HierarTree,data){
  ord <- colnames(HierarTree)
}

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

rm(list = ls())

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


