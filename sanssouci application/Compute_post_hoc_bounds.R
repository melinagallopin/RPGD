load("results/hierarchy.RData")
load("results/tab_pval_length.RData")
library(sanssouci)

## Debug functions ---------------

zetas.tree.no.extension <- function(C, leaf_list, method, pvalues, alpha, refine=FALSE, verbose=FALSE) {
  H <- length(C)
  K <- nb.elements.no.extension(C)
  ZL <- list()
  new_K <- K
  continue <- TRUE
  nb_loop <- 0
  while (continue) {
    usage_K <- new_K
    new_K <- K
    for (h in H:1) {
      Ch <- C[[h]]
      len <- length(Ch)
      zeta_inter <- numeric(len)
      for (j in 1:len) {
        Chj <- Ch[[j]]
        pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
        if (typeof(method) == "list") {
          if (typeof(method[[h]]) == "list") {
            zeta_method <- method[[h]][[j]]
          } else {
            zeta_method <- method[[h]]
          }
        } else {
          zeta_method <- method
        }
        zeta_inter[j] <- zeta_method(pvals, alpha / usage_K)
        if (refine && (zeta_inter[j] == 0) )
          new_K <- new_K - 1
      }
      ZL[[h]] <- zeta_inter
    }
    if (verbose) {
      nb_loop <- nb_loop + 1
      print(paste0("loop number=", nb_loop,", usage_K=",usage_K,", new_K=",new_K))
    }
    continue <- refine && (new_K < usage_K)
  }
  return(ZL)
}

nb.elements.no.extension <- function(C) {
  H <- length(C)
  count <- 0
  for (h in H:1) {
    count <- count + length(C[[h]])
  }
  return(count)
}

zeta.HB.no.extension <- function(pval, lambda) {
  m <- length(pval)
  sorted.pval <- sort(pval)
  
  thresholds <- lambda / (m:1)
  v <- sorted.pval - thresholds
  indexes <- which(v > 0)
  if (! length(indexes)) {
    return(0)
  }
  else{
    return(m - indexes[1] + 1)
  }
}

V.star.no.extension <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves) 
  for (i in 1:nb_leaves) {
    Vec[i] <- sum(S %in% leaf_list[[i]])
  }
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in C and ZL
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions>0) {
      for (j in 1:nb_regions) {
        Chj <- C[[h]][[j]]
        if (Chj[1]==Chj[2]) { # means this is an atom, no need to compute 
          # len_inter given that we already did it during initialization,
          # furthermore there are no successors
          len_inter <- Vec[Chj[1]]
          res <- min(ZL[[h]][j], len_inter)
        } else {
          region_vector <- unlist(leaf_list[Chj[1]:Chj[2]])
          len_inter <- sum(S %in% region_vector)
          sum_succ <- sum(Vec[Chj[1]:Chj[2]]) 
          res <- min(ZL[[h]][j], len_inter, sum_succ)
        }
        Vec[Chj[1]:Chj[2]] <- 0
        Vec[Chj[1]] <- res
      }
    }
  }
  return(sum(Vec))
}

## Functions for tree usage ------------

is_in <- function(x,y) {
  if (y[1]<=x[1] && y[2]>=x[length(x)] ){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

get_mole_tree_leaves <- function(HierarTree,cog_p_values) {
  mole_tree <- list(list(c(1,length(HierarTree))))
  mole_leaves <- list(c(1,length(HierarTree)))
  for (i in 2:length(HierarTree)){
    for (j in 1:i){
      all_p_val <- which(unname(HierarTree[[i]])==j)
      region <- c(all_p_val[1],all_p_val[length(all_p_val)])
      H <- length(mole_tree)
      h <- H
      continue <- TRUE
      if (region[1] < region[2]) {
        while (continue && h !=0){
          treeh <- mole_tree[[h]]
          for (k in 1:length(treeh)) {
            treehk <- treeh[[k]]
            if (treehk[1] <= region[1] && treehk[2] >= region[2]){
              if (treehk[1] != region[1] || treehk[2] != region[2]) {
                if (h == H){
                  mole_tree <- append(mole_tree, list(list(region)))
                }
                else {
                  mole_tree[[h+1]] <- append(mole_tree[[h+1]], list(region))
                }
                ind <- which(unlist(lapply(mole_leaves, identical, y = treehk)))
                if (identical(ind, integer(0))){
                  mole_leaves <- append(mole_leaves,list(region))
                }
                else {
                  mole_leaves[[ind]]<-region
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
      else {
        if (sum(sapply(unlist(mole_leaves), identical, y = region[1])) == 0){
          mole_leaves <- append(mole_leaves,list(region[1]))
        }
      }
    }
  }
  mole_leaves <- mole_leaves[order(sapply(mole_leaves,function(x){x[1]}))]
  mole_tree_bis <- mole_tree
  for (h in 1:length(mole_tree_bis)){
    for (i in 1:length(mole_tree_bis[[h]])){
      indices <- which(unlist(lapply(mole_leaves,is_in,y = mole_tree_bis[[h]][[i]])))
      mole_tree_bis[[h]][[i]] <- c(indices[1],indices[length(indices)])
    }
  }
  return(list(mole_tree_bis, mole_leaves))
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
  return(new_tree)
}

next_leaves <- function(leaves) {
  nb_p_values <- length(unlist(leaves))
  new_leaves <- lapply(leaves,function(vect){sapply(vect,"+",y = nb_p_values)})
  return(new_leaves)
}

colnames(tab_pval)[colnames(tab_pval) == "Pcle_3"] <- "Prev"
tab_pval_no_na <- tab_pval[!rowSums(is.na(tab_pval))>0,]

get_full_forest <- function(HierarTree,cog_p_values){
  mole_tree_leaves <- get_mole_tree_leaves(HierarTree)
  mole_tree <- mole_tree_leaves[[1]]
  order_species <- names(HierarTree[[1]])
  mole_leaves <- lapply(mole_tree_leaves[[2]],function(x){
    for (i in 1:length(x)){x <- replace(x,i, which(grepl( order_species[x[i]] , colnames(cog_p_values[2:length(cog_p_values)]))))
    }
  return(x)})
  nb_row <- length(cog_p_values[[1]])
  forest <- mole_tree
  leaves <- mole_leaves
  H <- length(mole_tree)
  for (i in 2:nb_row) {
    mole_tree <- next_tree(mole_tree)
    mole_leaves <- next_leaves(mole_leaves)
    leaves <- append(leaves,mole_leaves)
    for (h in 1:H) {
      forest[[h]] <- append(forest[[h]], mole_tree[[h]])
    }
  }
  return(list(forest,leaves))
}

get_zetas_phylotree_no_na <- function(HierarTree,cog_p_values,method,alpha){
  forest_leaves <- get_full_forest(HierarTree,cog_p_values)
  forest <- forest_leaves[[1]]
  leaves <- forest_leaves[[2]]
  reduced_tab <- cog_p_values[2:length(cog_p_values)]
  p_values <- t(data.frame(p_values=matrix(t(reduced_tab))))
  bounds <- zetas.tree.no.extension(forest,leaves,method,p_values,alpha)
  return(bounds)
}

get_Vstar_cogs_phylotree_no_na <- function(HierarTree,cog_p_values,method,alpha){
  forest_leaves <- get_full_forest(HierarTree,cog_p_values)
  forest <- forest_leaves[[1]]
  leaves <- forest_leaves[[2]]
  reduced_tab <- cog_p_values[2:length(cog_p_values)]
  p_values <- t(data.frame(p_values=matrix(t(reduced_tab))))
  zetas <- zetas.tree.no.extension(forest,leaves,method,p_values,alpha)
  cogs <- lapply(forest[[1]],function(x){return(unlist(leaves[x[1]:x[2]]))})
  Vstar_cog <- lapply(cogs,function(S){print(S[1]) 
    zeta <-V.star.no.extension(S,forest,zetas,leaves)
    print(zeta)
    return(zeta)})
  Vstar_cog <- unlist(Vstar_cog)
  nb_diff <- sum(Vstar_cog != zetas[[1]])
  return(list(data.frame(cbind(cog_p_values[1],Vstar_cog)),nb_diff)) ## delete nb_diff for usage
}

## Functions for HB global procedure -------

HB_procedure <- function(p_values, alpha) {
  order_p_val <- order(p_values)
  ordered_p_values <- p_values[order_p_val]
  tresholds <- alpha/length(p_values):1
  indices <- which(ordered_p_values-tresholds > 0)
  if (!length(indices)){
    return(order_p_val)
  }
  else{
    return(order_p_val[1:(indices[1]-1)])
  }
}

global_HB <- function(cog_p_values,alpha){
  p_values <- t(data.frame(p_values=matrix(t(cog_p_values[colnames(cog_p_values)[-1]]))))
  indices_HB <- HB_procedure(p_values,alpha)
  nb_p_values_by_cog <- length(cog_p_values)-1
  nb_cog <- length(cog_p_values[,1])
  bound_false_positives_by_cog <- rep(nb_p_values_by_cog,each = nb_cog) 
  for (index in indices_HB){
    bound_false_positives_by_cog[index %/% nb_p_values_by_cog +1] <- bound_false_positives_by_cog[index %/% nb_p_values_by_cog +1] - 1
  }
  return(cbind.data.frame(cog_p_values[,1],bound_false_positives_by_cog))
}

## Naive refinement -------------

temp_0 <- sum(unlist(get_zetas_phylotree_no_na(HierarTree,tab_pval_no_na,zeta.HB.no.extension,0.05))==0)
number_of_regions_at_0 <- temp_0
K <- nb.elements.no.extension(get_full_forest(HierarTree,tab_pval_no_na))
while (temp_0 != 0L) {
temp_0 <- sum(unlist(get_zetas_phylotree_no_na(HierarTree,tab_pval_no_na,zeta.HB.no.extension,0.05*(K/(K - number_of_regions_at_0))))==0) - number_of_regions_at_0
number_of_regions_at_0 <- number_of_regions_at_0 + temp_0
print(temp_0)
}
zetas_refinement <- get_zetas_phylotree_no_na(HierarTree,tab_pval_no_na,zeta.HB.no.extension,0.05*(K/(K - number_of_regions_at_0)))
zetas <- get_zetas_phylotree_no_na(HierarTree,tab_pval_no_na,zeta.HB.no.extension,0.05)
nb_change <- sum(unlist(zetas) - unlist(zetas_refinement) != 0)

## Simes --------------------

simes_procedure <- function(cog_p_values,alpha){
  p_values <- t(data.frame(p_values=matrix(t(cog_p_values[colnames(cog_p_values)[-1]]))))
  nb_p_values_by_cog <- length(cog_p_values)-1
  nb_cog <- length(cog_p_values[,1])
  bound_false_positives_by_cog <- rep(nb_p_values_by_cog,each = nb_cog)
  for (i in 1:nb_cog){
    bound_false_positives_by_cog[i] <- posthocBySimes(p_values,((i-1)*nb_p_values_by_cog +1):(i*nb_p_values_by_cog +1),alpha,Rcpp = TRUE)
    print(i)
    print(((i-1)*nb_p_values_by_cog +1):(i*nb_p_values_by_cog +1))
  }
  return(cbind.data.frame(COG_ID = cog_p_values[,1],Bound_Simes = bound_false_positives_by_cog))
}

## Cherry ----------

library(cherry)
library(here)

par(mfrow=c(2,2))
for(i in 2:8){
  hist(tab_pval[,i])
}

reduced_tab <- tab_pval_no_na[2:length(tab_pval_no_na)]
p_values <- t(data.frame(p_values=matrix(t(reduced_tab))))

#correction d'hommel
hom <- hommelFast(as.numeric(p_values),simes=TRUE)

#donne les vrais positifs
cogs_H1 <- sapply(tab_pval_no_na$COG_ID, FUN=function(i) pickSimes(hom, select=tab_pval_no_na$COG_ID == i))

# cogs_H1_list <- list()
# for(ut in tab_pval$COG_ID){
#   print(ut)
#   cogs_H1_list[ut] <- pickSimes(hom, select=tab_pval_no_na$COG_ID == ut)
# }
# 
# tab_pval_no_na[tab_pval_no_na$COG_ID == "3564",]

#tableau de résultat
bound_tp_cherry <- data.frame(COG_ID=tab_pval$COG_ID,  TP_pred=cogs_H1)


## Tests --------------------

results <- get_Vstar_cogs_phylotree_no_na(HierarTree,tab_pval_no_na,zeta.HB.no.extension,0.05)

results_bis <- results[[1]]

results_HB_glob <- global_HB(tab_pval_no_na,0.05)

nb_HB_glob_better <- sum(results_bis[,2] - results_HB_glob[,2] > 0)

nb_HB_loc_better <- sum(results_bis[,2] - results_HB_glob [,2] < 0)

nb_cherry_better_than_hb_tree <- sum(results_bis[,2] + bound_tp_cherry[,2] > 7)

nb_cherry_better_than_hb_glob <- sum(bound_tp_cherry[,2] + results_HB_glob[,2] > 7)

nb_hb_tree_better_than_cherry <- sum(results_bis[,2] + bound_tp_cherry[,2] < 7)

nb_hb_glob_better_than_cherry <- sum(bound_tp_cherry[,2] + results_HB_glob[,2] < 7)

nb_positives_HB_glob <- 7*5895 - sum(results_HB_glob[,2])

nb_positives_HB_loc <- 7*5895 - sum(results_bis[,2])


m <- 10000
m1 <- 200
p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
R <- union(1:10, sample(m, 10))
alpha <- 0.10
if (require("cherry")) {
  hom <- hommelFast(p)
  pickSimes(hom, R, silent=TRUE, alpha = alpha)
}
posthocBySimes(p, R, alpha=alpha)
