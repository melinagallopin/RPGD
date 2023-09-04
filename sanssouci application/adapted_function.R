organized_p_values<- function(HierarTree,cog_p_values){  ## need to change get_mole_forest to get a re-organized which match with the dataset?
  order_species <- names(HierarTree[[1]])
  swap_p_values <- data.frame()
  for (i in  1:length(HierarTree)){
    swap_p_values <- rbind(swap_p_values,cog_p_values[, grepl( order_species[i] , names(cog_p_values))])
  }
  df <- stack(swap_p_values)
  return(data.frame(df[,1]))
} ## ok considering we want to swap p_values
colnames(tab_pval)[colnames(tab_pval) == "Pcle_3"] <- "Prev"
df <- organized_p_values(HierarTree,tab_pval)
get_full_forest <- function(HierarTree,cog_p_values){
  mole_tree <- get_mole_tree(HierarTree)
  nb_row <- cog_p_values[[1]]
  forest <- mole_tree
  H <- length(mole_tree)
  for (i in 2:nb_row) {
    mole_tree <- following_tree(mole_tree,length(HierarTree))
    for (h in 1:H) {
      forest[[h]] <- append(forest[[h]], mole_tree[[h]])
    }
  }
  return(forest)
}
get_bounds_DKWM <- function(HierarTree,cog_p_values){
  full_forest <- get_full_forest(HierarTree,cog_p_values)
  p_values <- organized_p_values(HierarTree,cog_p_values)
  
}