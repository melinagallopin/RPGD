################################################################################
## Tree
################################################################################
library(ape)
## Get the tree
tree <- read.tree(file = here("Arbres/TreeHyperMai2020.nwk"))
plot(tree)

## Species names
sample_annotation <- read.csv(here("data/samples_annotation.txt"), sep="\t")
correspondances <- unique(sample_annotation[, c("Nom.espece", "ID.espece")])
# Format
correspondances[, "Nom.espece"] <- sub(" ", "_", correspondances[, "Nom.espece"])
# Typos
correspondances[, "Nom.espece"] <- sub("Senocio_coronatus", "Senecio_coronatus", correspondances[, "Nom.espece"])
correspondances[, "Nom.espece"] <- sub("Microthlaspi_perfoliatum", "Microthalspi_perfoliatum", correspondances[, "Nom.espece"])
correspondances[, "Nom.espece"] <- sub("Homalium_kanaliense", "Homalium_kanalense", correspondances[, "Nom.espece"])
correspondances[, "Nom.espece"] <- sub("Homalium_betulifoium", "Homalium_betulifolium", correspondances[, "Nom.espece"])
# Match
tree_data_cor <- match(tree$tip.label, correspondances[, "Nom.espece"])
data_tree_cor <- match(correspondances[, "Nom.espece"], tree$tip.label)
# Species in the tree NOT in data
tree$tip.label[is.na(tree_data_cor)]
# Species in data NOT in the tree
correspondances[is.na(data_tree_cor), "Nom.espece"]

## Format Tree
# Rename species with ID
tree$tip.label <- as.character(correspondances[match(tree$tip.label, correspondances[, "Nom.espece"]), "ID.espece"])
plot(tree)
## this tree contains all species 

## we load the pvalue table
load(here("results/tab_pval_length.RData"))
## get list species in pvalue table
list_species <- sapply(colnames(tab_pval)[-1],function(i) strsplit(i,"_")[[1]][1])
## correct species names Pcle => Prev
list_species <- sub("Pcle","Prev",list_species)
## we drop all species not in the list of the "pvalue table" species 
tree <- drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% list_species)])
plot(tree)

## we get the hierarchy of interest based on the phylo tree
thres=unique(unlist(as.list(vcv(tree))))[order(unique(unlist(as.list(vcv(tree)))))]
HierarTree=lapply(thres,function(i) components(graph_from_adjacency_matrix(vcv(tree)>=i,mode = "undirected"))$membership)
HierarTree
# save hierarchy 
save(HierarTree,file="results/hierarchy.RData")

