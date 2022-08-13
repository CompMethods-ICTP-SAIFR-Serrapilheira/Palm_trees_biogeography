####### Writing began in August 2022 by Thales
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script imports and inspect the phylogenetic tree used
####### during the rest of the work

## The tree used here was inferred outside R, using other freewares
## I took this path because the tools available to infer phylogenies in R are not the best
## For the documentation on how the phylogenetic tree was inferred can be read by running
## the commented code line bellow
#file.edit('./docs/tree_inference_documentation.Rmd')

# Library
library(phytools)
library(ape)

# Importing the tree
tree <- ape::read.nexus("./output/final_consensus_tree_run1.nex")
ape::plot.phylo(tree, type = "phylogram") # Checking how the tree is

# Seeing the structure of the object tree of 'phylo' class
str(tree)

# Plotting tree with node and tip labels
plot.phylo(tree, type = "phylogram")
nodelabels(cex=1, frame = "none")

# The tree is rooted in the node 40, but it should be rooted in node 45, let us fix this
r_tree <- ape::root(tree, node = 45)
plot.phylo(r_tree)

# Saving the re-rooted tree
writeNexus(r_tree, "./output/rerooted_tree.nex")

# Importing the re-rooted tree
rm(r_tree)
tree <- read.nexus("./output/rerooted_tree.nex")
plot.phylo(tree)
