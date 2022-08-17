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

# Checking the root of the tree
is.rooted(tree)

# The tree is unrooted, but it should be, let us fix this defining who are the outgroups
r_tree <- ape::root(tree,
                    outgroup = c("Nypa_fruticans",
                                 "Pseudophoenix_vinifera",
                                 "Trachycarpus_martianus"),
                    resolve.root = T)
plot.phylo(r_tree)

is.rooted(r_tree)

# Saving the re-rooted tree
writeNexus(r_tree, "./output/rerooted_tree.nex")

# In the ancestral area reconstruction we are going use the trees must be ultra metric
# Let us transform this tree in ultrametric using to different methods
tree.Ult.MPL <- ape::chronoMPL(r_tree) # Mean path length method
tree.Ult.S <- ape::chronos(r_tree) # Maximum likelihood method

is.ultrametric(tree.Ult.MPL)
is.ultrametric(tree.Ult.S)

plot(tree.Ult.MPL) # This one has negative branch values, so it will be disregarded
plot(tree.Ult.S)

# Saving the ultrametric re-rooted tree
rm(tree.Ult.S)
writeNexus(tree.Ult.S, "./output/ultra_rerooted_tree.nex")
