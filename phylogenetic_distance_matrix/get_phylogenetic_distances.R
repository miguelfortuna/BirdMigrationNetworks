########################################
### get phylogenetic distance matrix ###
########################################

### load required packages ----
library("motmot.2.0")
library("ape")
library("geiger")
library("phytools")
library("PCIT")
library("tidyverse") # use dplyr::select (instead of select) if library MASS is loaded

### print package version ----
packageVersion("motmot.2.0") # 1.1.2
packageVersion("ape") # 5.3
packageVersion("geiger") # 2.0.6.1
packageVersion("phytools") # 0.6.60
packageVersion("PCIT") # 1.5.3
packageVersion("tidyverse") # 1.2.1

### read phylogeny (no branch lengths) ----
tree <- read.tree("phylogeny.tre")


###
### set arbitrary branch lengths resulting in different levels of hierarchy ----
###

# Grafen ----
# branch lengths set such that the depth of a node is equal to the number of descendant tips minus one
n <- length(tree$tip.label)
depth <- node.depth(tree, method = 1) - 1
tree$edge.length <- (depth[tree$edge[, 1]]/(n - 1)) - (depth[tree$edge[, 2]]/(n - 1))
grafen <- tree
plot(grafen, type="fan", no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)
#plot(grafen, no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)

# Pagel ----
# branch lengths set such that the depth of a node is equal to the maximum of the number of bifurcation levels of the left child and the right child plus one
n <- length(tree$tip.label)
depth <- node.depth(tree, method = 2)
tree$edge.length <- (depth[tree$edge[, 1]]/(n - 1)) - (depth[tree$edge[, 2]]/(n - 1))
pagel <- tree
plot(pagel, type="fan", no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)
#plot(pagel, no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)

# Nee ----
# branch lengths set such that the depth of a node is equal to the log 10 of the number of tips descending from that node
n <- length(tree$tip.label)
depth <- log(node.depth(tree, method = 1), 10)
tree$edge.length <- (depth[tree$edge[, 1]]/(n - 1)) - (depth[tree$edge[, 2]]/(n - 1))
nee <- tree
plot(nee, type="fan", no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)
#plot(nee, no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)

### decompose the phylogeny to its essential components to show the branch transformations ----
#phy <- grafen
phy <- pagel
#phy <- nee
plot(phy, no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)
tiplabels(cex = 0.5, frame="circle", width = 3, offset = 1)
nodelabels(cex = 0.5, frame = "circle")


###
### find the branch lengths that best fit the distribution of body mass (log) ----
###

### read bodysize of the 84 bird species ----
bodysize <- read_csv("body_mass.csv") %>%
  dplyr::select(species_name, log_mean) %>%
  rename(size = log_mean) %>%
  column_to_rownames("species_name")
bodysize <- matrix(bodysize$size, dimnames=list(rownames(bodysize)))

### sort bodysize data to match phylogeny ----
bodysize <- as.matrix(bodysize[match(tree$tip.label, rownames(bodysize)), ])
### plot our tree and data showing species body mass (log) at the tips ----
#startingPhylo <- grafen
startingPhylo <- pagel
#startingPhylo <- nee
traitData.plot(phy = startingPhylo, y = bodysize)

### compute maximum likelihood (log-likelihood) for different values of the Ornstein-Uhlenbeck transformation parameter (alpha) ---
# starting from the three initial branch lengths (Grafen, Pagel, and Nee)
# check if bifurcating tree
# is.binary.tree(startingPhylo) # FALSE (it contains polytomies)
# ou.ml <- transformPhylo.ML(phy = startingPhylo, y = bodysize, model = "OU", returnPhy = TRUE)
# it does not work because function "pic" (from geiger function "transformPhylo.ML") needs a binary tree as input
# "multi2di" transforms polytomies into a series of dichotomies with one (or several) branch(es) of length zero
# plot likelihood profile for the branch-transformation parameter, in this case OU's alpha ----

# Grafen ----
### sort bodysize data to match phylogeny ----
bodysize <- as.matrix(bodysize[match(grafen$tip.label, rownames(bodysize)), ])
startingPhylo <- grafen
ou.ml <- transformPhylo.ML(
  phy = multi2di(startingPhylo),
  y = bodysize, model = "OU",
  returnPhy = TRUE,
  profilePlot = TRUE,
  lowerBound = 0.1,
  upperBound = 10)
# get value of the parameter alpha that maximizes log-likelihood (and hence, minimizes MSE)
ou.ml$Alpha[1] # 4.044429
# get maximum log-likelihood
ou.ml$MaximumLikelihood # -35.52056
# get the phylogeny with branch lengths transformed by the ML model parameters and plot it
transformedGrafen <- ou.ml$ouPhy
# plot the phylogeny with branch lengths transformed and showing species body mass (log) at the tips
traitData.plot(phy = transformedGrafen, y = bodysize)

# Pagel ----
### sort bodysize data to match phylogeny ----
bodysize <- as.matrix(bodysize[match(pagel$tip.label, rownames(bodysize)), ])
startingPhylo <- pagel
ou.ml <- transformPhylo.ML(
  phy = multi2di(startingPhylo),
  y = bodysize, model = "OU",
  returnPhy = TRUE,
  profilePlot = TRUE,
  lowerBound = 0.1,
  upperBound = 5)
# get value of the parameter alpha that maximizes log-likelihood (and hence, minimizes MSE)
ou.ml$Alpha[1] # 2.113702
# get maximum log-likelihood
ou.ml$MaximumLikelihood # -31.97973
# get the phylogeny with branch lengths transformed by the ML model parameters and plot it
transformedPagel <- ou.ml$ouPhy
# plot the phylogeny with branch lengths transformed and showing species body mass (log) at the tips
traitData.plot(phy = transformedPagel, y = bodysize)

# Nee ----
### sort bodysize data to match phylogeny ----
bodysize <- as.matrix(bodysize[match(nee$tip.label, rownames(bodysize)), ])
startingPhylo <- nee
ou.ml <- transformPhylo.ML(
  phy = multi2di(startingPhylo),
  y = bodysize, model = "OU",
  returnPhy = TRUE,
  profilePlot = TRUE,
  lowerBound = 0.001,
  upperBound = 10)
# get value of the parameter alpha that maximizes log-likelihood (and hence, minimizes MSE)
ou.ml$Alpha[1] # 0.001
# get maximum log-likelihood
ou.ml$MaximumLikelihood # -34.31021
# get the phylogeny with branch lengths transformed by the ML model parameters and plot it
transformedNee <- ou.ml$ouPhy
# plot the phylogeny with branch lengths transformed and showing species body mass (log) at the tips
traitData.plot(phy = transformedNee, y = bodysize)

### best phylogeny ----
#the branch lengths that best fit the distribution of log body mass (maximum log-likelihood = -31.97973)
# is the OU-transformation with alpha = 2.113702 of the initial tree obtained by the Pagel algorithm
# plot
traitData.plot(phy = transformedPagel, y = bodysize)
plot(transformedPagel, type="fan", no.margin = TRUE, edge.width = 1, show.tip.label = FALSE)

### phylogenetic distances between species represented by the phylogenetic variance-covariance matrix ----
phy_dist <- cophenetic(transformedPagel)

### sort alphabetically so that rownames and colnames follow the same order as in the ecological distances ----
phy_dist <- phy_dist[sort(rownames(phy_dist)), sort(colnames(phy_dist))]

### phylogenetic distances as an edge list ----
phy_dist_edges <- getEdgeList(phy_dist, rm.zero = FALSE) # getEdgeList assumes the matrix is symmetric
phy_dist_edges <- phy_dist_edges %>% rename(species_1 = From, species_2 = To, distance = Weight)

### save output to file ----
#write.tree(phy = transformedPagel, "transformedPagel.tre")
#write.csv(phy_dist, "phylogenetic_distance_matrix.csv")
#write.csv(phy_dist_edges, "phylogenetic_distance_edgelist.csv", row.names = FALSE)


### plot phylogeny ----
# read transformed phylogeny
tree <- read.tree("transformedPagel.tre")
plot(tree, type = "fan",
     label.offset = 0.01,
     edge.width = 1,
     edge.lty = 1, 
     cex = 0.75)

### save phylogeny ----
svg(filename = "fig_3_phylo.svg", 
    width = 32, height = 32, pointsize = 36)
plot(tree, type = "fan", label.offset = 0.005, edge.width = 1, edge.lty = 1, cex = 0.75)
dev.off()