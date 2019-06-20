#########################################################################
### correlation between phylogenetic and ecological distance matrices ###
#########################################################################

### load required packages ----
library("ape")
library("ade4")
library("energy")
library("tidyverse") # use dplyr::select (instead of select) if library MASS is loaded

### print package version ----
packageVersion("ape") # 5.3
packageVersion("ade4") # 1.7.13
packageVersion("energy") # 1.7.5
packageVersion("tidyverse") # 1.2.1

### read phylogenetic distance_matrix ----
phylo_dist <- read_csv("../phylogenetic_distance_matrix/phylogenetic_distance_matrix.csv") # species names are separated by hyphen

### read ecological distance_matrices ----
eco_dist_WN <- read_csv("../ecological_distance_matrix/ecological_distance_matrix_WN.csv") # species names are separated by blank space
eco_dist_ST <- read_csv("../ecological_distance_matrix/ecological_distance_matrix_ST.csv") # species names are separated by blank space
eco_dist_OS <- read_csv("../ecological_distance_matrix/ecological_distance_matrix_OS.csv") # species names are separated by blank space
eco_dist_rel_ST <- read_csv("../ecological_distance_matrix/ecological_distance_matrix_rel_ST.csv") # species names are separated by blank space
eco_dist_rel_OS <- read_csv("../ecological_distance_matrix/ecological_distance_matrix_rel_OS.csv") # species names are separated by blank space
eco_dist_S <- read_csv("../ecological_distance_matrix/ecological_distance_matrix_S.csv") # species names are separated by blank space


### Distance covariance/correlation between the ecological and phylogenetic matrices ----

# phylogenetic distance matrix
phylo_dist_matrix <- as.matrix(phylo_dist %>% select(-X1))
phylo_dist_matrix <- as.dist(phylo_dist_matrix)
# Euclidean properties
summary(phylo_dist_matrix) # using "ape" package
is.euclid(phylo_dist_matrix, print = TRUE) # the matrix is Euclidean by the Gower's theorem (i.e, all the eigenvalues are positive or equal to zero)

# WN
eco_dist_matrix_WN <- as.matrix(eco_dist_WN %>% select(-X1))
eco_dist_matrix_WN <- as.dist(eco_dist_matrix_WN)
summary(eco_dist_matrix_WN) # using "ape" package
is.euclid(eco_dist_matrix_WN, print = TRUE) # the matrix is Euclidean by the Gower's theorem (i.e, all the eigenvalues are positive or equal to zero)
# distance covariance/correlation (using "energy" package) 
dcov.test(phylo_dist_matrix, eco_dist_matrix_WN, R = 1000)
# nV^2 = 0.21098; p-value = 0.000999
dcor.test(phylo_dist_matrix, eco_dist_matrix_WN, R = 1000)
# dCor = 0.6976; p-value = 0.000999
# Mantel test (using "ade4" package)
r_WN <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_WN)
# Pearson correlation = 0.2806331
# p-value < 0.001
plot(r_WN <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_WN))
#the observed correlation suggests that the matrix entries are positively associated

# ST
eco_dist_matrix_ST <- as.matrix(eco_dist_ST %>% select(-X1))
eco_dist_matrix_ST <- as.dist(eco_dist_matrix_ST)
summary(eco_dist_matrix_ST) # using "ape" package
sum(eco_dist_matrix_ST == 0) # 880 zero values
# replace them by 0.0000001
eco_dist_matrix_ST[eco_dist_matrix_ST == 0] <- 0.0000001
is.euclid(eco_dist_matrix_ST, print = TRUE) # the matrix is not Euclidean
# make the matrix Euclidean (using "ade4" package)
eco_dist_matrix_ST <- quasieuclid(eco_dist_matrix_ST)
is.euclid(eco_dist_matrix_ST, print = TRUE) # the matrix is now Euclidean
# distance covariance/correlation (using "energy" package) 
dcov.test(phylo_dist_matrix, eco_dist_matrix_ST, R = 1000)
# nV^2 = 0.068832; p-value = 0.1269
dcor.test(phylo_dist_matrix, eco_dist_matrix_ST, R = 1000)
# dCor = 0.55976; p-value = 0.1349
# Mantel test (using "ade4" package)
r_ST <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_ST) # using "ade4" package
# correlation = 0.05031181
# p-value = 0.055
plot(r_ST <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_ST))
#the observed correlation suggests that the matrix entries are not associated

# OS
eco_dist_matrix_OS <- as.matrix(eco_dist_OS %>% select(-X1))
eco_dist_matrix_OS <- as.dist(eco_dist_matrix_OS)
summary(eco_dist_matrix_OS) # using "ape" package
sum(eco_dist_matrix_OS == 0) # 155 zero values
# replace them by 0.0000001
eco_dist_matrix_OS[eco_dist_matrix_OS == 0] <- 0.0000001
is.euclid(eco_dist_matrix_OS, print = TRUE) # the matrix is not Euclidean
# make the matrix Euclidean (using "ade4" package)
eco_dist_matrix_OS <- quasieuclid(eco_dist_matrix_OS)
is.euclid(eco_dist_matrix_OS, print = TRUE) # the matrix is now Euclidean
# distance covariance/correlation (using "energy" package) 
dcov.test(phylo_dist_matrix, eco_dist_matrix_OS, R = 1000)
# nV^2 = 0.23612; p-value = 0.000999
dcor.test(phylo_dist_matrix, eco_dist_matrix_OS, R = 1000)
# dCor = 0.66074; p-value = 0.000999
# Mantel test (using "ade4" package)
r_OS <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_OS) # using "ade4" package
# correlation = 0.09858485
# p-value = 0.001
plot(r_OS <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_OS))
#the observed correlation suggests that the matrix entries are positively associated

# rel_ST
eco_dist_matrix_rel_ST <- as.matrix(eco_dist_rel_ST %>% select(-X1))
eco_dist_matrix_rel_ST <- as.dist(eco_dist_matrix_rel_ST)
summary(eco_dist_matrix_rel_ST) # using "ape" package
sum(eco_dist_matrix_rel_ST == 0) # 880 zero values
# replace them by 0.0000001
eco_dist_matrix_rel_ST[eco_dist_matrix_rel_ST == 0] <- 0.0000001
is.euclid(eco_dist_matrix_rel_ST, print = TRUE) # the matrix is not Euclidean
# make the matrix Euclidean (using "ade4" package)
eco_dist_matrix_rel_ST <- quasieuclid(eco_dist_matrix_rel_ST)
is.euclid(eco_dist_matrix_rel_ST, print = TRUE) # the matrix is now Euclidean
# distance covariance/correlation (using "energy" package) 
dcov.test(phylo_dist_matrix, eco_dist_matrix_rel_ST, R = 1000)
# nV^2 = 0.069879; p-value = 0.1868
dcor.test(phylo_dist_matrix, eco_dist_matrix_rel_ST, R = 1000)
# dCor = 0.55936; p-value = 0.1818
# Mantel test (using "ade4" package)
r_rel_ST <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_rel_ST) # using "ade4" package
# correlation = 0.04866855
# p-value = 0.054
plot(r_rel_ST <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_rel_ST))
#the observed correlation suggests that the matrix entries are not associated

# rel_OS
eco_dist_matrix_rel_OS <- as.matrix(eco_dist_rel_OS %>% select(-X1))
eco_dist_matrix_rel_OS <- as.dist(eco_dist_matrix_rel_OS)
summary(eco_dist_matrix_rel_OS) # using "ape" package
sum(eco_dist_matrix_rel_OS == 0) # 155 zero values
# replace them by 0.0000001
eco_dist_matrix_rel_OS[eco_dist_matrix_rel_OS == 0] <- 0.0000001
is.euclid(eco_dist_matrix_rel_OS, print = TRUE) # the matrix is not Euclidean
# make the matrix Euclidean (using "ade4" package)
eco_dist_matrix_rel_OS <- quasieuclid(eco_dist_matrix_rel_OS)
is.euclid(eco_dist_matrix_rel_OS, print = TRUE) # the matrix is now Euclidean
# distance covariance/correlation (using "energy" package) 
dcov.test(phylo_dist_matrix, eco_dist_matrix_rel_OS, R = 1000)
# nV^2 = 0.2224; p-value = 0.000999
dcor.test(phylo_dist_matrix, eco_dist_matrix_rel_OS, R = 1000)
# dCor = 0.64512; p-value = 0.000999
# Mantel test (using "ade4" package)
r_rel_OS <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_rel_OS) # using "ade4" package
# correlation = 0.06817952
# p-value = 0.01
plot(r_rel_OS <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_rel_OS))
#the observed correlation suggests that the matrix entries are positively associated

# S
eco_dist_matrix_S <- as.matrix(eco_dist_S %>% select(-X1))
eco_dist_matrix_S <- as.dist(eco_dist_matrix_S)
summary(eco_dist_matrix_S) # using "ape" package
is.euclid(eco_dist_matrix_S, print = TRUE) # the matrix is not Euclidean
# make the matrix Euclidean (using "ade4" package)
eco_dist_matrix_S <- quasieuclid(eco_dist_matrix_S)
is.euclid(eco_dist_matrix_S, print = TRUE) # the matrix is now Euclidean
# distance covariance/correlation (using "energy" package) 
dcov.test(phylo_dist_matrix, eco_dist_matrix_S, R = 1000)
# nV^2 = 0.21323; p-value = 0.000999
dcor.test(phylo_dist_matrix, eco_dist_matrix_S, R = 1000)
# dCor = 0.68018; p-value = 0.000999
# Mantel test (using "ade4" package)
r_S <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_S) # using "ade4" package
# correlation = 0.1837736
# p-value = < 0.001
plot(r_S <- mantel.randtest(phylo_dist_matrix, eco_dist_matrix_S))
#the observed correlation suggests that the matrix entries are positively associated

### Partial distance correlation and covariance (using "energy" package) 
# it measures the association of two random variables X, Y , controlling for a third random variable Z
pdcov.test(phylo_dist_matrix, eco_dist_matrix_OS, eco_dist_matrix_S, R = 1000)
# nV^2  = 0.023177; p-value = 0.000999
pdcor.test(phylo_dist_matrix, eco_dist_matrix_OS, eco_dist_matrix_S, R = 1000)
# dCor = 0.094525; p-value = 0.000999
#the observed correlation suggests that the matrix entries are positively associated
#even when controlling for differences in geographic ranges among species
