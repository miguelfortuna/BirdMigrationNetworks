######################################
### get ecological distance matrix ###
######################################

### load required packages ----
library("betalink")
library("igraph")
library("PCIT")
library("DCG")
library("tidyverse") # use dplyr::select (instead of select) if library MASS is loaded

### print package version ----
packageVersion("betalink") # 2.2.1
packageVersion("igraph") # 1.2.4
packageVersion("PCIT") # 1.5.3
packageVersion("DCG") # 0.9.2
packageVersion("tidyverse") # 1.2.1

### read all links among IBAs ----
data <- read_csv("../working_data/links_84_species.csv")

### get the 84 species names ----
species <- read_csv("../working_data/species_84.csv") %>% select(species_label, species_name)

### store the 84 networks into a list (same order as in "species")----
networks <- list()
for(i in 1: nrow(species)){
  networks[[i]] <- data %>% filter(species_label == species$species_label[[i]]) %>% select(from, to)
  print (species$species_label[[i]])
}

### convert each network in the list from tibble/dataframe to matrix
networks <- sapply(networks, as.matrix)

### convert networks into igraph objects (requirement for betalink to work) ----
networks_igraph <- list()
  for(i in 1: length(networks)){
    networks_igraph[[i]] <- graph_from_edgelist(networks[[i]], directed = TRUE)
}

### add the name of the species to each network stored into the list ----
names(networks_igraph) <- species$species_name

### plot a single network
plot(networks_igraph[["Branta bernicla"]],
     layout = layout.fruchterman.reingold,
     vertex.size = 4,
     vertex.label.dist = 1.25,
     vertex.color = "red",
     vertex.label.color = "black",
     edge.width = 1,
     edge.color = "black",
     edge.arrow.size = 0.05)

### calculate dissimilarity (i.e., beta-diversity) for the list of networks ----
betadiversity <- network_betadiversity(networks_igraph) # it might take 30 minutes
betadiversity <- betadiversity %>% rename(species_1 = i, species_2 = j)

### report OS and ST as % of the total dissimilarity (WN) ----
betadiv_rel <- betadiversity %>%
  mutate(relative_OS = OS/WN, relative_ST = ST/WN) %>%
  select(relative_OS, relative_ST)

### ecological distances as edge list (based on dissimilarity in their migration networks) ----
eco_dist_edges <- cbind(betadiversity, betadiv_rel)

# show rows with missing values (i.e., WN=1 and ST=NaN, and OS=NaN; there are no IBAs in common)
eco_dist_edges %>% filter(!complete.cases(.)) %>% head

# replace NaN values with 0
eco_dist_edges <- eco_dist_edges %>% replace(is.na(.), 0)

### sort alphabetically so that species_1 and species_2 follow the same order as in the phylogenetic distances ----
eco_dist_edges <- eco_dist_edges %>% arrange(species_1, species_2)

### ecological distances as an adjacency matrix ----
# S = dissimilarity in IBAs composition
S <- eco_dist_edges %>% select(species_1, species_2, S)
S$species_1 <- as.character(S$species_1)
S$species_2 <- as.character(S$species_2)
eco_dist_matrix_S <- as.symmetricAdjacencyMatrix(S, weighted = TRUE, rule = "weak")

# WN = total dissimilarity between the two networks ----
WN <- eco_dist_edges %>% select(species_1, species_2, WN)
WN$species_1 <- as.character(WN$species_1)
WN$species_2 <- as.character(WN$species_2)
eco_dist_matrix_WN <- as.symmetricAdjacencyMatrix(WN, weighted = TRUE, rule = "weak")

# ST = dissimilarity due to species turnover
ST <- eco_dist_edges %>% select(species_1, species_2, ST)
ST$species_1 <- as.character(ST$species_1)
ST$species_2 <- as.character(ST$species_2)
eco_dist_matrix_ST <- as.symmetricAdjacencyMatrix(ST, weighted = TRUE, rule = "weak")

# OS = dissimilarity due to rewiring
OS <- eco_dist_edges %>% select(species_1, species_2, OS)
OS$species_1 <- as.character(OS$species_1)
OS$species_2 <- as.character(OS$species_2)
eco_dist_matrix_OS <- as.symmetricAdjacencyMatrix(OS, weighted = TRUE, rule = "weak")

# rel_ST = dissimilarity due to species turnover as a fraction of total dissimilarity
rel_ST <- eco_dist_edges %>% select(species_1, species_2, relative_ST)
rel_ST$species_1 <- as.character(rel_ST$species_1)
rel_ST$species_2 <- as.character(rel_ST$species_2)
eco_dist_matrix_rel_ST <- as.symmetricAdjacencyMatrix(rel_ST, weighted = TRUE, rule = "weak")

# rel_OS = dissimilarity due to rewiring as a fraction of total dissimilarity
rel_OS <- eco_dist_edges %>% select(species_1, species_2, relative_OS)
rel_OS$species_1 <- as.character(rel_OS$species_1)
rel_OS$species_2 <- as.character(rel_OS$species_2)
eco_dist_matrix_rel_OS <- as.symmetricAdjacencyMatrix(rel_OS, weighted = TRUE, rule = "weak")

# save to file ----
#write.csv(eco_dist_edges, "ecological_distance_edgelist.csv", row.names = FALSE)
#write.csv(eco_dist_matrix_S, "ecological_distance_matrix_S.csv")
#write.csv(eco_dist_matrix_WN, "ecological_distance_matrix_WN.csv")
#write.csv(eco_dist_matrix_ST, "ecological_distance_matrix_ST.csv")
#write.csv(eco_dist_matrix_OS, "ecological_distance_matrix_OS.csv")
#write.csv(eco_dist_matrix_rel_ST, "ecological_distance_matrix_rel_ST.csv")
#write.csv(eco_dist_matrix_rel_OS, "ecological_distance_matrix_rel_OS.csv")
