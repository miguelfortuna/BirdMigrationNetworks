### load required libraries ----
library("igraph")
library("mapproj")
library("tidyverse")

### print package version ----
packageVersion("igraph") # 1.2.4
packageVersion("mapproj") # 1.2.6
packageVersion("tidyverse") # 1.2.1

### read data ----
data <- read_csv("../working_data/links_84_species.csv")
species <- read_csv("../working_data/species_84.csv") %>% select(species_label, species_name)
ibas <- read_csv("../working_data/ibas_84_species.csv")

### function to compute the degree, strength and betweenness of the nodes of a network ----
comp_node_traits <- function(speciesName) {
  speciesName <- speciesName
  speciesNetwork <- str_replace(speciesName, "_", " ")
  # read data
  web <- inner_join(data, species, by = "species_label") %>%
    filter(species_name == speciesNetwork) %>%
    select(from, to, individuals) %>% rename(weight = individuals)
  # from dataframe to igraph object
  web_igraph <- graph_from_data_frame(web, directed = TRUE)
  # compute degree and strength (as vectors)
  degree_in <- degree(web_igraph, V(web_igraph), mode = "in")
  strength_in <- strength(web_igraph, V(web_igraph), mode = "in")
  # from vector to dataframe
  degree_in <- rownames_to_column(as.data.frame(degree_in)) %>% rename("iba_label" = "rowname")
  strength_in <- rownames_to_column(as.data.frame(strength_in)) %>% rename("iba_label" = "rowname")
  # compute betweenness
  betweenness <- as.data.frame(betweenness(web_igraph, V(web_igraph), directed = TRUE))
  betweenness <- rownames_to_column(betweenness)
  colnames(betweenness) <- c("iba_label", "betweenness")
  df <- inner_join(
    inner_join(
      degree_in,
      strength_in,
      by = "iba_label"),
    betweenness,
    by = "iba_label") %>%
    mutate(log.degree_in = log(degree_in + 1),
           log.strength_in = log(strength_in + 1),
           log.betweenness = log(betweenness + 1),
           species_name = speciesName)
  df
  write_csv(df, paste0("84_species_networks/", speciesName, ".csv"))
}

### compute degree, strength and betweenness of the nodes of the 84 species ----
data_list <- list()
dir.create("84_species_networks")
species_84 <- species %>% distinct(species_name)
for(i in 1:nrow(species_84)) {
  speciesN <- str_replace(species_84$species_name[i], " ", "_")
  #comp_node_traits(speciesN)
  data_list[[i]] <- comp_node_traits(speciesN)
  print(speciesN)
}
data_84_species <- bind_rows(data_list)

### plots to explore correlations ----
plot(degree_in ~ betweenness, data = data_84_species)
plot(strength_in ~ betweenness, data = data_84_species)
plot(strength_in ~ degree_in, data = data_84_species)

# PCA analysis after centring and scaling the predictor variables ----
# (i.e., center variables to have mean zero and scale them to have standard deviation one)
PCA <- prcomp(data_84_species[c("log.degree_in", "log.strength_in", "log.betweenness")],
              center = TRUE, scale = TRUE)
summary(PCA)
#Importance of components:
#  PC1    PC2     PC3
#Standard deviation     1.6188 0.5808 0.20508
#Proportion of Variance 0.8736 0.1124 0.01402
#Cumulative Proportion  0.8736 0.9860 1.00000

### plot PCA ----
# create a data fame with Principal Component scores ----
PC <- data.frame(PCA$x[,1], PCA$x[,2])
p <- ggplot(PC, aes(PCA$x[,1], PCA$x[,2])) +
  geom_point(aes(),
             shape = 21,
             color = "gray50",
             #size = data_84_species$log.degree_in,
             #size = data_84_species$log.strength_in,
             size = data_84_species$log.betweenness,
             stroke = 1) +
  labs(x = "Principal Component 1 (87%)") +
  labs(y = "Principal Component 2 (11%)") +
  theme_classic() +
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))
p

### alternative plot ----
biplot(PCA, scale = 0, xlabs=rep("●", nrow(data_84_species)))

### get list of 100 IBAs predicted by PC1 ----
data <- data_84_species %>%
  mutate(pc1 = PCA$x[,1])

# write output to file
write.csv(data, "data_pc1.csv", row.names = FALSE)

df <- list()
for(i in 1:nrow(species)) {
  speciesName <- str_replace(species$species_name[i], " ", "_")
  top_10 <- data %>% 
    filter(species_name == speciesName) %>%
    nrow * 0.1 # top 10% 
  df[[i]] <- data %>% 
    filter(species_name == speciesName) %>%
    select(species_name, iba_label, pc1) %>%
    arrange(pc1) %>%
    head(round(top_10))
  print(speciesName)
}
top_list <- bind_rows(df)

top_100 <- top_list %>%
  arrange(pc1) %>%
  select(iba_label) %>%
  distinct %>%
  head(100) 
id <- rownames(top_100)
top_100 <- top_100 %>%
  mutate(priority = as.numeric(id))

### get number of species for which each IBA on the global list was on the species list ----
data_top_100 <- inner_join(top_100,
                           inner_join(top_100, top_list, by = "iba_label") %>%
                             group_by(iba_label) %>% summarise(n = n()),
                           by = "iba_label")

### get longitude and latitude ----
nodes_top_100_all <- inner_join(data_top_100, ibas, by = "iba_label")

### save list to file ----
write_csv(nodes_top_100_all, "nodes_100_list_84_species.csv")






##########################

### compare flu list to this list ----
flu_100_nodes <- read_csv("nodes_100_list.csv")

shared_nodes <- inner_join(nodes_top_100_all %>% select(iba_label),
                           flu_100_nodes %>% select(iba_label),
                           by = "iba_label")
# 65 IBAs were selected when using data from the 22 species
# and also when using data from the 84 species, so
# then, we get 35 additional IBAs from the 84 species

### get list of 100 IBAs predicted by PC1 ----
data <- data_84_species %>%
  mutate(pc1 = PCA$x[,1])

df <- list()
for(i in 1:nrow(species)) {
  speciesName <- str_replace(species$species_name[i], " ", "_")
  top_10 <- data %>% 
    filter(species_name == speciesName) %>%
    nrow * 0.1 # top 10% 
  df[[i]] <- data %>% 
    filter(species_name == speciesName) %>%
    select(species_name, iba_label, pc1) %>%
    arrange(pc1) %>%
    head(round(top_10))
  print(speciesName)
}
top_list <- bind_rows(df)

top_165 <- top_list %>%
  arrange(pc1) %>%
  select(iba_label) %>%
  distinct %>%
  head(165) # 100 + additional 65 IBAs
id <- rownames(top_165)
top_165 <- top_165 %>%
  mutate(priority = as.numeric(id))

### get number of species for which each IBA on the global list was on the species list ----
data_top_165 <- inner_join(top_165,
                           inner_join(top_165, top_list, by = "iba_label") %>%
                             group_by(iba_label) %>% summarise(n = n()),
                           by = "iba_label")

### get longitude and latitude ----
nodes_top_165_all <- inner_join(data_top_165, ibas, by = "iba_label")

list_IBAs <- anti_join(nodes_top_165_all, shared_nodes, by = "iba_label")


### plot the nodes and save the map ----
svg(filename = "fig_4a.svg",
#svg(filename = "fig_4b.svg", 
    width=8, 
    height=6, 
    pointsize=12)

# plot background map
maps::map("world",
          border = NA,
          fill = T,
          col = "gray90",
          bg = "white",
          projection = "lagrange",
          orientation = c(90,0,0),
          xlim = c(-25, 48),
          ylim = c(35.5, 65),
          lforce = "l")

# plot nodes from the 84 species
points(mapproject(nodes_top_100_all$longitude,
                  nodes_top_100_all$latitude,
                  projection = "lagrange",
                  orientation = c(90, 0, 0)),
       pch= 21,
       col = "white",
       bg = "orange",
       lwd = 1,
       cex = log(nodes_top_100_all$area_ha + 1) * 0.1)
       #cex = (nodes_top_100_all$area_ha / max(nodes_top_100_all$area_ha)) * 5)
       #cex = log(nodes_top_100_all$n + 1))

# plot nodes from the 22 species
points(mapproject(flu_100_nodes$longitude,
                  flu_100_nodes$latitude,
                  projection = "lagrange",
                  orientation = c(90, 0, 0)),
       pch= 21,
       col = "white",
       bg = "red",
       lwd = 1,
       cex = log(flu_100_nodes$area_ha + 1) * 0.1)
       #cex = (flu_100_nodes$area_ha / max(flu_100_nodes$area_ha)) * 5)
       #cex = log(flu_100_nodes$n + 1))

dev.off()
