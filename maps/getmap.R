### load required libraries ----
library("maps")
library("mapproj")
library("geosphere")
library("tidyverse")


### print package version ----
packageVersion("maps") # 3.3.0
packageVersion("mapproj") # 1.2.6
packageVersion("geosphere") # 1.5.7
packageVersion("tidyverse") # 1.2.1


### read data ----
iba_nodes <- read_csv("../working_data/ibas_84_species.csv")
global_network <- read_csv("../working_data/links_84_species.csv")
species <- read_csv("../working_data/species_84.csv")
df_flu <- read_csv("../glmer/df_flu.csv")

### replace species label by species name in the global network
global_web <- inner_join(global_network, species, by = "species_label") %>%
  select(-species_label, -individuals.y, -records) %>%
  rename(individuals = individuals.x)

### function to build species networks ----
get_speciesNetwork <- function(speciesName){
  speciesName <- speciesName
  df <- inner_join(inner_join(global_web %>%
                                filter(species_name == speciesName) %>%
                                rename(iba_label = from, iba_to = to),
                              iba_nodes,
                              by = "iba_label") %>%
                     rename(iba_from = iba_label, longitude_from = longitude, latitude_from = latitude),
                   inner_join(global_web %>%
                                filter(species_name == speciesName) %>%
                                rename(iba_label = to, iba_from = from),
                              iba_nodes,
                              by = "iba_label") %>%
                     rename(iba_to = iba_label, longitude_to = longitude, latitude_to = latitude),
                   by = c("iba_from", "iba_to", "species_name")) %>%
    rename(individuals = individuals.x, species = species_name, area_from = area_ha.x, area_to = area_ha.y) %>%
    select(-individuals.y, -iba_name.x, -iba_name.y)
  species_name <- str_replace(speciesName, " ", "_")
  write_csv(df, paste0("species_networks/", species_name, ".csv"))
}

# build and save the species networks ----
dir.create("species_networks")
for(i in 1:nrow(species)) {
  get_speciesNetwork(species$species_name[i])
  print(species$species_name[i])
}



##############
### Fig. 1 ###
##############
### save figure ----
svg(filename = "maps/fig_1.svg", 
    width=8, 
    height=6, 
    pointsize=12)
    # plot the map
    maps::map("world",
              border = NA,
              fill = T,
              col = "gray90",
              bg = "white",
              projection = "lagrange",
              orientation = c(90,0,0),
              xlim = range(iba_nodes$longitude) + 7,
              ylim = range(iba_nodes$latitude),
              lforce = "s"
    )
    # plot the links
    species <- species_84 %>% select(species_name) %>% distinct
    # get unique links from the networks of the 84 species and then choose randomly 10% and plot them
    df_links <- list()
    for(i in 1:nrow(species)) {
      speciesName <- str_replace(species$species_name[i], " ", "_")
      df_links[[i]] <- read_csv(paste0("species_networks/", speciesName, ".csv")) %>%
        select(iba_from, longitude_from, latitude_from, iba_to, longitude_to, latitude_to)
    }
    all_links <- bind_rows(df_links) %>% distinct
    n_links <- round(nrow(all_links) * 0.1)
    all_10_links <- all_links %>% sample_n(n_links)
    for (i in (1:nrow(all_10_links))) {
      links <- gcIntermediate(c(all_10_links$longitude_from[i], all_10_links$latitude_from[i]),
                              c(all_10_links$longitude_to[i], all_10_links$latitude_to[i]),
                              n = 1000)
      project_links <- mapproject(links[,1], links[,2], projection = "lagrange", orientation = c(90,0,0)) 
      lines(project_links,
            lty = 19,
            lwd = 0.15,
            col = "deepskyblue")
    }
    # plot the nodes
    points(mapproject(iba_nodes$longitude,
                      iba_nodes$latitude,
                      projection = "lagrange",
                      orientation = c(90, 0, 0)),
           pch= 21,
           col = "white",
           bg = "deepskyblue",
           lwd = 0.5,
           cex = iba_nodes$area_ha / max(iba_nodes$area_ha) + 0.5)
dev.off()


###############
### Fig. 2a ###
###############
### function to plot background map ----
plot_map_flu <- function(){
  maps::map("world",
            border = NA,
            fill = T,
            col = "gray90",
            bg = "white",
            projection = "lagrange",
            orientation = c(90,0,0),
            xlim = c(-12, 48),
            ylim = c(35.5, 60),
            lforce = "l"
  )
}

### function to plot flu locations ----
plot_flu_locations <- function() {
  locations_flu <- inner_join(iba_nodes %>% select(-iba_name),
                   df_flu %>% select(iba_label) %>% distinct,
                   by = "iba_label")
  points(mapproject(locations_flu$longitude,
                    locations_flu$latitude,
                    projection = "lagrange",
                    orientation = c(90, 0, 0)),
         pch= 21,
         col = "white",
         bg = "red",
         lwd = 1.5,
         cex = 2)
}

### function to plot both the incoming and outgoing links of the flu locations of a species network ----
get_flu_links <- function(speciesNetworkName, flu) {
  speciesNetworkName <- speciesNetworkName
  flu <- flu
  speciesNetwork <- read_csv(paste0("species_networks/", speciesNetworkName, ".csv"))
  rbind(inner_join(speciesNetwork %>%
                     select(iba_from, longitude_from, latitude_from, iba_to, longitude_to, latitude_to),
                   df_flu %>%
                     filter(species_name == speciesNetworkName & is_flu == flu) %>%
                     select(iba_label) %>%
                     rename(iba_to = iba_label),
                   by = "iba_to"),
        inner_join(speciesNetwork %>%
                     select(iba_from, longitude_from, latitude_from, iba_to, longitude_to, latitude_to),
                   df_flu %>%
                     filter(species_name == speciesNetworkName & is_flu == flu) %>%
                     select(iba_label) %>%
                     rename(iba_from = iba_label),
                   by = "iba_from")) %>% distinct
}

### save figure ----
svg(filename = "maps/fig_2a.svg", 
    width=8, 
    height=6, 
    pointsize=12)

    ### plot flu locations and their links ----
    # plot background map
    plot_map_flu()
    # plot links
    flu_species <- df_flu %>% select(species_name) %>% distinct
    df_list <- list()
    for(i in 1:nrow(flu_species)) {
      speciesName <- flu_species$species_name[i]
      df_list[[i]] <- get_flu_links(speciesNetworkName = speciesName, flu = "1")
      print(speciesName)
    }
    df_flu_links <- bind_rows(df_list) %>% distinct
    for (i in (1:nrow(df_flu_links))) {
      links <- gcIntermediate(c(df_flu_links$longitude_from[i], df_flu_links$latitude_from[i]),
                              c(df_flu_links$longitude_to[i], df_flu_links$latitude_to[i]),
                              n = 1000,
                              addStartEnd = TRUE)
      project_links <- mapproject(links[,1], links[,2], projection = "lagrange", orientation = c(90,0,0)) 
      lines(project_links,
            lty = 19,
            lwd = 0.15,
            col = "red")
    }
    # plot nodes
    plot_flu_locations()
dev.off()


#################
### Fig. 2b-d ###
#################
### function to plot background map focused on the avian flu locations
plot_map_zoom <- function(){
  maps::map("world",
            border = NA,
            fill = T,
            #col = "gray90",
            col = "white",
            bg = "white",
            projection = "lagrange",
            orientation = c(90,0,0),
            xlim = c(5, 6),
            ylim = c(51, 53.5),
            lforce = "s")
}

### function to plot flu or non-flu locations for a specific species ----
plot_species_locations <- function(speciesName, isFlu, color) {
  speciesName <- speciesName
  isFlu <- isFlu
  locations_flu <- inner_join(iba_nodes %>% select(-iba_name), df_flu, by = "iba_label") %>%
    filter(species_name == speciesName & is_flu == isFlu)
  points(mapproject(locations_flu$longitude,
                    locations_flu$latitude,
                    projection = "lagrange",
                    orientation = c(90, 0, 0)),
         pch= 21,
         col = "white",
         bg = color,
         lwd = 1.5,
         cex = locations_flu$pc1 * (-1)) # revert the direction of pc1 
}

### species ----
speciesName <- "Larus_ridibundus"
speciesName <- "Anas_platyrhynchos"
speciesName <- "Anser_albifrons"

### save figure ----
svg(filename = "maps/fig_2_Anser.svg", 
    width=8, 
    height=6, 
    pointsize=12)

  ### plot flu locations and their links for a specific species ----
  # plot background map
  plot_map_zoom()
  # plot links from and to non-flu locations
  df_links <- get_flu_links(speciesNetworkName = speciesName, flu = "0")
  for (i in (1:nrow(df_links))) {
    links <- gcIntermediate(c(df_links$longitude_from[i], df_links$latitude_from[i]),
                            c(df_links$longitude_to[i], df_links$latitude_to[i]),
                            n = 1000,
                            addStartEnd = TRUE)
    project_links <- mapproject(links[,1], links[,2], projection = "lagrange", orientation = c(90,0,0)) 
    lines(project_links,
          lty = 19,
          lwd = 0.25,
          col = "forestgreen")
  }
  # plot links from and to flu locations
  df_links <- get_flu_links(speciesNetworkName = speciesName, flu = "1")
  for (i in (1:nrow(df_links))) {
    links <- gcIntermediate(c(df_links$longitude_from[i], df_links$latitude_from[i]),
                            c(df_links$longitude_to[i], df_links$latitude_to[i]),
                            n = 1000,
                            addStartEnd = TRUE)
    project_links <- mapproject(links[,1], links[,2], projection = "lagrange", orientation = c(90,0,0)) 
    lines(project_links,
          lty = 19,
          lwd = 0.25,
          col = "red")
  }
  # plot flu locations
  plot_species_locations(speciesName = speciesName, isFlu = "1", color = "red")
  # plot non-flu locations
  plot_species_locations(speciesName = speciesName, isFlu = "0", color = "forestgreen")

dev.off()


##############
### Fig. 3 ###
##############
### function to plot the links of a species network ----
get_links <- function(speciesNetworkName) {
  speciesNetworkName <- speciesNetworkName
  speciesNetwork <- read_csv(paste0("species_networks/", speciesNetworkName, ".csv"))
  rbind(inner_join(speciesNetwork %>%
                     select(iba_from, longitude_from, latitude_from, iba_to, longitude_to, latitude_to),
                   df %>%
                     filter(species_name == speciesNetworkName) %>%
                     select(iba_label) %>%
                     rename(iba_to = iba_label),
                   by = "iba_to"),
        inner_join(speciesNetwork %>%
                     select(iba_from, longitude_from, latitude_from, iba_to, longitude_to, latitude_to),
                   df %>%
                     filter(species_name == speciesNetworkName) %>%
                     select(iba_label) %>%
                     rename(iba_from = iba_label),
                   by = "iba_from")) %>% distinct
}

### function to plot the nodes for a specific species ----
plot_nodes <- function(speciesName) {
  speciesName <- speciesName
  nodes <- inner_join(iba_nodes %>% select(-iba_name), df, by = "iba_label") %>%
    filter(species_name == speciesName)
  points(mapproject(nodes$longitude,
                    nodes$latitude,
                    projection = "lagrange",
                    orientation = c(90, 0, 0)),
         pch= 21,
         col = "white",
         bg = "deepskyblue",
         lwd = 1.5,
         cex = nodes$area_ha / max(iba_nodes$area_ha) + 1)
}

### species ----
speciesName <- "Anas_penelope"
speciesName <- "Anas_acuta"
# distance W_OS = 0.73

speciesName <- "Anas_penelope"
speciesName <- "Cygnus_columbianus"
# distance W_OS = 0.85

speciesName <- "Cygnus_columbianus"
speciesName <- "Tringa_totanus"
# distance W_OS = 0.94



### save figure ----
svg(filename = "maps/fig_3_Tringa_totanus.svg", 
    width=8, 
    height=6, 
    pointsize=12)

# plot the map
maps::map("world",
          border = NA,
          fill = T,
          col = "gray90",
          bg = "white",
          projection = "lagrange",
          orientation = c(90,0,0),
          xlim = range(iba_nodes$longitude) + 7,
          ylim = range(iba_nodes$latitude),
          lforce = "s"
)
# plot links
df_links <- get_links(speciesNetworkName = speciesName)
for (i in (1:nrow(df_links))) {
  links <- gcIntermediate(c(df_links$longitude_from[i], df_links$latitude_from[i]),
                          c(df_links$longitude_to[i], df_links$latitude_to[i]),
                          n = 1000,
                          addStartEnd = TRUE)
  project_links <- mapproject(links[,1], links[,2], projection = "lagrange", orientation = c(90,0,0)) 
  lines(project_links,
        lty = 19,
        lwd = 0.5,
        col = "deepskyblue")
}
# plot nodes
plot_nodes(speciesName = speciesName)

dev.off()