#################################################################################
### mapping avian flu geographic locations to the global network of 2702 IBAs ###
#################################################################################

### load required packages ----
library("geosphere")
library("tidyverse")

### print package version ----
packageVersion("geosphere") # 1.5.7
packageVersion("tidyverse") # 1.2.1

### read data ----
ibas <- read_csv("../working_data/ibas_84_species.csv")
flu <- read_csv("avian_flu_working_dataset.csv")
species <- read_csv("../working_data/species_84.csv")
global_network <- read_csv("../working_data/links_84_species.csv")

flu_records <- flu %>%
  group_by(Host, Latitude, Longitude) %>%
  summarise(records = n()) %>%
  ungroup %>%
  group_by(Host) %>%
  summarise(n_records = sum(records), n_locations = n()) %>%
  ungroup %>%
  arrange(desc(n_locations))

### mapping longitude and latitude to any of our IBAs ----
flu_locations <- flu %>% distinct

ibas_flu <- data.frame(matrix(ncol = 2, nrow = nrow(flu_locations)))
colnames(ibas_flu) <- c("id", "iba_label")
for(i in 1: nrow(flu_locations)){
  for(j in 1: nrow(ibas)){
    radius <- sqrt((ibas[[j,5]] * 10000)/pi) + 10000
    d <- distm(
      c(ibas[[j,3]], ibas[[j,4]]), # longitude, latitude
      c(flu_locations[[i,4]], flu_locations[[i,3]]), # longitude, latitude
      fun = distHaversine) # it might take 5 minutes or longer
    if(d <= radius){ibas_flu[[i,1]] = i
    ibas_flu[[i,2]] = ibas[[j,1]]}
  }
}

flu_data <- as_tibble(cbind(flu_locations, ibas_flu))
# write_csv(flu_data, "flu_data.csv")

flu_locations_mapped <- left_join(flu_records,
                                  flu_data %>%
                                    na.omit %>%
                                    group_by(Host, Latitude, Longitude) %>%
                                    summarise(locations = n()) %>%
                                    ungroup %>%
                                    group_by(Host) %>%
                                    summarise(n_locations = sum(locations)) %>%
                                    ungroup %>%
                                    arrange(desc(n_locations)),
                                  by = "Host") %>%
  rename(n_locations = n_locations.x, n_locations_mapped = n_locations.y)
# write_csv(flu_locations_mapped, "flu_locations_mapped.csv")

### summary ----
flu %>% nrow
# 12210 records of individuals infected (i.e., positive samples for the avian flu)

flu_locations %>% nrow
# 297 total locations where the avian flu was detected
# (i.e., infected individuals of different species where found in the same location)

flu %>% select(Latitude, Longitude) %>% distinct %>% nrow
# 185 distinct locations where the avian flu was detected

inner_join(flu %>% select(Latitude, Longitude) %>% distinct,
           flu_data %>% na.omit %>% select(Latitude, Longitude),
           by = c("Latitude", "Longitude")) %>%
  distinct %>% nrow
# 143 out of 185 (78%) flu locations were mapped into ...

flu_data %>% select(iba_label) %>% na.omit %>% distinct %>% nrow
# ... 55 distinct IBAs, and ...

flu_locations_mapped %>% select(Host) %>% distinct %>% nrow
# ... 38 infected species were been found in at least one of those 55 IBAs

flu_mapped_species_list <- inner_join(flu_data %>% na.omit %>% select(Host) %>% rename(species_name = Host),
                                      species %>% select(species_name),
                                      by = "species_name") %>% distinct
flu_mapped_species_list %>% nrow
# 35 out of 38 species ---for which at least one flu location was mapped--- were included in our 84-species list

# save the 35 species for which the avian flu was detected and that are included in our 84-species list
#write_csv(flu_mapped_species_list, "flu_mapped_species_list.csv")

### mapping the flu locations of the 35 species to our species networks ----
### to check whether the mapped flu locations correspond, at least, to one IBA of the species networks
# (i.e., a flu location for that species might not be included in the IBAs of the migration network of that species)

# replace species label by species name (only for the 22 species with flu) in the global network
global_web <- inner_join(global_network, 
                         inner_join(species, flu_mapped_species_list, by = "species_name") %>%
                           select(species_label, species_name),
                         by = "species_label") %>%
  mutate(species_name = str_replace(species_name, " ", "_"))

# match flu locations to the nodes of the 35 species-specific networks
flu_mapped_species_list_working <- list()
for(i in 1:nrow(flu_mapped_species_list)) {
  species <- str_replace(flu_mapped_species_list$species_name[i], " ", "_")
  flu_mapped_species_list_working[[i]] <- 
    inner_join(rbind(global_web %>%
                       filter(species_name == species) %>%
                       select(from) %>%
                       distinct %>%
                       rename(iba_label = from),
                     global_web %>%
                       filter(species_name == species) %>%
                       select(to) %>%
                       distinct %>%
                       rename(iba_label = to)
                    ) %>% distinct,
               flu_data %>% filter(Host == flu_mapped_species_list$species_name[i]) %>%
                 select(iba_label) %>%
                 na.omit %>%
                 distinct,
               by = "iba_label") %>%
      mutate(species_name = species)
}
flu_data_working <- bind_rows(flu_mapped_species_list_working)
flu_data_working
# 86 rows
# 22 out of 35 species have at least one flu location mapped into an IBA of its own species-specific network
# 35 distinct flu locations out of 55 correspond to IBAs of the networks of the 22 species
#write_csv(flu_data_working, "flu_data_working.csv")
#write_csv(flu_data_working %>% distinct(species_name), "flu_mapped_species_list_working.csv")
#write_csv(flu_data_working %>% distinct(iba_label), "flu_nodes.csv")
