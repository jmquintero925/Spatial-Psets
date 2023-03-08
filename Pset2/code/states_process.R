### Turn degree into csv ###

library(tidyverse)
library(readxl)
library(sf) 

`%notin%` <- Negate(`%in%`)

states_shp <- st_read("data/degree_states", "degree_states")
states_shp <- cbind(states_shp, st_coordinates(states_shp))
states_shp <- states_shp %>% 
    select(PRIMARY_ST, X, Y, CELL_ID) %>%
    mutate(X = as.character(round(X, 1)), 
           Y = as.character(round(Y, 1)))

gecon <- read.csv("data/upload_us_xi_090110.csv") %>%
    rename(X = longitude, Y=lat, CELL_ID = cellid) %>%
    mutate(X = as.character(X + 0.5),
           Y = as.character(Y + 0.5))

states_shp <- left_join(gecon, states_shp, by=c("X", "Y")) %>%
    select(X, Y, PRIMARY_ST) %>%
    filter(PRIMARY_ST %in% state.name) 

# Merge in policy index 
state_crosswalk <- cbind(state_abb = state.abb, 
        PRIMARY_ST = state.name) %>%
    as.data.frame()
policy <- read_excel("data/state_policies.xlsx") %>%
    rename(state_abb = `State/ Territory`) %>%
    left_join(state_crosswalk, by="state_abb") %>%
    mutate(PRIMARY_ST = ifelse(state_abb == "DC", 
            "Virginia", PRIMARY_ST)) %>%
    group_by(PRIMARY_ST) %>%
    summarize(num_policy = n()) %>%
    filter(!is.na(PRIMARY_ST))
states_shp <- left_join(states_shp, policy, by="PRIMARY_ST") %>%
    mutate(num_policy = ifelse(is.na(num_policy), 0, num_policy),
           num_policy = num_policy / max(num_policy))
write.csv(states_shp, "data/states.csv")
