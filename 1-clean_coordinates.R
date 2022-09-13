#---------------------------#
# 1. Clean coordinates 
#---------------------------#

library(tidyverse)
library(CoordinateCleaner)
library(data.table)

# load data sets
dir.df <- list.files("../../../../Data/datasets", full.names = T, pattern = ".csv")

# load data sets
i = 6 # choose data set
col.long <- "decimalLongitude"
col.lat <- "decimalLatitude"

df <- fread(dir.df[[i]])

df <- df %>% select(species, year, col.long, col.lat) %>% na.omit()
flags <- clean_coordinates(
  df,
  lon = col.long, 
  lat = col.lat,
  tests = c("capitals", "centroids", "equal", "seas", "zeros","institutions"),
  capitals_rad = 1000,
  centroids_rad = 1000
)
df.clean <- df[flags$.summary,]

name.group <- gsub("../../../../Data/datasets/", "", dir.df[[i]])
name.group <- gsub(".csv", "", name.group)

df.clean <- df.clean %>% filter(species != "") %>% select(1:4)
names(df.clean)[3:4] <- c("x", "y")

fwrite(df.clean, paste0("../../../../Data/datasets/cleaned/", name.group, ".csv"))
