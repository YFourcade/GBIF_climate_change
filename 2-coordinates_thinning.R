#---------------------------#
# 2. Coordinates thinning 
#---------------------------#

library(tidyverse)
library(data.table)
library(sf)
sf::sf_use_s2(FALSE)
library(terra)
library(rnaturalearth)
library(doFuture)

setwd("C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle")

# load data sets
dir.df <- list.files("./Data/datasets/cleaned", full.names = T, pattern = ".csv")

# sous-échantillonnage + sélectionner occurrences en AmN et Europe
world <- st_as_sf(countries110) %>% st_buffer(., 1)
world <- world %>% group_by(continent) %>% summarise()
world <- world %>% filter(continent %in% c("Europe", "North America"))

dir.create("./Data/datasets/cleaned/thinned") 

registerDoFuture()
plan(multisession, workers = 5)
options(future.globals.maxSize = +Inf)

for(j in c(3,7)){
  
  name.group <- gsub("./Data/datasets/cleaned/", "", dir.df[[j]])
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  data_year <- fread(dir.df[[j]])
  
  data_year$year_5yrs <-as.factor(data_year$year)
  data_year$year_5yrs <- fct_recode(data_year$year_5yrs, 
                                    "1992" = "1990",
                                    "1992" = "1991",
                                    "1992" = "1992",
                                    "1992" = "1993",
                                    "1992" = "1994",
                                    "1997" = "1995",
                                    "1997" = "1996",
                                    "1997" = "1997",
                                    "1997" = "1998",
                                    "1997" = "1999",
                                    "2002" = "2000",
                                    "2002" = "2001",
                                    "2002" = "2002",
                                    "2002" = "2003",
                                    "2002" = "2004",
                                    "2007" = "2005",
                                    "2007" = "2006",
                                    "2007" = "2007",
                                    "2007" = "2008",
                                    "2007" = "2009",
                                    "2012" = "2010",
                                    "2012" = "2011",
                                    "2012" = "2012",
                                    "2012" = "2013",
                                    "2012" = "2014",
                                    "2017" = "2015",
                                    "2017" = "2016",
                                    "2017" = "2017",
                                    "2017" = "2018",
                                    "2017" = "2019")

  # créer un raster avec l'étendue des points et une resolution de 5 km
  r <- rast((ext(world)+1)*1.2)
  res(r) <- 0.04
  crs(r) <- "epsg:4326"
  r <- wrap(r)
  
  data_year$breaks <- cut_number(1:nrow(data_year), 5, labels = F)
  
  cells.pts <- foreach(
    i = 1:5,
    .combine = c
  ) %dopar% { 
    cn <- extract(rast(r), as.matrix(data_year[data_year$breaks == i,3:4]), cells = T)
    cn$cell
  }
  
  data_year <- data_year[,cn:=cells.pts]
  pts_thin_sp <- data_year[,.SD[sample(.N, min(1,.N))],by = .(species, year_5yrs, cn)]
  pts_thin_sp <- pts_thin_sp[,c(1,2,5,6)]
  names(pts_thin_sp)[2:4] <- c("year", "X","Y")
  
  # garder les occurrences d'AMnord et Europe
  pts_thin_sp_sf <- st_as_sf(pts_thin_sp, coords = c("X", "Y"), crs = 4326)
  
  pts_thin_sp <- pts_thin_sp[extract(vect(world["continent"]), vect(pts_thin_sp_sf))[,2] %in% 
                               c("Europe", "North America"),]
  
  fwrite(pts_thin_sp, paste0("./Data/datasets/cleaned/thinned/", name.group, ".csv"))
  
}
