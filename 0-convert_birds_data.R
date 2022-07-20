library(inborutils)
library(tidyverse)
library(sf)
sf::sf_use_s2(FALSE)
library(terra)
library(rnaturalearth)
library(doFuture)
library(data.table)
library(CoordinateCleaner)

registerDoFuture()
plan(multisession, workers = 10)
options(future.globals.maxSize = +Inf)

setwd("C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle")

#-----------------------------#
# 0. Convert bird data into
# database because too big :)
#-----------------------------#

sqlite_file1 <- "./Data/datasets/Aves_summer.sqlite"
table_name1 <- "birds_summer"

csv_to_sqlite(csv_file = "C:/Users/200597/Desktop/0357713-210914110416597.csv", 
              sqlite_file1, table_name1, pre_process_size = 1000, 
              delim = "\t",
              chunk_size = 50000, show_progress_bar = FALSE)

sqlite_file2 <- "./Data/datasets/Aves_winter.sqlite"
table_name2 <- "birds_winter"

csv_to_sqlite(csv_file = "C:/Users/200597/Desktop/0357715-210914110416597.csv", 
              sqlite_file2, table_name2, pre_process_size = 1000, 
              delim = "\t",
              chunk_size = 50000, show_progress_bar = FALSE)

#-----------------------------#
# 1-2. Clean bird coordinates
# and spatial thinning
#-----------------------------#
world <- st_as_sf(countries110) %>% st_buffer(., 1)
world <- world %>% group_by(continent) %>% summarise()
world <- world %>% filter(continent %in% c("Europe", "North America"))

## summer communities ##

# load database
# data <- DBI::dbConnect(RSQLite::SQLite(), "./Data/datasets/Aves_summer.sqlite")
# dbplyr::src_dbi(data)
# df <- tbl(data, "birds_summer", sep = "\t")

# extract all unique species
# sp <- df %>% pull(species) %>% unique %>% na.omit
# write_csv(cbind.data.frame(species = as.vector(sp)), "sp_birds.summer.csv")
sp <- read_csv("sp_birds.summer.csv") %>% pull(species)
sp <- gsub("/", "-", sp)

# find remaining species
dir <- "./Data/datasets/cleaned/thinned/birds_summer/"
dir.create(dir)

sp.done <- list.files(dir); sp.done <- gsub(".csv", "", sp.done)
sp.done.index <- match(sp.done, sp)

db.path <- "C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle/Data/datasets/Aves_summer.sqlite"

# loop for 
foreach(k = setdiff(1:length(sp), sp.done.index)) %dopar% {
  # select species
  sp.to.select <- sp[[k]]
  
  data <- DBI::dbConnect(RSQLite::SQLite(), db.path)
  # dbplyr::src_dbi(data)
  df.temp <- tbl(data, "birds_summer", sep = "\t") %>% filter(species == sp.to.select)
  # RSQLite::dbClearResult(RSQLite::dbSendQuery(data, "PRAGMA busy_timeout=5000;"));
  
  col.long <- "decimalLongitude"
  col.lat <- "decimalLatitude"
  
  # clean coordinates
  df.temp <- df.temp %>% select(species, year, all_of(c(col.long, col.lat)))
  df.temp <- as_tibble(df.temp, .name_repair = "minimal") %>% na.omit
  flags <- clean_coordinates(
    df.temp,
    lon = col.long, 
    lat = col.lat,
    tests = c("capitals", "centroids", "equal","institutions"),
    capitals_rad = 1000,
    centroids_rad = 1000
  )
  df.clean <- df.temp %>% filter(flags$.summary)
  df.clean <- df.clean %>% filter(species != "") %>% select(1:4)
  names(df.clean)[3:4] <- c("x", "y")
  
  # group by 5-year temporal window
  df.clean$year_5yrs <-as.factor(df.clean$year)
  df.clean$year_5yrs <- fct_recode(df.clean$year_5yrs, 
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
  
  if(nrow(df.clean) > 1){
    
    pts_thin_sp <- foreach(
      i = unique(df.clean$year_5yrs),
      .combine = rbind.data.frame
    ) %do% { 
      
      df.temp <- df.clean[df.clean$year_5yrs == i, c("x", "y")] %>% 
        distinct %>% 
        as.data.frame
      
      if(nrow(df.temp) > 1){
        # projeter occurrences
        pts <- st_as_sf(df.temp, coords = c("x","y"), crs ="epsg:4326")
        
        # créer un raster avec l'étendue des points et une resolution de 5 km
        r <- rast((ext(pts)+1)*1.2)
        res(r) <- 0.04
        
        # obtenir le numéro de pixel de chaque occurrence
        cell <- cellFromXY(r, st_coordinates(pts))
        dup <- duplicated(cell)
        
        # selectionner points non dupliqués
        pts_thinned <- st_coordinates(pts)[!dup, ]
        
      } else {
        pts_thinned <- df.temp %>% rename(X = x, Y = y)
      }
      
      if(is.null(dim(pts_thinned))){
        pts_thinned <- bind_cols(species = unique(df.clean$species), year = i, X = pts_thinned[1], Y = pts_thinned[2])
      } else {
        pts_thinned <- bind_cols(species = unique(df.clean$species), year = i, X = pts_thinned[,1], Y = pts_thinned[,2])  }
      
      return(pts_thinned)
    }
    
    # garder les occurrences d'AMnord et Europe
    pts_thin_sp_sf <- st_as_sf(pts_thin_sp, coords = c("X", "Y"), crs = 4326)
    
    pts_thin_sp <- pts_thin_sp[extract(vect(world["continent"]), vect(pts_thin_sp_sf)) %>% 
                                 group_by(id.y) %>% summarise_all(first) %>% pull(continent) %in% 
                                 c("Europe", "North America"),]
    
    # write results
    fwrite(pts_thin_sp, paste0(dir, gsub("/", "-", sp[[k]]), ".csv"))
  } else {
    fwrite(
      cbind.data.frame(
        species = sp[[k]],
        year = NA,
        X = NA,
        Y = NA
      ), paste0(dir, gsub("/", "-", sp[[k]]), ".csv")
    )
  }
}

birds_summer_ls <- lapply(list.files(dir, full.names = T), fread)
birds_summer_ls <- rbindlist(birds_summer_ls)
fwrite(birds_summer_ls, "./Data/datasets/cleaned/thinned/birds_summer.csv")

## winter communities ##
# load database
data <- DBI::dbConnect(RSQLite::SQLite(), "./Data/datasets/Aves_winter.sqlite")
dbplyr::src_dbi(data)
df <- tbl(data, "birds_winter", sep = "\t")

# extract all unique species
# sp <- df %>% pull(species) %>% unique %>% na.omit
# write_csv(cbind.data.frame(species = as.vector(sp)), "sp_birds.winter.csv")
sp <- read_csv("sp_birds.winter.csv") %>% pull(species)
sp <- gsub("/", "-", sp)

# find remaining species
dir <- "./Data/datasets/cleaned/thinned/birds_winter/"
dir.create(dir)

sp.done <- list.files(dir); sp.done <- gsub(".csv", "", sp.done)
sp.done.index <- match(sp.done, sp)

db.path <- "C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle/Data/datasets/Aves_winter.sqlite"

# loop
foreach(k = setdiff(1:length(sp), sp.done.index)) %dopar% {
  # select species
  sp.to.select <- sp[[k]]
  
  data <- DBI::dbConnect(RSQLite::SQLite(), db.path)
  # dbplyr::src_dbi(data)
  df.temp <- tbl(data, "birds_winter", sep = "\t") %>% filter(species == sp.to.select)
  
  col.long <- "decimalLongitude"
  col.lat <- "decimalLatitude"
  
  # clean coordinates
  df.temp <- df.temp %>% select(species, year, all_of(c(col.long, col.lat)), month)
  df.temp <- as_tibble(df.temp) %>% na.omit
  flags <- clean_coordinates(
    df.temp,
    lon = col.long, 
    lat = col.lat,
    tests = c("capitals", "centroids", "equal","institutions"),
    capitals_rad = 1000,
    centroids_rad = 1000
  )
  df.clean <- df.temp %>% filter(flags$.summary)
  df.clean <- df.clean %>% filter(species != "") %>% select(1:5)
  names(df.clean)[3:4] <- c("x", "y")
  
  # group by 5-year temporal window
  df.clean$year_5yrs <- df.clean$year
  df.clean <- df.clean %>% mutate(
    year_5yrs = ifelse(
      month == 12, year_5yrs + 1, year_5yrs
    )
  )
  df.clean$year_5yrs <- as.factor(df.clean$year_5yrs)
  df.clean$year_5yrs <- fct_recode(df.clean$year_5yrs, 
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
  
  if(nrow(df.clean) > 1){
    
    pts_thin_sp <- foreach(
      i = unique(df.clean$year_5yrs),
      .combine = rbind.data.frame
    ) %do% { 
      
      df.temp <- df.clean[df.clean$year_5yrs == i, c("x", "y")] %>% 
        distinct %>% 
        as.data.frame
      
      if(nrow(df.temp) > 1){
        # projeter occurrences
        pts <- st_as_sf(df.temp, coords = c("x","y"), crs ="epsg:4326")
        
        # créer un raster avec l'étendue des points et une resolution de 5 km
        r <- rast((ext(pts)+1)*1.2)
        res(r) <- 0.04
        
        # obtenir le numéro de pixel de chaque occurrence
        cell <- cellFromXY(r, st_coordinates(pts))
        dup <- duplicated(cell)
        
        # selectionner points non dupliqués
        pts_thinned <- st_coordinates(pts)[!dup, ]
        
      } else {
        pts_thinned <- df.temp %>% rename(X = x, Y = y)
      }
      
      if(is.null(dim(pts_thinned))){
        pts_thinned <- bind_cols(species = unique(df.clean$species), year = i, X = pts_thinned[1], Y = pts_thinned[2])
      } else {
        pts_thinned <- bind_cols(species = unique(df.clean$species), year = i, X = pts_thinned[,1], Y = pts_thinned[,2])  }
      
      return(pts_thinned)
    }
    
    # garder les occurrences d'AMnord et Europe
    pts_thin_sp_sf <- st_as_sf(pts_thin_sp, coords = c("X", "Y"), crs = 4326)
    
    pts_thin_sp <- pts_thin_sp[extract(vect(world["continent"]), vect(pts_thin_sp_sf)) %>% 
                                 group_by(id.y) %>% summarise_all(first) %>% pull(continent) %in% 
                                 c("Europe", "North America"),]
    
    # write results
    fwrite(pts_thin_sp, paste0(dir, gsub("/", "-", sp[[k]]), ".csv"))
  } else {
    fwrite(
      cbind.data.frame(
        species = sp[[k]],
        year = NA,
        X = NA,
        Y = NA
      ), paste0(dir, gsub("/", "-", sp[[k]]), ".csv")
    )
  }
}

birds_winter_ls <- lapply(list.files(dir, full.names = T), fread)
birds_winter_ls <- rbindlist(birds_winter_ls)
fwrite(birds_winter_ls, "./Data/datasets/cleaned/thinned/birds_winter.csv")

