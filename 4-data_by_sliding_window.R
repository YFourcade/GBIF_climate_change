#-----------------------------------#
# 4. Extract data by sliding window 
#-----------------------------------#

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(data.table)
library(rnaturalearth)
sf::sf_use_s2(FALSE)

setwd("C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle")

### Data ###
# load data sets
dir.df <- list.files("./Data/datasets/cleaned/thinned", full.names = T, pattern = ".csv")
# load STI
dir.sti <- list.files("./Data/STI/", full.names = T, pattern = ".csv")

# Human influence index

# hii <- rast("./Data/hii_v2geo/hdr.adf")
# gdalUtils::gdalwarp(srcfile = "./Data/hii_v2geo/hdr.adf", 
#          dstfile = "./Data/hii_v2geo/hii_proj.tiff",
#          t_srs = crs("ESRI:54009"),
#          overwrite = T)

hii <- rast("./Data/hii_v2geo/hii_proj.tiff")

# load temperature
tmp <- rast("C:/Users/200597/Desktop/cru_ts4.04.1901.2019.tmp.dat.nc")
tmp9019 <-subset(tmp, 1068:1428)
yrs <- unique(year(time(tmp9019)))[-1]

# annual temperature
annual_means <- lapply(yrs,function(x) mean(tmp9019[[grep(x,time(tmp9019))]],na.rm=TRUE))
tmp9019_mean <-do.call(c, annual_means)
names(tmp9019_mean) <- yrs
tmp9019_mean_5years <- tapp(tmp9019_mean, rep(1:6,each=5), fun = mean)
names(tmp9019_mean_5years) <- c(1992,1997,2002,2007,2012,2017)

# spring/summer temperature
summer_temp <- tmp9019[[grep("-03-|-04-|-05-|-06-|-07-|-08-|-09-",time(tmp9019))]]
summer_means <- lapply(yrs,function(x) mean(summer_temp[[grep(x,time(summer_temp))]],na.rm=TRUE))
tmp9019_summer_mean <-do.call(c, summer_means)
names(tmp9019_summer_mean) <- yrs
tmp9019_summer_mean_5years <- tapp(tmp9019_summer_mean, rep(1:6,each=5), fun = mean)
names(tmp9019_summer_mean_5years) <- c(1992,1997,2002,2007,2012,2017)

# winter temperature
winter_temp <- tmp9019[[grep("-12-|-01-",time(tmp9019))]]
winter_temp <- winter_temp[[-nlyr(winter_temp)]]
tmp9019_winter_mean <- tapp(winter_temp, rep(1:30,each=2), fun = mean)
tmp9019_winter_mean_5years <- tapp(tmp9019_winter_mean, rep(1:6,each=5), fun = mean)
names(tmp9019_winter_mean_5years) <- c(1992,1997,2002,2007,2012,2017)

### Sliding windows ###
# création d'un raster de résolution 1*1
r2 <-rast(nrows=180, ncols=360, xmin=-180, xmax=180, ymin=-90, ymax=90,
          crs = "epsg:4326", resolution = 1)

# Conversion du raster en polygone
pol2 <- as.polygons(r2)
pol2 <- st_as_sf(pol2)
pol2$id_polygons <- 1:nrow(pol2)
pol2$id_polygons <- paste0("ID_", 1:nrow(pol2))

# centroïdes
cent2 <- st_centroid(pol2)
cent2$id_polygons <- pol2$id_polygons

# incorporate continent and coordinates
world <- st_as_sf(countries110) %>% st_buffer(., 1)
world <- world %>% group_by(continent) %>% summarise()
world <- world %>% filter(continent %in% c("Europe", "North America"))

cent2_cont <- extract(vect(world["continent"]), vect(cent2)) %>% 
  group_by(id.y) %>% summarise_all(first) %>% pull(continent)

cent_unproj <- bind_cols(id_polygons = pol2$id_polygons, 
                         continent = cent2_cont, 
                         X = st_coordinates(cent2)[,1],
                         Y = st_coordinates(cent2)[,2])
cent_unproj <- na.omit(cent_unproj)

# projection
cent_proj <- st_transform(cent2[cent2$id_polygons %in% cent_unproj$id_polygons,], "ESRI:54009")

# création de buffers autour de chaque centroïde
Buffer <- st_buffer(cent_proj, dist = 200000) %>% 
  mutate(X = st_coordinates(cent_proj)[,1], Y = st_coordinates(cent_proj)[,2]) 

### temperature per year ###
tmp9019_mean_5years.proj <- project(tmp9019_mean_5years, "ESRI:54009")
tmp9019_summer_mean_5years.proj <- project(tmp9019_summer_mean_5years, "ESRI:54009")
tmp9019_winter_mean_5years.proj <- project(tmp9019_winter_mean_5years, "ESRI:54009")

# annual
tmp_buf.annual <- exact_extract(tmp9019_mean_5years.proj, Buffer, fun = "mean")
tmp_buf.annual$id_polygons <- cent_proj$id_polygons
names(tmp_buf.annual)[1:6] <- c(1992, 1997, 2002, 2007, 2012, 2017)
tmp_buf.annual <- na.omit(tmp_buf.annual)
tmp_buf.annual <- tmp_buf.annual %>% pivot_longer(cols = 1:6, names_to = "Year", values_to = "Temperature")

# summer
tmp_buf.summer <- exact_extract(tmp9019_summer_mean_5years.proj, Buffer, fun = "mean")
tmp_buf.summer$id_polygons <- cent_proj$id_polygons
names(tmp_buf.summer)[1:6] <- c(1992, 1997, 2002, 2007, 2012, 2017)
tmp_buf.summer <- na.omit(tmp_buf.summer)
tmp_buf.summer <- tmp_buf.summer %>% pivot_longer(cols = 1:6, names_to = "Year", values_to = "Temperature")

# winter
tmp_buf.winter <- exact_extract(tmp9019_winter_mean_5years.proj, Buffer, fun = "mean")
tmp_buf.winter$id_polygons <- cent_proj$id_polygons
names(tmp_buf.winter)[1:6] <- c(1992, 1997, 2002, 2007, 2012, 2017)
tmp_buf.winter <- na.omit(tmp_buf.winter)
tmp_buf.winter <- tmp_buf.winter %>% pivot_longer(cols = 1:6, names_to = "Year", values_to = "Temperature")

# trend per window
tmp_trend_annual <- tmp_buf.annual %>% group_by(id_polygons) %>% 
  do(broom::tidy(lm(Temperature ~ as.numeric(Year), .))) %>% 
  filter(term == "as.numeric(Year)") %>% select(1,3) %>% rename(trend = estimate) %>% 
  left_join(
    tmp_buf.annual %>% group_by(id_polygons) %>% 
      summarise(mean = mean(Temperature))
  )

tmp_trend_summer <- tmp_buf.summer %>% group_by(id_polygons) %>% 
  do(broom::tidy(lm(Temperature ~ as.numeric(Year), .))) %>% 
  filter(term == "as.numeric(Year)") %>% select(1,3) %>% rename(trend = estimate) %>% 
  left_join(
    tmp_buf.summer %>% group_by(id_polygons) %>% 
      summarise(mean = mean(Temperature))
  )

tmp_trend_winter <- tmp_buf.winter %>% group_by(id_polygons) %>% 
  do(broom::tidy(lm(Temperature ~ as.numeric(Year), .))) %>% 
  filter(term == "as.numeric(Year)") %>% select(1,3) %>% rename(trend = estimate) %>% 
  left_join(
    tmp_buf.winter %>% group_by(id_polygons) %>% 
      summarise(mean = mean(Temperature))
  )

tmp_trend <- bind_rows(
  Annual = tmp_trend_annual,
  Summer = tmp_trend_summer,
  Winter = tmp_trend_winter,
  .id = "Season"
) %>% ungroup

tmp_trend <- tmp_trend %>% left_join(cent_unproj)

write_csv(tmp_trend, "./Results/temp_windows.csv")

### HII ###
hii.buf <- exact_extract(hii, Buffer, fun = "mean")
hii.buf <- cbind.data.frame(id_polygons = cent_proj$id_polygons, hii = hii.buf)

write_csv(hii.buf, "./Results/hii_windows.csv")


### CTI ###
for(i in 1:length(dir.df)){
  name.group <- gsub("./Data/datasets/cleaned/thinned/", "", dir.df[[i]])
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  # filter wrong dates and remove Hawaii + French Guiana
  df.sf.tmp <- read_csv(dir.df[[i]]) %>% filter(year %in% c(1992,1997,2002,2007,2012,2017))
  df.sf.tmp <- df.sf.tmp %>% filter(!(X < -150 & Y < 40), !(X > -63 & Y < 14))
  
  # remove rats, mice and domestic bees
  if(i == 2){
    df.sf.tmp <- df.sf.tmp %>% filter(!species %in% "Apis mellifera")
  } else if (i == 10){
    df.sf.tmp <- df.sf.tmp %>% filter(!species %in% c("Rattus norvegicus", "Mus musculus"))
  }

  df.sf <- st_as_sf(df.sf.tmp, coords = c("X", "Y"), crs ="epsg:4326")
  df.sf_proj <- st_transform(df.sf, crs = "ESRI:54009")
  
  # intersection occurrences et polygones
  intersect.tmp <-st_intersects(df.sf_proj, Buffer)
  intersect.tmp <-as.data.frame(intersect.tmp) %>% as.data.table
  df.sf.dt <- df.sf %>% st_drop_geometry %>% as.data.table()
  Buffer.dt <- Buffer %>% st_drop_geometry %>% as.data.table()
  
  df.Buffer <- cbind(
    df.sf.dt[intersect.tmp$row.id,],
    Buffer.dt[intersect.tmp$col.id,]
  )
  
  ## describe data
  # nb species and occurrences per windows per year
  Buffer_yr_occ_sp <- df.Buffer %>% 
    dplyr::group_by(id_polygons, year) %>%
    dplyr::summarise(n.sp.yr = length(unique(species)), n.occ = n()) 
  
  write_csv(Buffer_yr_occ_sp, paste0("./Results/Descr_per_buffer_", name.group ,".csv"))
  
  # nb species/window/occurrences per year
  yr_occ_sp <- df.Buffer %>% 
    dplyr::group_by(year) %>%
    dplyr::summarise(n.sp = length(unique(species)),
                     n.window = length(unique(id_polygons)),
                     n.occ = n())
  
  write_csv(yr_occ_sp, paste0("./Results/Descr_per_year_", name.group ,".csv"))
  
  # associer les STI
  sti <- fread(dir.sti[[i]])
  
  setkey(sti, "species")
  setkey(df.Buffer, "species")
  
  df.Buffer <-df.Buffer[sti]
  
  # Calcul du CTI
  CTI_anu <- df.Buffer %>%
    group_by(id_polygons, year) %>%
    summarise(cti = mean(sti, na.rm = TRUE))
  
  CTI_anu <- CTI_anu %>% left_join(
    df.Buffer %>%
      group_by(id_polygons, year) %>%
      summarise(n.sp = length(unique(species))) 
  )
  
  CTI_anu <- CTI_anu %>%
    filter(!n.sp == 1) %>% na.omit()
  
  write_csv(CTI_anu, paste0("./Results/CTI_annual_", name.group ,".csv"))
  
  # CTI trend
  CTI_trend <- CTI_anu %>% 
    group_by(id_polygons) %>% 
    do(broom::tidy(lm(cti ~ year, .))) %>% 
    filter(term == "year") %>% select(1,3,4) %>% rename(trend = estimate, trend.se = std.error) %>% 
    left_join(
      CTI_anu %>% group_by(id_polygons) %>% 
        summarise(mean = mean(cti))
    )
  
  write_csv(CTI_trend, paste0("./Results/CTI_trend_", name.group ,".csv"))
}
