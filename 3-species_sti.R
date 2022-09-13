#--------------------------#
# 3. Compute species' STI 
#--------------------------#

library(tidyverse)
library(terra)
library(data.table)

## leod data ##
# load data sets
dir.df <- list.files("../../../../Data/datasets/cleaned/thinned", full.names = T, pattern = ".csv")

# climate data CRU.TS4.04
tmp <- rast("cru_ts4.04.1901.2019.tmp.dat.nc")

# Compute mean annual/winter/summer temperature from monthly temperature data
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


## Calculate species' STI ##

dir.create("./Data/STI/")

for(i in 1:length(dir.df)){
  if(i %in% c(3,8)){
    temperature_rasters = tmp9019_summer_mean_5years
  } else if(i == 4) {
    temperature_rasters = tmp9019_winter_mean_5years
  } else {
    temperature_rasters = tmp9019_mean_5years
  }
  
  name.group <- gsub("../../../../Data/datasets/cleaned/thinned/", "", dir.df[[i]])
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  df <- fread(dir.df[[i]]) %>% as_tibble()
  anu_period <-split(df, df$year)
  
  # extract temperature at each occurrence
  temperature <-list()
  for (l in 1:nlyr(temperature_rasters)){
    temperature[[l]] <- extract(subset(temperature_rasters, l), anu_period[[l]][,3:4])[,2]
  }
  
  anu_sti <-lapply(X = 1:length(temperature), 
                   FUN = function(i)cbind(anu_period[[i]], 
                                          temperature[[i]]))
  
  anu_sti <-do.call(rbind, anu_sti)
  names(anu_sti)[names(anu_sti) == "temperature[[i]]"] <- "temperature"
  anu_sti <-anu_sti[!is.na(anu_sti$temperature),]
  
  # Mean temperature by species
  STI_anu <- anu_sti %>% 
    dplyr::group_by(species) %>%
    dplyr::summarise(sti = mean(temperature))
  
  fwrite(STI_anu, paste0("../../../../Data/STI/", name.group, ".csv"))
}
