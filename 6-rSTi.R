#---------------------------#
# 6. compute species' rSTI 
#---------------------------#

library(tidyverse)
library(terra)
library(rnaturalearth)
library(data.table)
library(mgcv)
library(sf)
library(foreach)
sf::sf_use_s2(FALSE)

setwd("C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle")

### load data sets directories ###
dir.df <- list.files("./Data/datasets/cleaned/thinned", full.names = T, pattern = ".csv")
dir.sti <- list.files("./Data/STI/", full.names = T, pattern = ".csv")
dir.cti <- list.files("./Results/", pattern = "CTI_annual",full.names = T)

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
Buffer.dt <- Buffer %>% st_drop_geometry %>% as.data.table()
Buffer.dt$col.id <- 1:nrow(Buffer.dt)

# taxa to do
sp <- gsub(".csv", "", dir.df); sp <- gsub("./Data/datasets/cleaned/thinned/", "", sp)

# find remaining taxa
sp.done <- list.files("./Results/", pattern = "gained_"); sp.done <- gsub(".csv", "", sp.done); sp.done <- gsub("gained_", "", sp.done)
sp.done.index <- match(sp.done, sp)

foreach(j = setdiff(1:length(sp), sp.done.index)) %do% {
  
  name.group <- gsub("./Data/datasets/cleaned/thinned/", "", dir.df[[j]])
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  anu_cti_Ck <- vector('list', 6)
  names(anu_cti_Ck) <- c(1992,1997,2002,2007,2012,2017)
  for(z in c(1992,1997,2002,2007,2012,2017)){
    print(z)
    df.t <- fread(dir.df[[j]])[year == z,] %>% na.omit
    
    # remove Hawaii + French Guiana
    df.t <- df.t[!(X < -150 & Y < 40) & !(X > -63 & Y < 14),]
    
    # remove rats, mice and domestic bees
    if(j == 2){
      df.t <- df.t[!species %in% "Apis mellifera",]
    } else if (j == 10){
      df.t <- df.t[!species %in% c("Rattus norvegicus", "Mus musculus"),]
    }
    
    df.sf <- st_as_sf(df.t, coords = c("X", "Y"), crs ="epsg:4326")
    df.sf_proj <- st_transform(df.sf, crs = "ESRI:54009")
    
    ### merge occurrences, CTI and STI ###
    intersect.tmp <-st_intersects(df.sf_proj, Buffer)
    intersect.tmp <-as.data.frame(intersect.tmp) %>% as.data.table
    df.sf.dt <- df.sf %>% st_drop_geometry %>% as.data.table
    df.sf.dt$row.id <- 1:nrow(df.sf.dt)
    
    intersect.tmp[df.sf.dt, on = "row.id", `:=`(species = species, year = year)]
    intersect.tmp[Buffer.dt, on = "col.id", `:=`(id_polygons = id_polygons, X = X, Y = Y)]
    
    STI <- fread(dir.sti[[j]])
    CTI <- fread(dir.cti[[j]])
    
    intersect.tmp[STI, on = 'species', sti := i.sti]
    intersect.tmp[CTI, on = .(id_polygons, year), cti := i.cti]
    
    anu_cti_Ck[[as.character(z)]] <- intersect.tmp[,3:9]
    
    rm(intersect.tmp, df.sf.dt, df.sf_proj, df.sf); gc()
  }
  anu_cti_Ck <- rbindlist(anu_cti_Ck)
  
  ### find gained and lost species ###
  new_anu <-data.table()
  ext_anu <-data.table()
  # memory.limit(size=50000)
  for(i in unique(anu_cti_Ck$id_polygons)){
    d <- anu_cti_Ck[id_polygons == i,] %>% mutate(year = as.factor(year)) %>% droplevels
    yr <- sort(levels(d$year)) %>% as.factor()
    for(j in yr[!yr %in% yr[1]]){
      sp_t1 <- d[year == j,"species"] %>% unique
      t1 <-d[year == j,"species"]
      sp_t0 <- d[year == yr[which(yr == j)-1],"species"] %>% unique
      t0 <- d[year == yr[which(yr == j)-1],]
      new_sp2 <- as.data.table(setdiff(sp_t1, sp_t0))
      names(new_sp2) <-"species"
      if(nrow(new_sp2) >0)
        new_anu <- rbind(
          new_anu,
          cbind(
            new_sp2, id_polygons = i, 
            interval = paste(yr[which(yr == j)-1], "-", j),
            year = yr[which(yr == j)-1]
          )
        )
      ext_sp2 <- as.data.table(setdiff(sp_t0, sp_t1))
      names(ext_sp2) <-"species"
      if(nrow(ext_sp2) >0)
        ext_anu <- rbind(
          ext_anu, 
          cbind(
            ext_sp2, id_polygons = i, 
            interval = paste(yr[which(yr == j)-1], "-", j),
            year = yr[which(yr == j)-1]
          )
        )
      # print(i)
      # print(j)
      
      rm(new_sp2, ext_sp2); gc()
    }
    # rm(i)
  }
  
  # gained species
  new_spanu <- distinct(new_anu) %>% as_tibble %>% 
    mutate(year = as.numeric(as.character(year))) %>% 
    left_join(STI) %>% left_join(CTI) %>% select(-year)
  
  new_spanu$rSTI <-new_spanu$sti - new_spanu$cti
  
  write_csv(new_spanu, paste0("./Results/gained_", name.group ,".csv"))
  
  # lost species
  ext_spanu <-distinct(ext_anu) %>% as_tibble %>% 
    mutate(year = as.numeric(as.character(year))) %>% 
    left_join(STI) %>% left_join(CTI) %>% select(-year)
  
  ext_spanu$rSTI <-ext_spanu$sti - ext_spanu$cti
  
  write_csv(ext_spanu, paste0("./Results/lost_", name.group ,".csv"))
  
}


### compute average rSTI ###
sites <- read_csv("./Results/temp_windows.csv")

coefs.rsti <- foreach(j = 1:length(dir.df), .combine = rbind.data.frame) %do% {
  
  name.group <- gsub("./Data/datasets/cleaned/thinned/", "", dir.df[[j]])
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  new_spanu  <- read_csv(paste0("./Results/gained_", name.group ,".csv"))
  ext_spanu  <- read_csv(paste0("./Results/lost_", name.group ,".csv"))
  
  new_spanu <- new_spanu %>% mutate(year = as.numeric(substr(interval, 1,4))) %>% 
    left_join(read_csv(paste0(list.files("./Results/", pattern = "Descr_per_buffer",full.names = T)[[j]]))) %>% 
    left_join(sites %>% group_by(id_polygons) %>% 
                summarise_at(.vars = c("X","Y"), .funs = unique)) %>% 
    mutate_at(.vars = c("interval", "species", "id_polygons"), .funs = as.factor) %>% na.omit
  
  Mgain <- new_spanu %>% summarise(estimate = t.test(rSTI)$estimate,
                          lwr.ci = t.test(rSTI)$conf.int[1],
                          upr.ci = t.test(rSTI)$conf.int[2],
                          n.pos = length(rSTI[rSTI > 0]) / n(),
                          n.neg = length(rSTI[rSTI < 0]) / n())
  
  ext_spanu <- ext_spanu %>% mutate(year = as.numeric(substr(interval, 1,4))) %>% 
    left_join(read_csv(paste0(list.files("./Results/", pattern = "Descr_per_buffer",full.names = T)[[j]]))) %>% 
    left_join(sites %>% group_by(id_polygons) %>% 
                summarise_at(.vars = c("X","Y"), .funs = unique)) %>% 
    mutate_at(.vars = c("interval", "species", "id_polygons"), .funs = as.factor) %>% na.omit
  
  Mlost <- ext_spanu %>% summarise(estimate = t.test(rSTI)$estimate,
                          lwr.ci = t.test(rSTI)$conf.int[1],
                          upr.ci = t.test(rSTI)$conf.int[2],
                          n.pos = length(rSTI[rSTI > 0]) / n(),
                          n.neg = length(rSTI[rSTI < 0]) / n())
  
  # Mgain <-bam(rSTI ~ 1 + 
  #               s(id_polygons, bs = "re") + 
  #               s(species, bs = "re") + 
  #               s(X, Y, bs = "sos") + 
  #               s(interval, bs = "re"),
  #             weights = n.occ, 
  #             discrete = T, 
  #             nthreads = 2,
  #             data = new_spanu)
  # 
  # coefs.rsti.gain <- summary(Mgain, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
  # rm(Mgain, new_spanu); gc()
  # 
  # Mlost <-bam(rSTI ~ 1 + 
  #               s(id_polygons, bs = "re") + 
  #               s(species, bs = "re") + 
  #               s(X, Y, bs = "sos") + 
  #               s(interval, bs = "re"),
  #             weights = n.occ, 
  #             discrete = T, 
  #             nthreads = 2,
  #             data = ext_spanu)
  # 
  # coefs.rsti.lost <- summary(Mlost, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
  
  rr <- bind_cols(Taxon = name.group, bind_rows(
    Gain = Mgain, 
    Loss = Mlost,
    .id = "Type"))

  # write_csv(rr, paste0("./Results/rsti_coefs_", name.group ,".csv"))
  return(rr)
}
write_csv(coefs.rsti, paste0("./Results/rsti_coefs" ,".csv"))

# plot
coefs.rsti <- read_csv("./Results/rsti_coefs.csv")

coefs.rsti <- coefs.rsti %>% 
  mutate(
    Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
    Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = as.factor(Taxon)
  )

taxa.col.order <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)
names(taxa.col.order) <- unique(coefs.rsti$Taxon)

taxa.order <- cbind.data.frame(Taxon = c("Aves (winter)",
                                             "Aves (summer)",
                                             "Chiroptera",
                                             "Caudata",
                                             "Anura",
                                             "Lumbricidae",
                                             "Apidae",
                                             "Lepidoptera",
                                             "Formicidae",
                                             "Rodentia"))
coefs.rsti <- coefs.rsti %>% 
  mutate(Taxon = factor(Taxon, levels = c(as.character(rev(taxa.order$Taxon)))))

coefs.rsti %>% ggplot(aes(x = estimate, y = Taxon, xmin = lwr.ci, xmax = upr.ci, color = Type, shape = Type)) +
  geom_vline(xintercept = 0, col = 'red2', linetype = 2) + 
  geom_point(position = position_dodge(.5)) + 
  geom_errorbar(position = position_dodge(.5), width = 0) +
  scale_x_continuous("Relative STI (rSTI, mean +/- 95% CI)") +
  scale_y_discrete("") +
  scale_color_manual("", values = c("#f00d60", "#1e459c"), labels = c("Gained species", "Lost species")) + 
  scale_shape_discrete("", labels = c("Gained species", "Lost species")) + 
  theme_minimal() + 
  theme(legend.position = c(.8, .8),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.box.margin	= margin(0,0,0,0, unit = 'pt'),
        legend.title=element_blank(),
        panel.grid.minor.x = element_line(colour = NA),
        axis.text.y = element_text(colour = taxa.col.order[rev(taxa.order$Taxon)])
  )

ggsave(filename = "Fig4a.pdf", height = 4, width = 5)


coefs.rsti %>% 
  mutate(Type = ifelse(Type == "Gain", "Gained species", "Lost species")) %>% 
  pivot_longer(cols = c(6,7)) %>% 
  ggplot(aes(y = Taxon, x = value, fill = name)) +
  geom_bar(stat = 'identity') +
  geom_vline(xintercept = 0.5) +
  scale_y_discrete("") +
  scale_x_continuous("Proportion", labels = scales::percent_format(scale = 1)) +
  scale_fill_manual("", values = c("#A2B5CD", "#CD3333"), labels = c("rSTI < 0", "rSTI > 0")) + 
  facet_grid(~Type) +
  theme_minimal() +
  theme(
        axis.text.y = element_text(colour = taxa.col.order[rev(taxa.order$Taxon)]),
        panel.grid = element_line(colour = NA)
  )

ggsave(filename = "Fig4b.pdf", height = 4, width = 5.5)


coefs.rsti %>% pivot_wider(names_from = "Type", values_from = c(estimate,lwr.ci,upr.ci,n.pos,n.neg))
