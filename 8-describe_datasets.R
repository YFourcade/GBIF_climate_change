#-----------------------------------#
# 8. Describe datasets 
#-----------------------------------#

library(tidyverse)
library(data.table)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(viridis)
sf::sf_use_s2(FALSE)

# load data directories
dir.df.buff <- list.files("../../../../Results/", full.names = T, pattern = "Descr_per_buffer_")
dir.df.yr <- list.files("../../../../Results/", full.names = T, pattern = "Descr_per_year_")
temp_windows <- read_csv("../../../../Results/temp_windows.csv")

# find taxon names
name.group <- gsub("../../../../Results/Descr_per_year_", "", dir.df.yr)
name.group <- gsub(".csv", "", name.group)
print(name.group)

# read and merge with taxon names
df.yr <- lapply(
  1:length(dir.df.yr), 
  function(x){
    cbind.data.frame(Taxon = name.group[[x]], read_csv(dir.df.yr[[x]]))
  }
) %>% 
  rbindlist %>% 
  mutate(Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
         Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
         Taxon = gsub("Fourmis", "Formicidae", Taxon),
         Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
         Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
         Taxon = gsub("Rongeurs", "Rodentia", Taxon),
         Taxon = gsub("Caudata", "Urodela", Taxon)
  )
df.yr$Taxon <- as.factor(df.yr$Taxon)

df.buff <- lapply(
  1:length(dir.df.buff), 
  function(x){
    cbind.data.frame(Taxon = name.group[[x]], read_csv(dir.df.buff[[x]]))
  }
) %>% 
  rbindlist %>% 
  mutate(Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
         Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
         Taxon = gsub("Fourmis", "Formicidae", Taxon),
         Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
         Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
         Taxon = gsub("Rongeurs", "Rodentia", Taxon),
         Taxon = gsub("Caudata", "Urodela", Taxon)
  )
df.buff$Taxon <- as.factor(df.buff$Taxon)

# plots Fig 2a-b
p1 <- df.yr %>% mutate(name_lab = ifelse(year == 2017, as.character(Taxon), NA_character_)) %>% 
  ggplot(aes(x = year, y = n.occ, col = Taxon, group = Taxon)) +
  geom_point() +
  geom_path() +
  scale_x_continuous("Year", breaks = unique(df.yr$year),
                     expand = c(0, 0),
                     limits = c(1990, 2029)) + 
  scale_y_log10(
    "No. occurrences",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)) + 
  geom_text_repel(
    aes(color = Taxon, label = name_lab),
    direction = "y",
    xlim = c(2020, NA),
    hjust = 0,
    segment.linetype = "dotted",
    size = 3,
    box.padding = .3,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  theme_classic() +
  theme(legend.position = "none")


p2 <- df.yr %>% mutate(name_lab = ifelse(year == 2017, as.character(Taxon), NA_character_)) %>% 
  ggplot(aes(x = year, y = n.window, col = Taxon, group = Taxon)) +
  geom_point() +
  geom_path() +
  scale_x_continuous("Year", breaks = unique(df.yr$year),
                     expand = c(0, 0),
                     limits = c(1990, 2029)) + 
  scale_y_continuous("No. sliding windows") +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)) + 
  geom_text_repel(
    aes(color = Taxon, label = name_lab),
    direction = "y",
    xlim = c(2020, NA),
    hjust = 0,
    size = 3,
    segment.linetype = "dotted",
    box.padding = .3,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  theme_classic() +
  theme(legend.position = "none")

# plot windows ~ time
world <- rnaturalearth::ne_coastline(returnclass = "sf")

ggplot() +
  geom_raster(data = df.buff %>% left_join(temp_windows %>% filter(Season == "Annual")) %>% 
                rename(`No. of\noccurrences` = n.occ) %>% 
                filter(Taxon %in% levels(df.buff$Taxon)[1:5]),
              aes(x = X, y = Y, fill = `No. of\noccurrences`)) +
  geom_sf(data = world,
          color = "black", size = .01, fill = "white") +
  facet_grid(year ~ Taxon) +
  scale_fill_viridis(trans = "log", option = "inferno", breaks = c(1,100,10000)) +
  coord_sf(xlim = c(-170, 55), ylim = c(5, NA)) +
  theme_void() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.text.y = element_text(hjust = 0, margin = margin(1, 1, 1, 1, "lines")),
    strip.text.x = element_text(margin = margin(.5, .5, .5, .5, "lines")),
    legend.position = "bottom"
  )
ggsave("../../../../FigS1a.pdf", width = 10, height = 6)

ggplot() +
  geom_raster(data = df.buff %>% left_join(temp_windows %>% filter(Season == "Annual")) %>% 
                rename(`No. of\noccurrences` = n.occ) %>% 
                filter(Taxon %in% levels(df.buff$Taxon)[6:10]),
              aes(x = X, y = Y, fill = `No. of\noccurrences`)) +
  geom_sf(data = world,
          color = "black", size = .01, fill = "white") +
  facet_grid(year ~ Taxon) +
  scale_fill_viridis(trans = "log", option = "inferno", breaks = c(1,100,10000)) +
  coord_sf(xlim = c(-170, 55), ylim = c(5, NA)) +
  theme_void() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.text.y = element_text(hjust = 0, margin = margin(1, 1, 1, 1, "lines")),
    strip.text.x = element_text(margin = margin(.5, .5, .5, .5, "lines")),
    legend.position = "bottom"
  )
ggsave("../../../../FigS1b.pdf", width = 10, height = 6)


# data sources
dir.datasets <- list.files("../../../../Data/datasets/", full.names = T)[-7]

data_source <- NULL
for(i in 1:length(dir.datasets)){
  name.group <- gsub("../../../../Data/datasets/", "", dir.datasets)
  name.group <- gsub(".csv", "", name.group)
  name.group <- gsub(".sqlite", "", name.group)
  print(name.group[i])
  
  if(grepl(".csv", dir.datasets[i])){
    df <- fread(dir.datasets[i])[, c("species", "institutionCode")]
    df.sum <- df[,.N,by = institutionCode] %>% as.tibble()
  } else if(i == 3){
    df.db <- DBI::dbConnect(RSQLite::SQLite(), dir.datasets[i])
    df <- tbl(df.db, "birds_summer", sep = "\t") %>% select(species, institutionCode)
    df.sum <- df %>% group_by(institutionCode) %>% count() %>% collect() %>% rename(N = n)
  } else if(i == 4){
    df.db <- DBI::dbConnect(RSQLite::SQLite(), dir.datasets[i])
    df <- tbl(df.db, "birds_winter", sep = "\t") %>% select(species, institutionCode)
    df.sum <- df %>% group_by(institutionCode) %>% count() %>% collect() %>% rename(N = n)
  }
  df.sum <- df.sum %>% mutate(institutionCode = gsub("Non renseigné", "", institutionCode))
  df.sum <- df.sum %>% mutate(institutionCode = ifelse(is.na(institutionCode), "", institutionCode))
  df.sum <- df.sum %>% group_by(institutionCode) %>% summarise(N = sum(N)) %>% mutate(Taxon = name.group[i])
  data_source <- rbind.data.frame(data_source, df.sum)
}

write_csv(data_source, "../../../../Results/data_sources.csv")

# data_source <- read_csv("../../../../Results/data_sources.csv")

data_source <- data_source %>% 
  mutate(
    Taxon = gsub("Aves_summer", "Aves (summer)", Taxon),
    Taxon = gsub("Aves_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = gsub("Caudata", "Urodela", Taxon),
    Taxon = as.factor(Taxon)
  )

data_source <- data_source %>% group_by(Taxon) %>% 
  mutate(rang = rank(-N),
         n.taxon = sum(N)) %>% 
  ungroup %>% 
  mutate(proportion = N / n.taxon)


p3 <- data_source %>% filter(rang < 11) %>% 
  mutate(unknown = ifelse(is.na(institutionCode), "unknown", "known")) %>% 
  group_by(Taxon) %>% 
  mutate(name_lab = ifelse(rang == min(rang), as.character(Taxon), NA_character_)) %>% 
  ggplot(aes(x = rang, y = proportion, col = Taxon, group = Taxon)) +
  geom_point(aes(shape = unknown, size = unknown), alpha = .7) +
  geom_line(alpha = .7) +
  scale_x_continuous("Data provider rank",
                     expand = c(0, 0),
                     limits = c(.5, 11),
                     breaks = c(2,4,6,8,10)) + 
  scale_y_continuous("Proportion") +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)) + 
  scale_size_manual(values = c(1.5,3)) +
  geom_text_repel(
    aes(color = Taxon, label = name_lab),
    direction = "y",
    size = 3,
    xlim = c(4, NA),
    hjust = 0,
    segment.linetype = "dotted",
    box.padding = .3,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  theme_classic() +
  theme(legend.position = "none")

cowplot::plot_grid(p1,p2,p3, nrow = 3, labels = c("(a)","(b)","(c)"))
ggsave("../../../../Fig2.pdf", width = 4.5, height = 10)



### extract number of species for each taxon ###
dir.df <- list.files("../../../../Data/datasets/cleaned/thinned", full.names = T, pattern = ".csv")

sp <- c()
for(j in 1:length(dir.df)) {
  name.group <- gsub("../../../../Data/datasets/cleaned/thinned/", "", dir.df[[j]])
  name.group <- gsub("../../../.", "", name.group)
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  df <- fread(dir.df[[j]])
  sp.unique <- unique(df$species)
  
  sp <- rbind.data.frame(sp, cbind.data.frame(Taxon = name.group, sp.unique))
}

sp %>% group_by(Taxon) %>% summarise(n = n())
