#-------------------------------------#
# 5. Analyses at the community level 
#-------------------------------------#

library(tidyverse)
library(lmerTest)
library(doFuture)
library(doRNG)

# load data sets
dir.df <- "../../../../Results/"
temp_windows <- read_csv("../../../../Results/temp_windows.csv")
hii_windows <- read_csv("../../../../Results/hii_windows.csv")
eco.windows <- read_csv("../../../../Results/ecoregions.csv")

hii_windows %>% left_join(temp_windows) %>% 
  filter(Season == "Annual") %>% 
  ggplot(aes(x = X, y = Y, color = trend)) +
  geom_point()

hii_windows %>% left_join(temp_windows) %>% 
  filter(Season == "Annual") %>% 
  ggplot(aes(x = Y, y = trend, color = trend)) +
  geom_point()

hii_windows %>% left_join(temp_windows) %>% 
  filter(Season == "Annual") %>% 
  ggplot(aes(x = Y, y = hii, color = trend)) +
  geom_point()

# statistical analyses
registerDoFuture()
plan(multisession, workers = 5)

res <- foreach(i = 1:10, .combine = rbind.data.frame) %dorng% {
  name.group <- gsub("CTI_annual_", "", list.files(dir.df, pattern = "CTI_annual")[[i]])
  name.group <- gsub(".csv", "", name.group)
  print(name.group)
  
  # analysis CTI ~ Year
  if(i %in% c(3,8)){
    Season.df = "Summer"
  } else if(i == 4) {
    Season.df = "Winter"
  } else {
    Season.df = "Annual"
  }
  
  CTI_yr <- read_csv(list.files(dir.df, pattern = "CTI_annual",full.names = T)[[i]]) %>% 
    left_join(temp_windows %>% filter(Season == Season.df)) %>% 
    left_join(read_csv(paste0(list.files(dir.df, pattern = "Descr_per_buffer",full.names = T)[[i]]))) %>% 
    mutate(id_polygons = as.factor(id_polygons), year = year - 1992) %>% 
    left_join(eco.windows)

  m1 <- lmer(
    cti ~ year + mean + continent + 
      (year|id_polygons) + (1|ECO_NAME),
    weights = n.occ,
    data = CTI_yr
  )

  # resid.sampl <- sample_n(residuals(m1) %>% as.data.frame() %>% rownames_to_column(), 1000)
  # sp.corel <- ncf::spline.correlog(x = CTI_yr[resid.sampl$rowname,"X"]$X,
  #                                  y = CTI_yr[resid.sampl$rowname,"Y"]$Y, 
  #                                  z = resid.sampl[,2], 
  #                                  resamp = 100, latlon = T, xmax = 1000)
  # plot(sp.corel ,main = paste("Taxon :", name.group, "- Model : m1"))
  # summary(sp.corel)
  
  # library(DHARMa)
  # plot(simulateResiduals(m1), quantreg=T)
  # hist(ranef(m1)[[1]][,2] + fixef(m1)[2])
  
  
  coefs.m1 <- summary(m1)$coefficients[1:4,] %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m1 <- coefs.m1 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`),
    R2m = MuMIn::r.squaredGLMM(m1)[1],
    R2c = MuMIn::r.squaredGLMM(m1)[2]
  )
  
  CTI_trend.Rand <- cbind.data.frame(
    id_polygons = unique(CTI_yr$id_polygons),
    trend = ranef(m1)[[1]][,2] + fixef(m1)[2],
    n.occ = CTI_yr %>% group_by(id_polygons) %>% summarise(n.occ = sum(n.occ)) %>% pull(n.occ)
  )
  
  # analysis CTI trend ~ temperature trend

  CTI_trend <- CTI_trend.Rand %>% left_join(temp_windows %>% filter(Season == Season.df) %>%
                                              rename(trend.temp = trend, mean.temp = mean))%>%
    mutate(id_polygons = as.factor(id_polygons)) %>% 
    left_join(eco.windows)

  m2 <- lmer(
    trend ~ trend.temp + mean.temp + continent + (1|ECO_NAME),
    weights = n.occ,
    data = CTI_trend
  )
  
  # resid.sampl <- sample_n(residuals(m2) %>% as.data.frame() %>% rownames_to_column(), 1000)
  # sp.corel <- ncf::spline.correlog(x = CTI_yr[resid.sampl$rowname,"X"]$X,
  #                                  y = CTI_yr[resid.sampl$rowname,"Y"]$Y,
  #                                  z = resid.sampl[,2],
  #                                  resamp = 100, latlon = T, xmax = 1000)
  # plot(sp.corel ,main = paste("Taxon :", name.group, "- Model : m2"))
  
  coefs.m2 <- summary(m2)$coefficients[1:4,] %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m2 <- coefs.m2 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`),
    R2m = MuMIn::r.squaredGLMM(m2)[1],
    R2c = MuMIn::r.squaredGLMM(m2)[2]
  )
  
  # analysis CTM ~ hii
  CTA <- CTI_trend %>% mutate(CTA = abs(trend - trend.temp)) %>% left_join(hii_windows)

  m3 <- lmer(
    log(CTA) ~ hii + mean.temp + continent + (1|ECO_NAME),
      weights = n.occ,
    data = CTA
  )
  
  # resid.sampl <- sample_n(residuals(m3) %>% as.data.frame() %>% rownames_to_column(), 1000)
  # sp.corel <- ncf::spline.correlog(x = CTI_yr[resid.sampl$rowname,"X"]$X,
  #                                  y = CTI_yr[resid.sampl$rowname,"Y"]$Y,
  #                                  z = resid.sampl[,2],
  #                                  resamp = 100, latlon = T, xmax = 1000)
  # plot(sp.corel ,main = paste("Taxon :", name.group, "- Model : m2"))

  coefs.m3 <- summary(m3)$coefficients[1:4,] %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m3 <- coefs.m3 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`),
    R2m = MuMIn::r.squaredGLMM(m3)[1],
    R2c = MuMIn::r.squaredGLMM(m3)[2]
  )

  return(
    cbind.data.frame(
      Taxon = name.group,
      bind_rows(CTI_Year = coefs.m1,
                CTI_trend_Temp_trend = coefs.m2,
                CTA_HII = coefs.m3,
                .id = "Analysis"
      )
    )
  )
}

write_csv(res, "../../../../Results/results_cti_models.csv")


## mean + plots ##
library(metafor)

res <- read_csv("../../../../Results/results_cti_models.csv")

res <- res %>% 
  mutate(
    Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
    Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = as.factor(Taxon)
  ) %>% filter(Variable %in% c("year", "trend.temp", "hii")) 

taxa.col.order <- cbind.data.frame(Taxon = unique(res$Taxon))
taxa.col.order$col <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)

taxa.order <- res %>% filter(Analysis == "CTI_Year") %>% 
  arrange(Estimate) %>% pull(Taxon)

meta.res <- NULL
for(j in 1:length(unique(res$Analysis))){
  res.tmp <- res %>% filter(Analysis == unique(Analysis)[j])
  meta <- rma.uni(yi = res.tmp$Estimate, sei = res.tmp$`Std. Error`)
  meta.res <- rbind.data.frame(
    meta.res,
    cbind.data.frame(
      Analysis = unique(res$Analysis)[j],
      Variable = unique(res.tmp$Variable),
      Taxon = "MEAN", Estimate = meta$b, lwr.ci = meta$ci.lb, upr.ci = meta$ci.ub
    )
  )
}

res <- plyr::rbind.fill(res, meta.res)
res <- res %>% 
  mutate(Taxon = factor(Taxon, levels = c("MEAN",as.character(rev(taxa.order)))),
         Analysis = factor(Analysis, levels = c("CTI_Year", "CTI_trend_Temp_trend", "CTA_HII")))

res %>% 
  ggplot(aes(x = Estimate, xmin = lwr.ci, xmax = upr.ci, y = Taxon, color = Taxon)) +
  geom_vline(xintercept = 0, col = 'red3', linetype = 2) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_grid(.~Analysis, scale = 'free') +
  scale_y_discrete("") +
  scale_color_manual(values = c("brown",taxa.col.order[rev(taxa.order),"col"])) + 
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(colour = NA),
    legend.position = 'none',
    axis.text.y = element_text(colour = c("brown",taxa.col.order[rev(taxa.order),"col"]))
  )

ggsave(file = "../../../../Fig3.pdf", width = 12, height = 4)

### look at numbers ###
res %>% 
  filter(Taxon == "MEAN") %>% 
  select(1,2,4,9,10)

res %>% 
  arrange(desc(Estimate)) %>% 
  mutate_at(4:12, .funs = function(x){round(x,4)}) %>% 
  select(1,2,4,9,10) %>% group_by(Analysis) %>% group_split()

### test CTI trend in same countries as Lehikoinen et al  2021 ###
library(sf)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- world %>% filter(admin %in% c("Belgium", "United States of America", "Canada", "Spain", 
                                       "Denmark", "Estonia", "Finland", "Hungary", "Netherlands",
                                       "sweden"))

# winter
CTI_yr.winter <- read_csv(list.files(dir.df, pattern = "CTI_annual",full.names = T)[[4]]) %>% 
  left_join(temp_windows %>% filter(Season == "Winter")) %>% 
  left_join(read_csv(paste0(list.files(dir.df, pattern = "Descr_per_buffer",full.names = T)[[4]]))) %>% 
  mutate(id_polygons = as.factor(id_polygons), year = year - 1992) %>% 
  left_join(eco.windows)

CTI_yr.winter.sf <- st_as_sf(CTI_yr.winter, coords = c("X", "Y"), crs = st_crs(world)) 

CTI_yr.winter.sel <- st_intersection(CTI_yr.winter.sf, world)

m.w <- lmer(
  cti ~ year + mean + continent + 
    (year|id_polygons) + (1|ECO_NAME),
  weights = n.occ,
  data = CTI_yr.winter.sel
)
summary(m.w)
confint(m.w, method = "Wald")

# summer
CTI_yr.summer <- read_csv(list.files(dir.df, pattern = "CTI_annual",full.names = T)[[3]]) %>% 
  left_join(temp_windows %>% filter(Season == "Summer")) %>% 
  left_join(read_csv(paste0(list.files(dir.df, pattern = "Descr_per_buffer",full.names = T)[[3]]))) %>% 
  mutate(id_polygons = as.factor(id_polygons), year = year - 1992) %>% 
  left_join(eco.windows)

CTI_yr.summer.sf <- st_as_sf(CTI_yr.summer, coords = c("X", "Y"), crs = st_crs(world)) 

CTI_yr.summer.sel <- st_intersection(CTI_yr.summer.sf, world)

m.s <- lmer(
  cti ~ year + mean + continent + 
    (year|id_polygons) + (1|ECO_NAME),
  weights = n.occ,
  data = CTI_yr.summer.sel
)
summary(m.s)
confint(m.s, method = "Wald")

