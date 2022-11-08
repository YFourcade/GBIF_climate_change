#-------------------------------------#
# 5. Analyses at the community level 
#-------------------------------------#

library(tidyverse)
library(lmerTest)
library(doFuture)
library(doRNG)
library(emmeans)

# load data sets
dir.df <- "../../../../Results/"
temp_windows <- read_csv("../../../../Results/temp_windows.csv")
hii_windows <- read_csv("../../../../Results/hii_windows.csv")
# quantile(hii_windows$hii, na.rm = T, probs = c(.1,.5,.9))
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

res <- foreach(i = 1:10) %dorng% {
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
    left_join(eco.windows) %>% mutate(CTA = abs(trend - trend.temp)) %>% left_join(hii_windows)
  
  m2 <- lmer(
    trend ~ trend.temp * (mean.temp + hii) + continent + (1|ECO_NAME),
    weights = n.occ,
    data = CTI_trend
  )
  
  cti_VS_temp_hii <- emtrends(
    m2, ~ trend.temp*hii, 
    var = "trend.temp", 
    at = list(hii = c(quantile(CTI_trend$hii, probs = .25), 
                      median(CTI_trend$hii), 
                      quantile(CTI_trend$hii, probs = .75)))
  ) %>% data.frame
  cti_VS_temp_hii$hii <- c("low", "median", "high")
  
  # resid.sampl <- sample_n(residuals(m2) %>% as.data.frame() %>% rownames_to_column(), 1000)
  # sp.corel <- ncf::spline.correlog(x = CTI_yr[resid.sampl$rowname,"X"]$X,
  #                                  y = CTI_yr[resid.sampl$rowname,"Y"]$Y,
  #                                  z = resid.sampl[,2],
  #                                  resamp = 100, latlon = T, xmax = 1000)
  # plot(sp.corel ,main = paste("Taxon :", name.group, "- Model : m2"))
  
  coefs.m2 <- summary(m2)$coefficients %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m2 <- coefs.m2 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`),
    R2m = MuMIn::r.squaredGLMM(m2)[1],
    R2c = MuMIn::r.squaredGLMM(m2)[2],
  )
  
  # # analysis CTM ~ hii
  # CTA <- CTI_trend %>% mutate(CTA = abs(trend - trend.temp)) %>% left_join(hii_windows)
  # 
  # m3 <- lmer(
  #   log(CTA) ~ hii + mean.temp + continent + (1|ECO_NAME),
  #     weights = n.occ,
  #   data = CTA
  # )
  # 
  # resid.sampl <- sample_n(residuals(m3) %>% as.data.frame() %>% rownames_to_column(), 1000)
  # sp.corel <- ncf::spline.correlog(x = CTI_yr[resid.sampl$rowname,"X"]$X,
  #                                  y = CTI_yr[resid.sampl$rowname,"Y"]$Y,
  #                                  z = resid.sampl[,2],
  #                                  resamp = 100, latlon = T, xmax = 1000)
  # plot(sp.corel ,main = paste("Taxon :", name.group, "- Model : m2"))
  # 
  #   coefs.m3 <- summary(m3)$coefficients[1:4,] %>% as.data.frame() %>% rownames_to_column("Variable")
  #   coefs.m3 <- coefs.m3 %>% mutate(
  #     lwr.ci = Estimate - (1.96*`Std. Error`),
  #     upr.ci = Estimate + (1.96*`Std. Error`),
  #     R2m = MuMIn::r.squaredGLMM(m3)[1],
  #     R2c = MuMIn::r.squaredGLMM(m3)[2]
  #   )
  
  return(
    list(
      model_results = cbind.data.frame(
        Taxon = name.group,
        bind_rows(CTI_Year = coefs.m1,
                  CTI_trend_Temp_trend = coefs.m2,
                  .id = "Analysis"
        )
      ),
      CTI_trend_var = cbind.data.frame(
        Taxon = name.group,
        CTI_trend.Rand),
      cti_VS_temp_hii = cbind.data.frame(
        Taxon = name.group,
        cti_VS_temp_hii)
    )
  )
}

res.models <- data.table::rbindlist(lapply(res, `[[`, 1))
cti.trends <- data.table::rbindlist(lapply(res, `[[`, 2))
trends.hii <- data.table::rbindlist(lapply(res, `[[`, 3))

write_csv(res.models, "../../../../Results/results_cti_models.csv")
write_csv(cti.trends, "../../../../Results/results_cti_ranef.csv")
write_csv(trends.hii, "../../../../Results/results_cti_temp_hii.csv")

## mean + plots ##
library(metafor)

res.all <- read_csv("../../../../Results/results_cti_models.csv") %>% 
  select(1,2,3,4,5,9,10) %>% filter(Variable == "year") %>% select(-Variable) %>% 
  mutate(hii = "median") %>% 
  rename(SE = `Std. Error`)
res.byHii <- read_csv("../../../../Results/results_cti_temp_hii.csv") %>% 
  select(-6,-2) %>% mutate(Analysis = "CTI_trend_Temp_trend") %>% 
  rename(Estimate = trend.temp.trend,
         lwr.ci = asymp.LCL,
         upr.ci = asymp.UCL)
res <- bind_rows(res.all, res.byHii)

res <- res %>% 
  mutate(
    Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
    Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = as.factor(Taxon)
  )

taxa.col.order <- cbind.data.frame(Taxon = unique(res$Taxon))
taxa.col.order$col <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)

taxa.order <- res %>% filter(Analysis == "CTI_Year") %>% 
  arrange(Estimate) %>% pull(Taxon)

meta.res <- NULL
for(j in unique(paste(res$Analysis, res$hii))){
  res.tmp <- res %>% filter(Analysis == strsplit(j, " ")[[1]][1], hii == strsplit(j, " ")[[1]][2])
  meta <- rma.uni(yi = res.tmp$Estimate, sei = res.tmp$SE)
  meta.res <- rbind.data.frame(
    meta.res,
    cbind.data.frame(
      Analysis = strsplit(j, " ")[[1]][1], hii = strsplit(j, " ")[[1]][2],
      Taxon = "MEAN", Estimate = meta$b, lwr.ci = meta$ci.lb, upr.ci = meta$ci.ub
    )
  )
}


res <- plyr::rbind.fill(res, meta.res)
res <- res %>% 
  mutate(Taxon = factor(Taxon, levels = c("MEAN",as.character(rev(taxa.order)))),
         Analysis = factor(Analysis, levels = c("CTI_Year", "CTI_trend_Temp_trend")),
         hii = factor(hii, levels = c("low", "median", "high")))

res %>% 
  ggplot(aes(x = Estimate, xmin = lwr.ci, xmax = upr.ci, y = Taxon, color = Taxon, shape = hii)) +
  geom_vline(xintercept = 0, col = 'red3', linetype = 2) +
  geom_point(position = position_dodge(.4)) +
  geom_errorbar(width = 0, position = position_dodge(.4)) +
  facet_grid(.~Analysis, scale = 'free') +
  scale_y_discrete("") +
  scale_color_manual(values = c("brown",taxa.col.order[rev(taxa.order),"col"])) + 
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(colour = NA),
    axis.text.y = element_text(colour = c("brown",taxa.col.order[rev(taxa.order),"col"]))
  )

ggsave(file = "../../../../Fig3.pdf", width = 8, height = 4)

### look at numbers ###
res.models <- read_csv("../../../../Results/results_cti_models.csv") 
  
res %>% 
  filter(Taxon == "MEAN")

res %>% filter(hii == 'median') %>% 
  arrange(desc(Estimate)) %>% 
  mutate_at(3:6, .funs = function(x){round(x,4)}) %>% 
  group_by(Analysis) %>% group_split()

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


### Display CTI trend by region ###
library(viridis)
world <- rnaturalearth::ne_coastline(returnclass = "sf")

trends.hii <- read_csv("../../../../Results/results_cti_ranef.csv") %>% 
  mutate(Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
         Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
         Taxon = gsub("Fourmis", "Formicidae", Taxon),
         Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
         Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
         Taxon = gsub("Rongeurs", "Rodentia", Taxon)
  )

ggplot() +
  geom_raster(data = trends.hii %>% left_join(temp_windows %>% filter(Season == "Annual") %>% select(2,6,7)) %>% 
                rename(`CTI trend` = trend),
              aes(x = X, y = Y, fill = `CTI trend`)) +
  geom_sf(data = world,
          color = "black", size = .01, fill = "white") +
  facet_wrap(. ~ Taxon, ncol = 3) +
  scale_fill_viridis(option = "inferno") +
  coord_sf(xlim = c(-170, 55), ylim = c(5, NA)) +
  theme_void() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.text.y = element_text(hjust = 0, margin = margin(1, 1, 1, 1, "lines")),
    strip.text.x = element_text(margin = margin(.5, .5, .5, .5, "lines")),
    legend.position = c(.7,.1),legend.direction = 'horizontal', legend.key.height = unit(.03, "npc")
  )

ggsave(file = "../../../../FigS2.pdf", width = 6, height = 5)


# trends.hii %>% left_join(temp_windows %>% filter(Season == "Annual") %>% select(2,4)) %>% 
#   rename(`CTI trend` = trend) %>% 
#   ggplot(aes(x = mean, y = `CTI trend`)) +
#   geom_point(size = .5) +
#   geom_smooth() +
#   facet_wrap(.~Taxon, ncol = 5) +
#   scale_x_continuous("Mean temperature (°C)") +
#   scale_y_continuous("CTI trend") +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA),
#   )
# 
# ggsave(file = "../../../../FigS2b.pdf", width = 6, height = 4)

