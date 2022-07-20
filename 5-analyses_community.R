#-------------------------------------#
# 5. Analyses at the community level 
#-------------------------------------#

library(tidyverse)
library(mgcv)
library(doFuture)
library(doRNG)

setwd("C:/Users/200597/OneDrive - UPEC/Recherche/Students/Projets étudiants 2021/Armelle")

# load data sets
dir.df <- "./Results/"
temp_windows <- read_csv("./Results/temp_windows.csv")
hii_windows <- read_csv("./Results/hii_windows.csv")

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
    mutate(id_polygons = as.factor(id_polygons))
  
    m1 <- bam(
      cti ~ year + mean + continent + 
        s(year, id_polygons, bs = "re") + 
        s(X, Y, bs = "sos"),
      weights = I(1/n.occ),
      data = CTI_yr,
    )

  coefs.m1 <- summary(m1, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m1 <- coefs.m1 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`)
  )
  
  # analysis CTI trend ~ temperature trend
  CTI_trend <- read_csv(list.files(dir.df, pattern = "CTI_trend",full.names = T)[[i]]) %>% 
    left_join(temp_windows %>% filter(Season == Season.df) %>% 
                rename(trend.temp = trend, mean.temp = mean)) %>% 
    mutate(id_polygons = as.factor(id_polygons))
  
    m2 <- bam(
      trend ~ trend.temp + mean + continent +
        s(X, Y, bs = "sos"),
      weights = I(1/trend.se), 
      discrete = T, 
      nthreads = 4,
      data = CTI_trend %>% filter(trend.se > 0),
    )

  coefs.m2 <- summary(m2, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m2 <- coefs.m2 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`)
  )
  
  # analysis CTM ~ hii
  CTA <- CTI_trend %>% mutate(CTA = abs(trend - trend.temp)) %>% left_join(hii_windows)
  
    m3 <-bam(
      log(CTA) ~ hii + mean.temp + continent + 
        s(X, Y, bs = "sos"),
      data = CTA %>% filter(!is.na(CTA), trend.se > 0),
      weights = I(1/trend.se),
      discrete = T, 
      nthreads = 4
    )
  coefs.m3 <- summary(m3, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
  coefs.m3 <- coefs.m3 %>% mutate(
    lwr.ci = Estimate - (1.96*`Std. Error`),
    upr.ci = Estimate + (1.96*`Std. Error`)
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

write_csv(res, "./Results/results_cti_models.csv")


## mean + plots ##
library(metafor)

res <- read_csv("./Results/results_cti_models.csv")

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

ggsave(file = "Fig3.pdf", width = 12, height = 4)




