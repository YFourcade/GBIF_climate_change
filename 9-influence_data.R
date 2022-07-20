#-------------------------------------#
# 9. Test influence of data quantity 
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

param <- expand.grid(n.occ.min = c(50, 75, 100, 125, 150), n.yr.min = c(2,3,4,5))

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
  
  res.temp <- c()
  for(j in 1:nrow(param)){
    
    CTI_yr <- read_csv(list.files(dir.df, pattern = "CTI_annual",full.names = T)[[i]]) %>% 
      left_join(temp_windows %>% filter(Season == Season.df)) %>% 
      left_join(read_csv(paste0(list.files(dir.df, pattern = "Descr_per_buffer",full.names = T)[[i]]))) %>% 
      mutate(id_polygons = as.factor(id_polygons))
    
    data.qty.ori <- nrow(CTI_yr)
    n.polygons.ori <- length(unique(CTI_yr$id_polygons))
    
    CTI_yr <- CTI_yr %>% group_by(id_polygons) %>% filter(n.occ > param[j,1]) %>% 
      mutate(n.yr = length(unique(year))) %>% 
      filter(n.yr > param[j,2])
    
    if(length(unique(CTI_yr$continent)) > 1){
      m1 <- bam(
        cti ~ year + mean + 
          s(year, id_polygons, bs = "re") + 
          s(X, Y, bs = "sos"),
        weights = I(1/n.occ),
        data = CTI_yr,
      )
    } else {
      m1 <- bam(
        cti ~ year + mean +
          s(year, id_polygons, bs = "re"),
        weights = I(1/n.occ),
        data = CTI_yr,
      )
    }
    
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
    
    CTI_trend <- CTI_trend %>% filter(id_polygons %in% CTI_yr$id_polygons)
    
    
    if(length(unique(CTI_trend$continent)) > 1){
      m2 <- bam(
        trend ~ trend.temp + mean + continent +
          s(X, Y, bs = "sos"),
        weights = I(1/trend.se), 
        discrete = T, 
        nthreads = 4,
        data = CTI_trend %>% filter(trend.se > 0),
      )
    }else {
      m2 <- bam(
        trend ~ trend.temp + mean,
        weights = I(1/trend.se), 
        discrete = T, 
        nthreads = 4,
        data = CTI_trend %>% filter(trend.se > 0),
      )
    }
    coefs.m2 <- summary(m2, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
    coefs.m2 <- coefs.m2 %>% mutate(
      lwr.ci = Estimate - (1.96*`Std. Error`),
      upr.ci = Estimate + (1.96*`Std. Error`)
    )
    
    # analysis CTM ~ hii
    CTA <- CTI_trend %>% mutate(CTA = abs(trend - trend.temp)) %>% left_join(hii_windows)
    
    if(length(unique(CTI_trend$continent)) > 1){
      m3 <-bam(
        log(CTA) ~ hii + mean.temp + continent + 
          s(X, Y, bs = "sos"),
        data = CTA %>% filter(!is.na(CTA), trend.se > 0),
        weights = I(1/trend.se),
        discrete = T, 
        nthreads = 4
      )
    }else {
      m3 <-bam(
        log(CTA) ~ hii + mean.temp,
        data = CTA %>% filter(!is.na(CTA), trend.se > 0),
        weights = I(1/trend.se),
        discrete = T, 
        nthreads = 4
      )
    }
    coefs.m3 <- summary(m3, re.test = FALSE)$p.table %>% as.data.frame() %>% rownames_to_column("Variable")
    coefs.m3 <- coefs.m3 %>% mutate(
      lwr.ci = Estimate - (1.96*`Std. Error`),
      upr.ci = Estimate + (1.96*`Std. Error`)
    )
    
    res.temp <- rbind.data.frame(
      res.temp,
      cbind.data.frame(
        Taxon = name.group,
        param[j,],
        data.qty.ori,
        data.qty.filtered = nrow(CTI_yr),
        n.polygons.ori,
        n.polygons.filtered = length(unique(CTI_yr$id_polygons)),
        bind_rows(
          CTI_Year = coefs.m1,
          CTI_trend_Temp_trend = coefs.m2,
          CTA_HII = coefs.m3,
          .id = "Analysis"
        )
      )
    )
  }
  return(res.temp)
}

write_csv(res, "./Results/results_cti_models_testData.csv")

## plot - analyses ##
res.rep <- read_csv("./Results/results_cti_models_testData.csv")
res <- read_csv("./Results/results_cti_models.csv")

res.merged <- res.rep %>% left_join(
  res %>% rename(est.full = Estimate, lwr.full = lwr.ci, upr.full = upr.ci, t.full = `t value`) %>% select(c(1:4,6,8,9))
) %>% filter(Variable %in% c("mean", "trend.temp", "hii"), !(Variable == "mean" & Analysis == "CTI_trend_Temp_trend"))

res.merged %>% filter(Analysis == "CTI_Year") %>% select(1,2,3,10,11,12,13,14,15)

# res.merged %>% mutate(perc_diff = (Estimate - est.full) / est.full) %>% 
#   ggplot(aes(y = perc_diff, x = as.factor(n.yr.min), fill = as.factor(n.occ.min))) +
#   geom_abline(slope = 0, intercept = 0) +
#   geom_boxplot() +
#   facet_grid(Analysis~ ., scale = "free") +
#   theme_classic()

p1 <- res.merged %>% filter(Analysis == "CTI_Year") %>%
  mutate(Significant = ifelse(`Pr(>|t|)` < 0.05, "Significant", "Not significant")) %>% 
  ggplot(aes(x = est.full, y = Estimate, color = Significant)) +
  facet_grid(fct_rev(as.factor(n.occ.min)) ~ n.yr.min) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0, col = "grey", linetype = 3) + 
  geom_vline(xintercept = 0, col = "grey", linetype = 3) + 
  geom_point() +
  geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci),width = 0) +
  geom_errorbarh(aes(xmin = lwr.full, xmax = upr.full), height = 0) +
  scale_x_continuous("Coefficient estimate (full dataset)", limits = range(c(res.merged %>% filter(Analysis == "CTI_Year") %>% pull(est.full),
                                       res.merged %>% filter(Analysis == "CTI_Year") %>% pull(Estimate)))) +
  scale_y_continuous("Coefficient estimate (filtered dataset)", limits = range(c(res.merged %>% filter(Analysis == "CTI_Year") %>% pull(est.full),
                                       res.merged %>% filter(Analysis == "CTI_Year") %>% pull(Estimate)) )) +
  coord_fixed() +
  scale_color_manual(values = c("red3", "blue2")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white"),
    legend.position = "bottom"
  )



p2 <- res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>%
  mutate(Significant = ifelse(`Pr(>|t|)` < 0.05, "Significant", "Not significant")) %>% 
  ggplot(aes(x = est.full, y = Estimate, color = Significant)) +
  facet_grid(fct_rev(as.factor(n.occ.min)) ~ n.yr.min) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0, col = "grey", linetype = 3) + 
  geom_vline(xintercept = 0, col = "grey", linetype = 3) + 
  geom_point() +
  geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci),width = 0) +
  geom_errorbarh(aes(xmin = lwr.full, xmax = upr.full), height = 0) +
  scale_x_continuous("Coefficient estimate (full dataset)") +
  scale_y_continuous("Coefficient estimate (filtered dataset)", limits = range(c(res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>% pull(est.full),
                                       res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>% pull(Estimate)) )) +
  coord_fixed(ylim = range(c(res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>% pull(est.full),
                             res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>% pull(Estimate))),
              xlim = range(c(res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>% pull(est.full),
                             res.merged %>% filter(Analysis == "CTI_trend_Temp_trend") %>% pull(Estimate)))) + 
  scale_color_manual(values = c("red3", "blue2")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white"),
    legend.position = "bottom"
  )



p3 <- res.merged %>% filter(Analysis == "CTA_HII") %>%
  mutate(Significant = ifelse(`Pr(>|t|)` < 0.05, "Significant", "Not significant")) %>% 
  ggplot(aes(x = est.full, y = Estimate, color = Significant)) +
  facet_grid(fct_rev(as.factor(n.occ.min)) ~ n.yr.min) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0, col = "grey", linetype = 3) + 
  geom_vline(xintercept = 0, col = "grey", linetype = 3) + 
  geom_point() +
  geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci),width = 0) +
  geom_errorbarh(aes(xmin = lwr.full, xmax = upr.full), height = 0) +
  scale_x_continuous("Coefficient estimate (full dataset)", limits = range(c(res.merged %>% filter(Analysis == "CTA_HII") %>% pull(est.full),
                                       res.merged %>% filter(Analysis == "CTA_HII") %>% pull(Estimate)))) +
  scale_y_continuous("Coefficient estimate (filtered dataset)", limits = range(c(res.merged %>% filter(Analysis == "CTA_HII") %>% pull(est.full),
                                       res.merged %>% filter(Analysis == "CTA_HII") %>% pull(Estimate)) )) +
  coord_fixed() + 
  scale_color_manual(values = c("red3", "blue2")) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white"),
    legend.position = "bottom"
  )

p.full <- cowplot::plot_grid(p1,p2,p3, ncol = 3, labels = c("(a) CTI ~ Time", "(b) CTI trend ~ T° trend", "(c) CTM ~ HII"))

cowplot::ggsave2(plot = p.full, filename = "Fig5.pdf", height = 6, width = 18)

