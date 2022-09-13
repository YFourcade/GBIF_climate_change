#-------------------------------------#
# 9. Test influence of data quantity 
#-------------------------------------#

library(tidyverse)
library(lmerTest)
library(doFuture)
library(doRNG)
library(cowplot)
library(emmeans)

# load data sets
dir.df <- "../../../../Results/"
temp_windows <- read_csv("../../../../Results/temp_windows.csv")
hii_windows <- read_csv("../../../../Results/hii_windows.csv")
eco.windows <- read_csv("../../../../Results/ecoregions.csv")

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
  
  CTI_yr <- read_csv(list.files(dir.df, pattern = "CTI_annual",full.names = T)[[i]]) %>% 
    left_join(temp_windows %>% filter(Season == Season.df)) %>% 
    left_join(read_csv(paste0(list.files(dir.df, pattern = "Descr_per_buffer",full.names = T)[[i]]))) %>% 
    mutate(id_polygons = as.factor(id_polygons), year = year - 1992) %>% 
    left_join(eco.windows)
  
  m1.full <- lmer(
    cti ~ year + mean + continent + 
      (year|id_polygons) + (1|ECO_NAME), 
    weights = n.occ,
    data = CTI_yr
  )
  
  slopes.full <- cbind.data.frame(id_polygons = unique(CTI_yr$id_polygons), cti.trend.full = ranef(m1.full)[[1]][,2] + fixef(m1.full)[2])
  
  res.temp <- c()
  for(j in 1:nrow(param)){
    
    # data.qty.ori <- nrow(CTI_yr)
    # n.polygons.ori <- length(unique(CTI_yr$id_polygons))
    
    CTI_yr.tmp <- CTI_yr %>% filter(n.occ >= param[j,1]) %>% 
      group_by(id_polygons) %>% 
      mutate(n.yr = length(unique(year))) %>% 
      filter(n.yr >= param[j,2])
    
    if(length(unique(CTI_yr.tmp$continent)) > 1){
      m1 <- lmer(
        cti ~ year + mean + continent + 
          (year|id_polygons) + (1|ECO_NAME), 
        weights = n.occ,
        data = CTI_yr.tmp
      )
    } else {
      m1 <- lmer(
        cti ~ year + mean + 
          (year|id_polygons) + (1|ECO_NAME), 
        weights = n.occ,
        data = CTI_yr.tmp
      )
    }
    
    # coefs.m1 <- cbind.data.frame(
    #   Taxon = name.group,
    #   coef.m1,
    #   param[j,],
    #   data.qty = nrow(CTI_yr),
    #   n.polygons = length(unique(CTI_yr$id_polygons))
    # )
    
    CTI_trend.Rand <- cbind.data.frame(
      Taxon = name.group,
      param[j,],
      id_polygons = unique(CTI_yr.tmp$id_polygons),
      cti.trend = ranef(m1)[[1]][,2] + fixef(m1)[2],
      data.qty = nrow(CTI_yr.tmp),
      n.polygons = length(unique(CTI_yr.tmp$id_polygons))
    ) %>% left_join(slopes.full)
    
    CTI_trend.Overall <- cbind.data.frame(
      Taxon = name.group,
      param[j,],
      id_polygons = "Full",
      cti.trend = fixef(m1)[2],
      data.qty = nrow(CTI_yr.tmp),
      n.polygons = length(unique(CTI_yr.tmp$id_polygons)),
      cti.trend.full = fixef(m1.full)[2]
    )
    
    CTI_trend.Rand <- rbind.data.frame(CTI_trend.Rand, CTI_trend.Overall)
    
    res.temp <- rbind.data.frame(
      res.temp,
      CTI_trend.Rand
    )
    
  }
  return(res.temp)
}

write_csv(res, "../../../../Results/results_cti_models_testData_re.csv")

## plot - analyses ##
res.rep <- read_csv("../../../../Results/results_cti_models_testData_re.csv")

res.rep %>% left_join(hii_windows) %>% 
  left_join(temp_windows %>% filter(Season == "Winter")) %>% 
  filter(Taxon == "birds_winter") %>% 
  ggplot(aes(y = cti.trend.full, x = Y)) +
  geom_point() +
  theme_classic()

# overall estimates
range.CTI_Year.values <- res.rep %>% filter(id_polygons == "Full") %>% 
  select(cti.trend.full, cti.trend) %>% 
  pivot_longer(cols = 1:2) %>% pull(value) %>% range

p1 <- res.rep %>% filter(id_polygons == "Full") %>% 
  ggplot(aes(x = cti.trend.full, y = cti.trend)) +
  facet_grid(fct_rev(as.factor(n.occ.min)) ~ n.yr.min) + 
  geom_abline(slope = 1, intercept = 0, lwd = .5, color = "grey50") + 
  geom_hline(yintercept = 0, col = "grey", linetype = 3) + 
  geom_vline(xintercept = 0, col = "grey", linetype = 3) + 
  geom_point(aes(color = Taxon)) +
  scale_x_continuous("Overall CTI trend (full dataset)") +
  scale_y_continuous("Overall CTI trend (filtered dataset)") +
  coord_fixed(ylim = range.CTI_Year.values,
              xlim = range.CTI_Year.values) + 
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)) + 
  scale_shape_manual(values = c(15,16)) +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(colour = "black", fill = NA),
    panel.spacing = unit(0, "points"),
    plot.margin = unit(c(.5, .5, .5, .5), "lines"),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white")
  )


cor.text <- res.rep %>% filter(id_polygons == "Full") %>% group_by(n.occ.min, n.yr.min) %>%
  summarise(cor=cor.test(cti.trend, cti.trend.full, method = "pearson")$estimate,
            p=cor.test(cti.trend, cti.trend.full, method = "pearson")$p.value) %>% print(n=100
            )


p1 <- p1 + geom_text(
  size    = 2,
  data    = cor.text,
  mapping = aes(x = -Inf, y = Inf, label = paste("r =", round(cor,2))),
  hjust   = -.2,
  vjust   = 1.5
)



# random slopes
res.rep.slopes <- res.rep %>% filter(id_polygons != "Full") %>% 
  mutate(abs.diff = abs(cti.trend - cti.trend.full),
         n.occ.min = as.factor(n.occ.min),
         n.yr.min = as.factor(n.yr.min)
  ) %>%   mutate(
    Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
    Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = as.factor(Taxon)
  )


m.data <- lmer(
  log(abs.diff) ~ as.factor(n.yr.min) + as.factor(n.occ.min) + 
    Taxon +
    (1|id_polygons),
  data = res.rep.slopes)
# library(DHARMa)
# plot(simulateResiduals(m.data), quantreg=T)
anov.text.chisq <- car::Anova(m.data)
anov.text.F <- sjstats::anova_stats(m.data)

taxa.order <- emmeans(m.data, ~ Taxon, type = "response") %>% as.data.frame() %>% 
  arrange(response) %>% pull(Taxon)

taxa.col.order <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10)
names(taxa.col.order) <- unique(res.rep.slopes$Taxon)

emm_options(lmerTest.limit = 202397)

p2 <- emmeans(m.data, ~ Taxon, type = "response") %>% as.data.frame() %>% 
  mutate(
    Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
    Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = as.factor(Taxon)
  ) %>% 
  mutate(Taxon = factor(Taxon, levels = c(as.character(rev(taxa.order))))) %>% 
  ggplot(aes(x = Taxon, y = response, color = Taxon)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) +
  scale_x_discrete("") +
  scale_y_continuous("Abs. difference in CTI trend\nbetween filtered and full dataset") +
  scale_color_manual(values = taxa.col.order[rev(taxa.order)]) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = taxa.col.order[rev(taxa.order)], angle = 45, vjust = .5),
  ) +
  geom_text(  size    = 3,
              data    = cbind.data.frame(anov.text.chisq[3,], partial.etasq = anov.text.F[3,"partial.etasq"]),
              mapping = aes(x = -Inf, y = -Inf, 
                            label = paste0("Chi-squared = ", round(Chisq,3),"\n",
                                          "P = ", round(`Pr(>Chisq)`,3),"\n",
                                          "Eta-squared = ", round(partial.etasq,3))),
              hjust   = -.1,
              vjust   = -.5,
              color = 'black'
  )


p3 <- emmeans(m.data, ~ n.yr.min, type = "response") %>% as.data.frame() %>% 
  ggplot(aes(x = n.yr.min, y = response)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) +
  scale_x_discrete("Minimum no. of years\nby sliding window") +
  scale_y_continuous("Abs. difference in CTI trend\nbetween filtered and full dataset",
                     limits = c(0.008,0.025)) +
  theme_classic() +
  geom_text(  size    = 3,
              data    = cbind.data.frame(anov.text.chisq[1,], partial.etasq = anov.text.F[1,"partial.etasq"]),
              mapping = aes(x = -Inf, y = -Inf, 
                            label = paste0("Chi-squared = ", round(Chisq,3),"\n",
                                          "P = ", round(`Pr(>Chisq)`,3),"\n",
                                          "Eta-squared = ", round(partial.etasq,3))),
              hjust   = -.1,
              vjust   = -.5,
              color = 'black'
  )

p4 <- emmeans(m.data, ~ n.occ.min, type = "response") %>% as.data.frame() %>% 
  ggplot(aes(x = n.occ.min, y = response)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) +
  scale_x_discrete("Minimum no. of occurrences\nby sliding window") +
  scale_y_continuous("Abs. difference in CTI trend\nbetween filtered and full dataset",
                     limits = c(0.008,0.025)) +
  theme_classic() +
  geom_text(  size    = 3,
              data    = cbind.data.frame(anov.text.chisq[2,], partial.etasq = anov.text.F[2,"partial.etasq"]),
              mapping = aes(x = -Inf, y = -Inf, 
                            label = paste0("Chi-squared = ", round(Chisq,3),"\n",
                                          "P =", round(`Pr(>Chisq)`,3),"\n",
                                          "Eta-squared = ", round(partial.etasq,3))),
              hjust   = -.1,
              vjust   = -.5,
              color = 'black'
  )

p5 <- plot_grid(p1, plot_grid(p2, plot_grid(p3, p4), nrow = 2, rel_heights = c(1,.75)), rel_widths = c(.9,1))

ggsave(p5, filename = "../../../../Fig5.pdf")
