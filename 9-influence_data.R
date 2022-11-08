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
    weights = log(n.occ),
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
        weights = log(n.occ),
        data = CTI_yr.tmp
      )
    } else {
      m1 <- lmer(
        cti ~ year + mean + 
          (year|id_polygons) + (1|ECO_NAME), 
        weights = log(n.occ),
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

# overall estimates
range.CTI_Year.values <- res.rep %>% filter(id_polygons == "Full") %>% 
  select(cti.trend.full, cti.trend) %>% 
  pivot_longer(cols = 1:2) %>% pull(value) %>% range

res.rep %>% filter(id_polygons == "Full") %>% mutate(diff = cti.trend - cti.trend.full) %>% 
  pull(diff) %>% hist

res.rep %>% filter(id_polygons != "Full") %>% mutate(diff = cti.trend - cti.trend.full) %>% mutate(
  Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
  Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
  Taxon = gsub("Fourmis", "Formicidae", Taxon),
  Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
  Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
  Taxon = gsub("Rongeurs", "Rodentia", Taxon),
  Taxon = gsub("Caudata", "Urodela", Taxon)
) %>% 
  ggplot(aes(x = diff)) +
  geom_histogram(color = 'white', fill = "lightgrey") +
  facet_wrap(.~Taxon, scales = 'free', ncol = 3) + 
  geom_vline(xintercept = 0, color = "red4") +
  scale_x_continuous("Difference in CTI trend between filtered and full dataset") +
  scale_y_continuous("Frequency") +
  theme_minimal()


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
  ) %>%   mutate(
    Taxon = gsub("birds_summer", "Aves (summer)", Taxon),
    Taxon = gsub("birds_winter", "Aves (winter)", Taxon),
    Taxon = gsub("Fourmis", "Formicidae", Taxon),
    Taxon = gsub("lepidoptera", "Lepidoptera", Taxon),
    Taxon = gsub("Lombrics", "Lumbricidae", Taxon),
    Taxon = gsub("Rongeurs", "Rodentia", Taxon),
    Taxon = gsub("Caudata", "Urodela", Taxon),
    Taxon = as.factor(Taxon)
  )


m.data <- lmer(
  log(abs.diff) ~ n.yr.min * n.occ.min * Taxon +
    (1|id_polygons),
  data = res.rep.slopes)
# library(DHARMa)
# plot(simulateResiduals(m.data), quantreg=T)
anov.text.chisq <- car::Anova(m.data, type = 3)
anov.text.F <- sjstats::anova_stats(m.data)

test(emtrends(m.data, spec = ~n.yr.min * Taxon, var = "n.yr.min"))
test(emtrends(m.data, spec = ~n.occ.min * Taxon, var = "n.occ.min"))

newdat <- expand.grid(Taxon = unique(res.rep.slopes$Taxon),
                      n.occ.min = seq(50,150,length.out = 5),
                      n.yr.min = seq(2,5,length.out = 5))

pred <- cbind.data.frame(
  newdat,
  pred = predict(m.data, newdata = newdat, re.form = NA)
) %>%  mutate(pred = exp(pred))


p2 <- pred %>% ggplot(aes(x = n.occ.min, y = n.yr.min, fill = pred)) +
  geom_raster() +
  facet_wrap(~Taxon, ncol = 3) +
  theme_classic() +
  scale_x_continuous("Minimum no. of occurrences by sliding window") +
  scale_y_continuous("Minimum no. of years by sliding window") +
  scale_fill_viridis_c(option = "magma") +
  theme(
    strip.background = element_blank(),
    legend.position = c(.7,.1),legend.direction = 'horizontal', legend.key.height = unit(.03, "npc")
  )

p3 <- plot_grid(p1, p2, align = "hv", axis = 'tb', labels = c("(a)", "(b)"))

ggsave(p3, filename = "../../../../Fig5.pdf")
