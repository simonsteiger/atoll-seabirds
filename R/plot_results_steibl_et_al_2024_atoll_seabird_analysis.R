# This script is part of the project associated with
# Article: Atolls are globally significant sites for tropical seabirds
# Authors: Steibl S, Steiger S, Wegmann AS, Holmes ND, Young, HS, Carr P, Russell JC 
# Last edited: 2024-03-26

# Load packages ----

library("dplyr")
library("ggplot2")
library("tidyr")
library("stringr")
library("sf")
library("purrr")
library("patchwork")
library("rnaturalearth")
library("rnaturalearthdata")
library("RColorBrewer")
library("ggh4x")
library("scales")
library("fields")
library("here")
library("hdrcde")
library("adespatial")


suffix <- "_steibl_et_al_2024_atoll_seabird_analysis"

# Read data ----

# Load metadata. see supplementary file S1 for data source references
envs <- read.csv(here(paste0("data/vars_atoll", suffix, ".csv")))
nest <- read.csv(here(paste0("data/vars_seabird", suffix, ".csv")))
forag <- read.csv(here(paste0("data/vars_foraging", suffix, ".csv")))

# Load results file from Bayesian model analysis (julia-pipeline)
preds <- read.csv(here(paste0("results/data/pred_and_obs_atolls", suffix, ".csv")))
atollw <- read.csv(here(paste0("results/data/summary_count_atollwise", suffix, ".csv"))) # lower and upper is 95% HDI
specw <- read.csv(here(paste0("results/data/summary_count_specieswise", suffix, ".csv"))) # lower and upper is 80% HDI
nutr.atollw <- read.csv(here(paste0("results/data/summary_nutrient_atollwise", suffix, ".csv"))) # lower and upper is 95% HDI
nutr.specw <- read.csv(here(paste0("results/data/summary_nutrient_specieswise", suffix, ".csv"))) # lower and upper is 95% HID

# Prepare data frames for plotting ----

# convert longitudes for 'sf' format
envs$long <- ifelse(envs$long < 0, envs$long + 360, envs$long)

# turn into factors
missing <- envs %>% mutate(across(1:3, as.factor))

envs <- envs %>%
  select(atoll, region, lat, long, land_area_sqkm) %>%
  mutate(across(1:2, as.factor))
preds <- preds %>%
  rename(nbirds = median) %>%
  select(atoll, species, nbirds, lower, upper) %>%
  mutate(across(1:2, as.factor))
nest <- mutate(nest, across(1:4, as.factor))
forag$species <- as.factor(forag$species)
atollw$atoll <- as.factor(atollw$atoll)
specw$species <- as.factor(specw$species)
nutr.atollw$atoll <- as.factor(nutr.atollw$atoll)
nutr.specw$species <- as.factor(nutr.specw$species)

# turn predicted counts into whole numbers
preds <- preds %>% mutate(across(c(nbirds, lower, upper), round))
atollw <- atollw %>% mutate(across(c(median, lower, upper), round))
specw <- specw %>% mutate(across(c(median, lower, upper), round))

# drop unneeded columns from foraging dataset
forag <- forag %>% select(!c(location, reference))

# add 'empty' atolls to count data frames
dat <- left_join(envs, preds, by = "atoll") %>% mutate_at(c("nbirds", "lower", "upper"), ~ replace_na(., 0))
atollw <- left_join(envs, atollw, by = "atoll") %>% 
                mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% 
                select(!c("region", "land_area_sqkm")) %>% 
                rename_with(~paste0(.x, "_nbirds"), c(median, lower, upper))
nutr.atollw <- left_join(envs, nutr.atollw, by = "atoll") %>% 
                    mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% 
                    select(!c("region", "lat", "long", "land_area_sqkm"))
                  

# convert long format into atoll x species matrix for beta-diversity computation
matrix <- dat %>%
  select(!c("lower", "upper")) %>%
  pivot_wider(names_from = species, values_from = nbirds, values_fill = 0) %>%
  select(!"NA") %>%
  as.data.frame()

rm(dat) # no longer needed, keep tidy

# calculate local contribution to beta diversity (LCBD) [and species contribution to beta diversity (SCBD)]
beta.atoll <- matrix[c(which(rowMeans(matrix[, -c(1:5)]) > 0)), -c(1:5)] %>%
  beta.div(., method = "percentdiff", sqrt.D = TRUE, samp = FALSE)
beta <- beta.atoll$LCBD %>% as.data.frame()
beta$atoll <- matrix[c(which(rowMeans(matrix[, -c(1:5)]) > 0)), ]$atoll
colnames(beta)[1] <- "lcbd"
rm(beta.atoll) # no longer needed, keep tidy

# add seabird metadata 
comb <- matrix %>% pivot_longer(!c(atoll, region, lat, long, land_area_sqkm),
  names_to = "species", values_to = "nbirds"
)
comb$species <- as.factor(comb$species)
comb <- comb %>%
  left_join(., nest, by = "species") %>%
  as.data.frame()

# calculate total N and P input in per year and per seabird and atoll using
# https://doi.org/10.1038/s41467-017-02446-8 and https://doi.org/0.1016/j.envpol.2004.02.008

F_nc <- 0.036 # average N content of seabird diet in g N g^-1
F_pc <- 0.006 # average P content of seabird diet in g P g^-1
F_ec <- 6.5 # average energy content of seabird diet kJ g^-1
A_eff <- 0.8 # kJ obtained by the bird per kJ consumed

# ADULTS:
# calculate basal metabolic rate for adults (kJ bird^-1 day^-1)
comb$BMR <- 2.3 * comb$bodymass^0.774

# calculate field metabolic rate for adults (kJ bird^-1 day^-1)
comb$AMR <- 4 * comb$BMR

# calculate total daily N and P excreted (g N or g P bird^-1 day^-1)
comb$dailyN <- (comb$AMR * F_nc) / (F_ec * A_eff)
comb$dailyP <- (comb$AMR * F_pc) / (F_ec * A_eff)

# extrapolate for breeding season and colony attendance
comb$seasonalN <- comb$dailyN * comb$days_at_colony * comb$time_at_colony
comb$seasonalP <- comb$dailyP * comb$days_at_colony * comb$time_at_colony

# extrapolate for colony size
comb$adultN <- comb$seasonalN * comb$nbirds
comb$adultP <- comb$seasonalP * comb$nbirds

# CHICKS:
# calculate metabolic rate for chicks (kJ bird^-1 year^-1)
comb$E_rearing <- 28.43 * comb$fledgling_mass^1.06

comb$chickN <- ((comb$E_rearing * F_nc) / (F_ec * A_eff)) * (comb$prod_per_pair / 2)
comb$chickP <- ((comb$E_rearing * F_pc) / (F_ec * A_eff)) * (comb$prod_per_pair / 2)

comb$chickN <- comb$chickN * comb$nbirds
comb$chickP <- comb$chickP * comb$nbirds

# Add up adult and chick excretion
comb$excretedN <- comb$adultN + comb$chickN
comb$excretedP <- comb$adultP + comb$chickP

# convert g to kg
comb$excretedN <- comb$excretedN / 1000
comb$excretedP <- comb$excretedP / 1000


# calculate total NH3 emissions
F_nv <- 0.6 # proportion of volatilised excreted N
mass.ratio <- 17 / 14 # mass ratio of NH3 to N
comb$F_hab <- ifelse(comb$nestingtype == "burrow", 0, 0.2) # re-adsorption of NH3 by nesting type

comb$NH3emit <- comb$excretedN * F_nv * mass.ratio * comb$F_hab

# calculate total bird biomass
comb$biomass <- comb$nbirds * comb$bodymass

# calculate total carbon stock for birds using
# https://doi.org/10.1073/pnas.1711842115

wat.cont <- 0.7 # % water wet weight
dry.carb <- 0.5 # % carbon dry weight

comb$birdC <- comb$biomass * wat.cont * dry.carb


# drop unneeded columns after computation
comb <- comb %>% select(!c(
  bodymass, nestingtype, BMR, AMR, dailyN, dailyP, seasonalN, seasonalP, days_at_colony, time_at_colony,
  adultN, adultP, chickN, chickP, E_rearing, prod_per_pair, fledgling_mass, F_hab, common_names
))


# summarize atoll-wise: merge data frames, add species richness per atoll
atoll <- left_join(atollw, nutr.atollw, by = "atoll")

atoll <- atoll %>%
  left_join(., beta, by = "atoll") %>%
  mutate(lcbd = replace_na(lcbd, 0))

spec.ric <- comb %>% group_by(atoll) %>% summarise(nspec = sum(nbirds > 0)) %>% as.data.frame()
atoll <- atoll %>% 
  left_join(., spec.ric, by = "atoll")

# no longer needed, keep tidy
rm(spec.ric)
rm(beta)

# summarize species-wise:
# calculate proportion contribution to global population for lower (2.5%), median, and upper (97.5%) estimate
specw <- specw %>% 
  left_join(., nest, by = "species") %>% 
  select(c(species, common_names, group, lower, median, upper, global_est)) %>% 
  as.data.frame()

specw$ratio.lower <- specw$lower/specw$global_est
specw$ratio.median <- specw$median/specw$global_est
specw$ratio.upper <- specw$upper/specw$global_est

# set upper limit to 1
specw$ratio.upper[specw$ratio.upper > 1] <- 0.999



# Base map ----

mapWorld <- ne_countries(scale = "medium", returnclass = "sf")
mapCoast <- ne_coastline(scale = "medium", returnclass = "sf")

# create map layer
map <- ggplot() +
  geom_sf(data = mapWorld, fill = "#EBE3D5", color = NA) +
  geom_sf(data = mapCoast, linewidth = 0.15) +
  coord_sf(
    xlim = c(46.5, 255),
    ylim = c(-30, 30)
  ) +
  theme_classic() +
  geom_segment(aes(x = 36.2, xend = 36.2, y = -30.1, yend = 30.1), linewidth = 0.8) +
  geom_segment(aes(x = 49.9, xend = 250.1, y = -33, yend = -33), linewidth = 0.8) +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.8),
    title = element_text(size = 8),
    legend.text = element_text(size = 5),
    legend.background = element_rect(color = "black", fill = "white"),
    plot.title = element_text(margin = margin(b = -7))
  )

# set up base structure aesthetics for boxplots
boxplotaesth <- theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8)
  )

# Analysis of missingness ----

missing_theme <- theme(
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11)
  )

spat.missing <- missing %>%
  arrange(rev(seabird_data)) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)
coords.missing <- spat.missing %>% st_coordinates(extract())

map.missing <- map +
  geom_point(
    data = spat.missing,
    aes(x = coords.missing[, 1], y = coords.missing[, 2], colour = seabird_data),
    size = 1.1, shape = 21, alpha = 0.9, stroke = 1.9
  ) +
  scale_colour_discrete(name = "seabird data \navailable?")

pmiss1 <- ggplot(data = missing, aes(x = seabird_data, y = number_islets, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Number of islands per atoll") +
  xlab("Seabird data available?") +
  ggtitle("a. Island per atoll") +
  theme_classic() +
  guides(x = "axis_truncated", y = "axis_truncated") +
  missing_theme

pmiss2 <- ggplot(data = missing, aes(x = seabird_data, y = annual_precipitation_mm, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Annual precipitation [mm]") +
  xlab("Seabird data available?") +
  ggtitle("b. Annual precipitation") +
  theme_classic() +
  guides(x = "axis_truncated", y = "axis_truncated") +
  missing_theme

pmiss3 <- ggplot(data = missing, aes(x = seabird_data, y = land_area_sqkm, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Total atoll land area [km²]") +
  xlab("Seabird data available?") +
  ggtitle("c. Atoll land area") +
  scale_y_continuous(
    trans = "log10", labels = scales::comma,
    breaks = c(0.1, 1, 10, 100)
  ) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  missing_theme

pmiss4 <- ggplot(data = missing, aes(x = seabird_data, y = lagoon_area_sqkm, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Total atoll lagoon area [km²]") +
  xlab("Seabird data available?") +
  ggtitle("d. Atoll lagoon area") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  missing_theme

pmiss5 <- ggplot(data = missing, aes(x = seabird_data, y = tropical_storms_50km, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Number of tropical storms") +
  xlab("Seabird data available?") +
  ggtitle("e. Tropical storms") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  missing_theme

pmiss6 <- ggplot(data = missing, aes(x = seabird_data, y = hurricanes_50km, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Number of cyclones") +
  xlab("Seabird data available?") +
  ggtitle("f. Cyclones") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  missing_theme

pmiss7 <- ggplot(data = missing, aes(x = seabird_data, y = distance_nearest_atoll_km, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Distance to nearest atoll [km]") +
  ggtitle("g. Isolation: nearest atoll") +
  scale_y_continuous(
    trans = "log10", labels = scales::comma,
    breaks = c(10, 100, 1000)
  ) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  missing_theme

pmiss8 <- ggplot(data = missing, aes(x = seabird_data, y = distance_nearest_high_island_km, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Distance to nearest high island [km]") +
  xlab("Seabird data available?") +
  ggtitle("h. Isolation: nearest high island") +
  scale_y_continuous(
    trans = "log10", labels = scales::comma,
    breaks = c(10, 100, 1000)
  ) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

pmiss9 <- ggplot(data = missing, aes(x = seabird_data, y = distance_continent_km, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Distance to nearest continent [km]") +
  xlab("Seabird data available?") +
  ggtitle("i. Isolation: continent") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

pmiss10 <- ggplot(data = na.omit(missing), aes(x = seabird_data, y = human_population + 1, fill = seabird_data)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5) +
  ylab("Human population") +
  xlab("Seabird data available?") +
  ggtitle("j. Human population") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_y_continuous(
    trans = "log10", labels = scales::comma,
    breaks = c(10, 100, 1000, 10000)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

plotmissing <- (pmiss1 + pmiss2 + pmiss3) / (pmiss4 + pmiss5 + pmiss6) / (pmiss7 + pmiss8 + pmiss9) / (pmiss10 + ggplot() + ggplot())

# Save and export plot for analysis of missingness
# ggsave(map.missing, path = here("results", "svg", "article"), filename = "figS3_analysis_missing_map.svg", dpi = 300, width = 210, height = 100, units = "mm")
# ggsave(plotmissing, path = here("results", "svg", "article"), filename = "figS4_analysis_missing_boxplot.svg", dpi = 300, width = 210, height = 250, units = "mm")


# FIGURE 1 PLOTTING: abundance and richness ----

# total bird numbers
spat.1a <- atoll %>%
  arrange(median_nbirds) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)
coords.1a <- spat.1a %>% st_coordinates(extract())

p1a <- map +
  geom_point(
    data = spat.1a,
    aes(x = coords.1a[, 1], y = coords.1a[, 2], colour = log10(median_nbirds + 1)),
    size = 1.1, shape = 21, alpha = 0.9, stroke = 1.9
  ) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1.0, 0.5, 0.1, 0),
    name = "log(abundance)"
  ) +
  ggtitle("a. population size") +
  theme(legend.position = "none")

# species richness
spat.2a <- atoll %>%
  arrange(nspec) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)
coords.2a <- spat.2a %>% st_coordinates(extract())

p2a <- map +
  geom_point(
    data = spat.2a,
    aes(x = coords.2a[, 1], y = coords.2a[, 2], colour = nspec),
    size = 1.1, shape = 21, alpha = 0.9, stroke = 1.9
  ) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1.0, 0.3, 0.05, 0),
    name = "Species richness"
  ) +
  ggtitle("b. diversity") +
  theme(legend.position = "none")

# dissimilarity
spat.3a <- atoll %>% st_as_sf(., coords = c("long", "lat"), crs = 4326)
coords.3a <- spat.3a %>% st_coordinates(extract())

p3a <- map +
  geom_point(
    data = spat.3a,
    aes(x = coords.3a[, 1], y = coords.3a[, 2], colour = lcbd),
    size = 1.1, shape = 21, alpha = 0.9, stroke = 1.9
  ) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1, 0.65, 0.001, 0),
    name = "LCBD"
  ) +
  ggtitle("c. compositional distinctness") +
  theme(legend.position = "none")


# boxplot abundance
p1b <- ggplot(data = atoll[which(atoll$median_nbirds > 0), ], aes(x = "identity", y = median_nbirds)) +
  geom_violin(fill = "#EBEBEB", alpha = 0.5) +
  geom_jitter(aes(colour = log(median_nbirds + 1)), width = 0.15, size = 1.5) +
  geom_boxplot(fill = "#EBEBEB", width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(
    trans = "log",
    breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000),
    labels = scales::comma
  ) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1.0, 0.5, 0.1, 0)
  ) +
  guides(x = "axis_truncated", y = guide_axis_truncated(trunc_lower = 10, trunc_upper = 1000000)) +
  ylab("abundance") +
  geom_hline(yintercept = 13200, linetype = "dotted", linewidth = 0.5) +
  annotate("text", label = "B3b (A4iii)", x = 1.4, y = 25000, size = 2, fontface = "italic") +
  boxplotaesth

# boxplot richness
p2b <- ggplot(data = atoll, aes(x = "identity", y = nspec)) +
  geom_violin(fill = "#EBEBEB", alpha = 0.5) +
  geom_jitter(aes(colour = nspec), width = 0.15, size = 1.5) +
  geom_boxplot(fill = "#EBEBEB", width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1.0, 0.3, 0.05, 0)
  ) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 19)) +
  guides(x = "axis_truncated", y = guide_axis_truncated(trunc_lower = 0, trunc_upper = 19)) +
  ylab("species richness") +
  boxplotaesth

# boxplot lcbd
p3b <- ggplot(data = atoll[which(atoll$median_nbirds > 0), ], aes(x = "identity", y = lcbd)) +
  geom_violin(fill = "#EBEBEB", alpha = 0.5) +
  geom_jitter(aes(colour = lcbd), width = 0.15, size = 1.5) +
  geom_boxplot(fill = "#EBEBEB", width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1, 0.65, 0.001, 0)
  ) +
  scale_y_continuous(breaks = c(0.0026, 0.0037, 0.0048, 0.0059)) +
  guides(x = "axis_truncated", y = guide_axis_truncated(trunc_lower = 0.0026, trunc_upper = 0.0059)) +
  ylab("LCBD") +
  boxplotaesth

pmapbox <- p1a + p1b + p2a + p2b + p3a + p3b + plot_layout(widths = c(3, 1, 3, 1, 3, 1), nrow = 3, ncol = 2)

# Save and export Figure 1
ggsave(pmapbox, path = here("results", "svg", "article"), filename = "fig01_pop-ric-lcbd.svg", dpi = 300, width = 180, height = 146, units = "mm")


# FIGURE 2 PLOTTING: bird populations ----

specw$range <- ifelse(specw$ratio.median >= 0.95, ">95%", 
                          ifelse(specw$ratio.median >= 0.75, ">75%", 
                                 ifelse(specw$ratio.median >= 0.50, ">50%",
                                        ifelse(specw$ratio.median >= 0.25, ">25%",
                                        "<25%")))) %>% as.factor()


pglobal <- ggplot(data = specw) +
  geom_bar(aes(x = ratio.median * 100, y = reorder(common_names, ratio.median), fill = range),
    stat = "identity", color = "black", width = 0.8, linewidth = 0.5
  ) +
  annotate("rect", xmin = 25, xmax = 50, ymin = 0, ymax = 37, fill = "#eeeeee") +
  annotate("rect", xmin = 50, xmax = 75, ymin = 0, ymax = 37, fill = "#dddddd") +
  annotate("rect", xmin = 75, xmax = 95, ymin = 0, ymax = 37, fill = "#cccccc") +
  annotate("rect", xmin = 95, xmax = 100, ymin = 0, ymax = 37, fill = "#bbbbbb") +
  geom_bar(aes(x = ratio.median * 100, y = reorder(common_names, ratio.median), fill = range),
    stat = "identity", color = "black", width = 0.8, linewidth = 0.5
  ) +
  geom_errorbarh(aes(xmin = ratio.lower * 100, xmax = ratio.upper * 100, y = reorder(common_names, ratio.median)),
    linewidth = 0.7, height = 0.25
  ) +
  scale_fill_viridis_d(
    aesthetics = "fill", option = "mako", direction = -1
  ) +
  guides(
    x = "axis_truncated", y = "axis_truncated",
    fill = guide_legend(byrow = TRUE)
  ) +
  xlab("% of global population") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(margin = margin(r = -6), size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    legend.position = "none"
  )

# Save and export Figure 2
ggsave(pglobal, path = here("results", "svg", "article"), filename = "fig02_proportion_species_global.svg", dpi = 300, width = 88, height = 165, units = "mm")



# FIGURE 3 AND S2 PLOTTING: nutrient inputs ----

pnitrog.a <- ggplot(data = atoll[which(atoll$median_nbirds > 0), ], aes(x = "identity", y = median_excretedN)) +
  geom_violin(fill = "#EBEBEB", alpha = 0.5) +
  geom_jitter(aes(colour = median_excretedN), width = 0.15, size = 1.5) +
  geom_boxplot(fill = "#EBEBEB", width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(
    trans = "log10", labels = scales::comma,
    breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)
  ) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1, 0.001, 0.0001, 0)
  ) +
  ylab(expression(paste("imported nitrogen per atoll [kg N  ", year^-1, "]"))) +
  ggtitle("a. Nitrogen import per atoll") +
  guides(y = guide_axis_truncated(trunc_lower = 1, trunc_upper = 1000000)) +
  boxplotaesth +
  theme(axis.line.x = element_blank(),
        plot.title = element_text(size = 8))

pphosph.a <- ggplot(data = atoll[which(atoll$median_nbirds > 0), ], aes(x = "identity", y = median_excretedP)) +
  geom_violin(fill = "#EBEBEB", alpha = 0.5) +
  geom_jitter(aes(colour = median_excretedP), width = 0.15, size = 1.5) +
  geom_boxplot(fill = "#EBEBEB", width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(
    trans = "log10", labels = scales::comma,
    breaks = c(1, 10, 100, 1000, 10000, 100000, 100000)
  ) +
  scale_colour_gradientn(
    colours = c("#fbe725", "#3dbc74", "#306b84", "black"),
    values = c(1, 0.001, 0.0001, 0)
  ) +
  ylab(expression(paste("imported phosphorous per atoll [kg P ", year^-1, "]"))) +
  ggtitle("a. phosphorous import per atoll") +
  guides(y = guide_axis_truncated(trunc_lower = 1, trunc_upper = 100000)) +
  boxplotaesth +
  theme(axis.line.x = element_blank(),
        plot.title = element_text(size = 8))


# Summarise nutrient-level input by species group
colonywise <- comb %>%
  group_by(atoll, group) %>%
  summarise(
    totN = sum(excretedN),
    totP = sum(excretedP)
  ) %>%
  as.data.frame()
colonywise$group <- factor(colonywise$group,
  levels = c("tropicbird", "tubenose", "frigatebird", "tern", "booby", "albatros")
)


pnitrog.b <- ggplot(
  data = colonywise[which(colonywise$totN > 0), ],
  aes(y = group, x = totN, fill = group)
) +
  geom_boxplot() +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000),
    labels = ~ ifelse(.x < 1, label_number(accuracy = .1)(.x),
      label_number(accuracy = 1)(.x)
    )
  ) +
  scale_y_discrete(
    labels = c("Tropicbirds", "Petrels &\n Shearwaters", "Frigatebirds", "Terns", "Boobies", "Albatrosses")
  ) +
  scale_fill_manual(values = c(
    "#DEF5E5FF", "#3497A9FF", "#395D9CFF",
    "#60CEACFF", "#382A54FF", "#0B0404FF"
  )) +
  theme_classic() +
  xlab(expression(paste("imported nitrogen per atoll-colony [kg N ", year^-1, "]"))) +
  guides(x = guide_axis_truncated(trunc_lower = 0.1, trunc_upper = 1000000), y = "axis_truncated") +
  ggtitle("b. keystone importers") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7, angle = 30, vjust = 0.8),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 8)
  )

pphosph.b <- ggplot(
  data = colonywise[which(colonywise$totP > 0), ],
  aes(y = group, x = totP, fill = group)
) +
  geom_boxplot() +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000),
    labels = ~ ifelse(.x < 1, label_number(accuracy = .1)(.x),
      label_number(accuracy = 1)(.x)
    )
  ) +
  scale_y_discrete(
    labels = c("Tropicbirds", "Petrels &\n Shearwaters", "Frigatebirds", "Terns", "Boobies", "Albatrosses")
  ) +
  scale_fill_manual(values = c(
    "#DEF5E5FF", "#3497A9FF", "#395D9CFF",
    "#60CEACFF", "#382A54FF", "#0B0404FF"
  )) +
  theme_classic() +
  xlab(expression(paste("imported phosphorous per atoll-colony [kg P ", year^-1, "]"))) +
  guides(x = guide_axis_truncated(trunc_lower = 0.1, trunc_upper = 100000), y = "axis_truncated") +
  ggtitle("b. keystone importers") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7, angle = 30, vjust = 0.8),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 8)
  )

pnutrientsN <- pnitrog.a + pnitrog.b + plot_layout(widths = c(1, 1.3)) # figure 3
pnutrientsP <- pphosph.a + pphosph.b + plot_layout(widths = c(1, 1.3)) # figure S2

# Save and export figure 3 and S2
 ggsave(pnutrientsN, path = here("results", "svg", "article"), filename = "fig03_nitrogen_import.svg", dpi = 300, width = 180, height = 72, units = "mm")
 ggsave(pnutrientsP, path = here("results", "svg", "article"), filename = "figS2_phosphorous_import.svg", dpi = 300, width = 180, height = 72, units = "mm")


# FIGURE S1: significant breeding sites ----

sitew <- comb %>% select(c(atoll, region, lat, long, species, nbirds, group, global_est))

sitew$atoll.ratios <- sitew$nbirds/sitew$global_est

sitew$species <- as.factor(sitew$species)


# subset only for largest value per atoll
atollmax <- sitew %>%
  group_by(atoll, lat, long) %>%
  summarise(atoll.ratios = max(atoll.ratios)) %>%
  as.data.frame()
atollmax$cuts <- cut(atollmax$atoll.ratios,
  breaks = c(0, 0.01, 0.10, 0.25, 0.50, 0.70, 1),
  labels = c("<1%", "1-10%", "10-25%", "25-50%", "50-70%", "70%")
)

atollmax[is.na(atollmax$cuts), ]$cuts <- "<1%"

atollmax <- atollmax %>% arrange(cuts)

spatw <- st_as_sf(atollmax, coords = c("long", "lat"), crs = 4326)
coordsw <- spatw %>% st_coordinates(extract())

psignif <- map +
  geom_point(
    data = spatw,
    aes(x = coordsw[, 1], y = coordsw[, 2], colour = cuts),
    size = 1.1, shape = 21, alpha = 0.9, stroke = 1.9
  ) +
  scale_colour_manual(values = c("#ECEE81", "#7ED7C1", "#BEADFA", "#DC8686", "#B06161", "#D10015")) +
  ggtitle("Single atoll contribution to global populations") +
  theme(
    legend.title = element_blank()
  )

# Save and export Figure S1
# ggsave(psignif, path = here("results", "svg", "article"), filename = "figS1_atollw_contribution_NO-ICONS.svg", dpi = 300, width = 210, height = 100, units = "mm")

# FIGURE 4 (PART): foraging distance plots ----

# summarize species foraging distance
foragdist <- forag %>%
  group_by(species) %>%
  summarise(
    dist = mean(maximum_distance_km),
    distvar = sd(maximum_distance_km)
  ) %>%
  as.data.frame()


# Take Hydrobates sp. values for Hydrobates tristrami
foragdist <- foragdist %>% mutate(species = str_replace(species, "Hydrobates_sp", "Hydrobates_tristrami"))


# Take Hydrobates sp. values for Nesofregetta fuliginosa
foragdist <- foragdist %>% add_row(
  species = "Nesofregetta_fuliginosa",
  dist = foragdist[which(foragdist$species == "Hydrobates_tristrami"), "dist"],
  distvar = foragdist[which(foragdist$species == "Hydrobates_tristrami"), "distvar"]
)


# replace Sternula sp. values for Sternula saundersi
foragdist <- foragdist %>% mutate(species = str_replace(species, "Sternula_sp", "Sternula_saundersi"))


# use genus-level variance for species with NA in vardist
foragdist[which(foragdist$species == "Anous_minutus"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Anous")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Hydroprogne_caspia"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Thalasseus")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Phoebastria_albatrus"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Phoebastria")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Pterodroma_alba"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Pterodroma")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Pterodroma_heraldica"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Pterodroma")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Puffinus_bailloni"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Puffinus")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Sterna_sumatrana"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Sterna")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Sternula_saundersi"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Sternula")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Sternula_nereis"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Sternula")), "maximum_distance_km"])

foragdist[which(foragdist$species == "Thalasseus_bengalensis"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Thalasseus")), "maximum_distance_km"])



# Add foraging distances to main data file for computation ('comb')
comb <- comb %>% left_join(., foragdist, by = "species")

# add missing species as mean of genus-level information
comb[which(comb$species == "Onychoprion_lunatus"), "dist"] <- mean(forag[which(str_detect(forag$species, "Onychoprion")), "maximum_distance_km"])
comb[which(comb$species == "Onychoprion_lunatus"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Onychoprion")), "maximum_distance_km"])

comb[which(comb$species == "Pterodroma_hypoleuca"), "dist"] <- mean(forag[which(str_detect(forag$species, "Pterodroma")), "maximum_distance_km"])
comb[which(comb$species == "Pterodroma_hypoleuca"), "distvar"] <- sd(forag[which(str_detect(forag$species, "Pterodroma")), "maximum_distance_km"])

# create foraging distance distribution
avgforag.atoll <- comb %>%
  select(atoll, nbirds, dist, distvar) %>%
  filter(nbirds > 0)

out <- map(seq_len(nrow(avgforag.atoll)), \(i) {
  data.frame(
    atoll = avgforag.atoll[i, ]$atoll,
    dist = rnorm(
      n = avgforag.atoll[i, ]$nbirds,
      mean = avgforag.atoll[i, ]$dist,
      sd = avgforag.atoll[i, ]$distvar
    )
  )
}) %>% list_rbind()

out$dist <- abs(out$dist)

mean.forag <- out %>%
  group_by(atoll) %>%
  summarise(
    mean = mean(dist),
    sd = sd(dist),
    nbirds = length(atoll),
    low.prctile = hdr(dist, prob = 50)$hdr[, 1] %>% as.vector(),
    upp.prctile = hdr(dist, prob = 50)$hdr[, 2] %>% as.vector()
  ) %>%
  as.data.frame()

mean.forag$area <- (pi * mean.forag$upp.prctile^2) - (pi * mean.forag$low.prctile^2)

p.forag <- ggplot() +
  geom_density(
    data = mean.forag,
    aes(x = area),
    fill = "#D2E0FB", linewidth = 0.3, alpha = 0.5
  ) +
  scale_x_continuous(
    trans = "log",
    breaks = c(100, 1000, 10000, 100000, 1000000, 10000000),
    labels = scales::comma
  ) +
  theme_classic() +
  xlab("km²") +
  ylab("Probability density") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, angle = 30, vjust = 0.8, hjust = 0.8),
    axis.title = element_text(size = 9),
    axis.text.y = element_blank()
  )

# Save and export Figure 4 (part)
# ggsave(p.forag, path = here("results", "svg", "article"), filename = "foraging.svg", dpi = 300, width = 575, height = 438, units = "px")


# Text-based results summary ----

# Individual atoll nesting populations ranging from ... to ...
min(atoll$median_nbirds)
max(atoll$median_nbirds)

# the mean nesting population per atoll is
mean(atoll$median_nbirds)
c(mean(atoll$lower_nbirds), mean(atoll$upper_nbirds)) # 95% HDI

# Number of atolls that house a colony for a given species above the B3b IBA threshold of 13,200 birds
# http://datazone.birdlife.org/site/ibacritreg
atoll %>%
  filter(median_nbirds >= 13200) %>%
  distinct(atoll, .keep_all = TRUE) %>%
  nrow()

# Number of atolls with a seabird colony that constitutes relative to the estimated global population more than
# >1%
atollmax %>%
  filter(cuts != "<1%") %>%
  nrow()
# >10%
atollmax %>%
  filter(cuts != "<1%" & cuts != "1-10%") %>%
  nrow()
# >70%
atollmax %>%
  filter(cuts == "70%") %>%
  nrow()

# Biomass of all nesting seabirds per atoll combined
# (divide by water content and dry carb content to extrapolate from C to tot biomass, see Bar-On et al. 2018)
mean(atoll$median_bird_C/(wat.cont * dry.carb))/1000
c(mean(atoll$lower_bird_C/(wat.cont * dry.carb))/1000,
  mean(atoll$upper_bird_C/(wat.cont * dry.carb))/1000)

# Atollwise seabird carbon stock
mean(atoll$median_bird_C) / 1000
c(mean(atoll$lower_bird_C/1000),
  mean(atoll$upper_bird_C/1000))

# Number of species where on atolls are nesting of the relative global population more than
# 25%
specw %>%
  filter(ratio.median > 0.25) %>%
  nrow()
# 50%
specw %>%
  filter(ratio.median > 0.50) %>%
  nrow()
# 75%
specw %>%
  filter(ratio.median > 0.75) %>%
  nrow()
# 95%
specw %>%
  filter(ratio.median > 0.95) %>%
  nrow()

# Seabird colonies on atolls import on average nitrogen in quantities of
mean(atoll$median_excretedN)
c(mean(atoll$lower_excretedN),
  mean(atoll$upper_excretedN))

# Seabird colonies on atolls import on average phosphorous in quantities of
mean(atoll$median_excretedP)
c(mean(atoll$lower_excretedP),
  mean(atoll$upper_excretedP))

# Estimated average global ammonia emissions from seabirds per atoll
mean(atoll$median_NH3emit)
c(mean(atoll$lower_NH3emit),
  mean(atoll$upper_NH3emit))
