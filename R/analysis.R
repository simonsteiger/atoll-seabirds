box::use(
  dp = dplyr,
  tdr = tidyr,
  tbl = tibble,
  ts = tidyselect,
  veg = vegan,
  gg = ggplot2,
  magrittr[`%>%`],
  ggd = ggdensity,
  str = stringr,
)

pop <- read.csv("data/atoll_seabird_populations_11Mar.csv")
envs_jicp <- read.csv("data/envs_jicp.csv")
conds <- read.csv("data/seabird_filterconditions_03Jul.csv", sep = ";")

pop_recode <- pop %>% 
  dp$select(-starts_with("X")) %>% 
  dp$mutate(
     dp$across(c(ts$where(is.character), -atoll), as.numeric)
   ) %>% 
  tdr$pivot_longer(!atoll, names_to = "species", values_to = "presence") %>% 
  dp$mutate(
    presence = ifelse(!is.na(presence), TRUE, FALSE),
  )

#' @export
joined <- dp$left_join(pop_recode, envs_jicp, by = "atoll") %>% 
  tdr$pivot_wider(names_from = "species", values_from = "presence") %>% 
  dp$select(-X)

rda.env <- veg$rda(joined[,c(7:12,15,17:25)], weight = FALSE, scale = TRUE)

env.scores <- summary(rda.env)[2]$sites %>% as.data.frame()
env.scores$atoll <- joined$atoll
env.scores$region <- joined$region

species_only <- joined %>% 
  dp$select(ts$where(is.logical), atoll)

env.scores <- dp$left_join(tbl$as_tibble(env.scores), species_only, by = "atoll") %>% 
  tdr$pivot_longer(!c(PC1:PC6, atoll, region), names_to = "species", values_to = "presence") %>% 
  dp$left_join(conds, by = "species") %>% 
  dp$mutate(
    cond = dp$case_match(
      filtercondition,
      "EXCLUDE" ~ "\\d",
      "none" ~ "\\w+",
      .default = filtercondition
    ),
    cond = str$str_replace_all(cond, ",\\s", "|")
  ) %>% 
  dp$select(-filtercondition)

write.csv(env.scores, "data/envscores.csv")

env.loadings <- veg$scores(rda.env)$species %>% as.data.frame()

gg$ggplot(na.omit(env.scores), gg$aes(x = PC1, y = PC3)) +
  ggd$geom_hdr(aes(fill = presence, colour = presence)) +
  gg$geom_point(gg$aes(colour = presence), alpha = 0.5, size = 0.8) +
  gg$scale_color_manual(values = c("red", "blue", "white")) +
  gg$geom_segment(data = env.loadings, gg$aes(x = 0, xend = PC1, y=0, yend = PC2), 
               color = "black", alpha = 0.4, arrow = gg$arrow(length = gg$unit(0.01,"npc"))) +

  gg$geom_hline(yintercept=0, linetype="dotted") +
  gg$geom_vline(xintercept=0, linetype="dotted") +
  gg$theme_classic() +
  gg$facet_wrap(gg$vars(species))

testspecies <- env.scores %>% 
  dp$filter(species == "Puffinus_bailloni") %>% 
  dp$rename(text = atoll)

plot_ly(data = testspecies, hoverinfo = "text") %>% 
  add_trace(
    type = "scatter3d",
    mode = "markers",
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    text = ~text,
    color = ~presence,
    colors = c("blue", "red"),
    alpha_stroke = 0.5,
    alpha = 0.5,
    size = 0.8
  )

plot_ly(
  data = testspecies, 
  x = testspecies$PC1, 
  y = testspecies$PC2, 
  z = testspecies$PC3, 
  type = "scatter3d",
  mode = "markers", 
  color = testspecies$presence,
  colors = c("blue", "red"),
  alpha_stroke = 0.5,
  alpha = 0.5,
  size = 0.8
)

