box::use(
  dp = dplyr,
  tdr = tidyr,
  tbl = tibble,
  ts = tidyselect,
  veg = vegan,
  gg = ggplot2,
  magrittr[`%>%`],
)

pop <- read.csv("data/atoll_seabird_populations_11Mar.csv")
envs_jicp <- read.csv("data/envs_jicp.csv")

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

env.scores <- veg$scores(rda.env)$sites %>% as.data.frame()
env.scores$atoll <- joined$atoll

species_only <- joined %>% 
  dp$mutate(
    dp$
    Anous_tenuirostris = ifelse(region == "Indian_Ocean", as.character(Anous_tenuirostris), NA),
    Anous_minutus = ifelse(region != "Indian_Ocean", as.character(Anous_minutus), NA)
  ) %>% 
  dp$select(ts$where(is.logical), atoll)

env.scores <- dp$left_join(tbl$as_tibble(env.scores), species_only, by = "atoll") %>% 
  tdr$pivot_longer(!c(PC1, PC2, atoll), names_to = "species", values_to = "presence")

# fix later
# hull_data <- env.scores %>% dp$group_by(Sula_leucogaster) %>% dp$slice(chull(PC1, PC2)) %>% as.data.frame()

env.loadings <- veg$scores(rda.env)$species %>% as.data.frame()

gg$ggplot(na.omit(env.scores), gg$aes(x = PC1, y = PC2)) +
  gg$geom_point(gg$aes(colour = presence), alpha = 0.5, size = 0.8) +
  gg$scale_color_manual(values = c("red", "blue", "white")) +
  gg$geom_segment(data = env.loadings, gg$aes(x = 0, xend = PC1, y=0, yend = PC2), 
               color = "black", alpha = 0.4, arrow = gg$arrow(length = gg$unit(0.01,"npc"))) +
  # gg$geom_text(data = env.loadings, 
  #           gg$aes(x = PC1, y = PC2, label = rownames(env.loadings),
  #               hjust = 0.1*(1-sign(PC1)), vjust = 0.5*(1-sign(PC2)))) +
  # gg$geom_polygon(data = hull_data,
  #              gg$aes(fill = presence, colour = presence), alpha = 0.3, show.legend = FALSE) +
  gg$geom_hline(yintercept=0, linetype="dotted") +
  gg$geom_vline(xintercept=0, linetype="dotted") +
  gg$theme_classic() +
  gg$facet_wrap(gg$vars(species))

