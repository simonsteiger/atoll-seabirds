library(dplyr)
library(tibble)
library(purrr)
library(tidyr)

data <- tibble(read.csv("/Users/simonsteiger/Desktop/other/atoll-seabirds/data/envscores.csv"))

data_as <- data %>% 
  filter(species == "Anous_stolidus")

N <- 1000

bootstrap <- function(df, N, target, r_out = 1) {
  df <- df %>% 
    group_by(.data[[target]]) %>% 
    mutate(prob = n()/nrow(df),
           revprob = ifelse(presence == 0, 95.6, 0.01))
  
  x <- c(row_number(df),
         sample(row_number(df), nrow(df), replace = TRUE, prob = df$revprob),
         sample(row_number(df), N-(nrow(df)*2), replace = TRUE, prob = r_out - df$prob)
        )
  
  res <- tibble()
  
  for (i in seq_along(x)) {
    temp <- df[x[i], ]
    res <- rbind(res, temp)
  }
  
  return(res)
}

x <- bootstrap(data_as, N, "presence")