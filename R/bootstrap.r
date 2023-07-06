library(dplyr)
library(tibble)
library(purrr)
library(tidyr)

data <- envscores

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

table(x$presence)


a <- as.numeric(summary(data_as$presence)[[2]])
n <- nrow(data_as)
N <- 1000

a/(n-a)


N_a <- (N/2)-a # this is how many additional (resampled) absence we need
N_p <- (N/2)-(n-a) # this is how many additional (resampled) presence we need

out <- c(
          row_number(data_as),
          sample(which(data_as$presence == 0), size = N_a, replace = TRUE),
          sample(which(data_as$presence == 1), size = N_p, replace = TRUE)
)

res <- tibble()

for (i in seq_along(out)) {
  temp <- data_as[out[i]]
  res <- rbind(res, temp)
}

table(res$presence)
