library(ncdf4)
library(lubridate)
library(ggplot2)

nc_file <- nc_open("data/copernicus_nppv.nc")
envs <- read.csv("data/seabird_atolls_envs_10Mar.csv")

# fix negative longitudes
envs <- envs %>% 
  mutate(
    long = ifelse(long < 0, long + 360, long)
  )

latitude <- as.vector(nc_file$dim[[3]]$vals)
longitude <- as.vector(nc_file$dim[[4]]$vals)

x <- ncvar_get(nc_file, nc_file$var[[1]])

matrix <- apply(x, c(1, 2), mean)

colnames(matrix) <- latitude

cp_data <- tibble::as_tibble(matrix) %>%
  bind_cols(longitude = longitude) %>% 
  tidyr::pivot_longer(-longitude, names_to = "latitude", values_to = "nppv") %>% 
  mutate(
    latitude = as.numeric(latitude)
    )

filter_copernicus <- function(cp_data, lat, long) {
  cp_data %>% 
    dplyr::filter(
      latitude >= lat - 1,
      latitude <= lat + 1,
      longitude >= long - 1,
      longitude <= long + 1
    ) %>% 
    summarise(
      mean_nppv = mean(nppv, na.rm = TRUE),
      sd_nppv = sd(nppv, na.rm = TRUE)
      )
}

out <- purrr::map(
  seq_len(nrow(envs)),
  ~ filter_copernicus(cp_data, envs[.x, ]$lat, envs[.x, ]$long)
  ) %>%
  purrr::list_rbind() %>% 
  mutate(
    atoll = envs$atoll, .before = mean_nppv
  )

ggplot(out, aes(atoll, mean_nppv)) +
  geom_line() +
  coord_flip()

res <- data %>% 
  dplyr::filter(
    latitude >= target_lat - 1,
    latitude <= target_lat + 1,
    longitude >= target_long - 1,
    longitude <= target_long + 1
  ) %>% 
  summarise(mean = mean(nppv, na.rm = TRUE))
