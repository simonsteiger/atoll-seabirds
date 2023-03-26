library(ncdf4)
library(lubridate)

nc_file <- nc_open("data/test_data.nc")
envs <- read.csv("data/seabird_atolls_envs_10Mar.csv")

time <- nc_file$dim[[1]]$vals
depth <- nc_file$dim[[2]]$vals
latitude <- nc_file$dim[[3]]$vals
longitude <- nc_file$dim[[4]]$vals

x <- ncvar_get(nc_file, nc_file$var[[1]])

matrix <- apply(x, c(1, 2), mean)

colnames(matrix) <- latitude
rownames(matrix) <- longitude

cp_data <- data.frame(longitude = longitude, matrix) %>% 
  tidyr::pivot_longer(-longitude, names_to = "latitude", values_to = "nppv") %>% 
  mutate(
    latitude = (stringr::str_remove(latitude, "X.")),
    latitude = as.numeric(latitude) * -1
  )

target_lat <- envs[6, ]$lat
target_long <- envs[6, ]$long

filter_copernicus <- function(cp_data, lat, long) {
  cp_data %>% 
    dplyr::filter(
      latitude >= lat - 1,
      latitude <= lat + 1,
      longitude >= long - 1,
      longitude <= long + 1
    ) %>% 
    summarise(mean = mean(nppv, na.rm = TRUE))
}

purrr::map(
  c(6, 7),
  ~ filter_copernicus(cp_data, envs[.x, ]$lat, envs[.x, ]$long)
  )

res <- data %>% 
  dplyr::filter(
    latitude >= target_lat - 1,
    latitude <= target_lat + 1,
    longitude >= target_long - 1,
    longitude <= target_long + 1
  ) %>% 
  summarise(mean = mean(nppv, na.rm = TRUE))
