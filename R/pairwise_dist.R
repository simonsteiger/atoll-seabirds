

library("tidyverse")
library("geosphere")

dat <- read.csv("C:/Users/sebas/OneDrive/Dokumente/Postdoc/manuscripts/atoll_seabird_conservation/seabird_atolls_envs_02May.csv")


pairwdist <- distm(dat[c("long","lat")]) / 1000 # default is in metres, divide by 1000 for km distance
rownames(pairwdist) <- dat$atoll
colnames(pairwdist) <- dat$atoll

pairwdist <- round(pairwdist)

write.csv(pairwdist, "C:/Users/sebas/OneDrive/Dokumente/Postdoc/manuscripts/atoll_seabird_conservation/pairwise_atoll_distance.csv")
