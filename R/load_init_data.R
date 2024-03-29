# 1. primate brain result
# read_rds_files("data/results_rdata")
load("data/primateBrainData.RData")

# 2. mayo RNA-seq
load("data/mayoTEKRABber_balance.RData")

# 3. age inference
kznf_infer <- utils::read.csv("tables/kznf_bucket.csv")
te_infer <- utils::read.csv("tables/Dfam_TE_simiiformes.csv", row.names = 1)

detach("package:dplyr", unload=TRUE)
library(dplyr)
conflicted::conflict_prefer("filter", "dplyr")
