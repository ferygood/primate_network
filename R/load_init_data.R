# 1. primate brain result
read_rds_files("data/results_rdata")

# 2. mayo RNA-seq
load("data/mayoTEKRABber_balance.RData")

# 3. age inference
kznf_infer <- read.csv("data/kznf_bucket.csv")
te_infer <- read.csv("data/Dfam_TE_simiiformes.csv", row.names = 1)

detach(package:dplyr)
library(dplyr)
conflicted::conflict_prefer("filter", "dplyr")