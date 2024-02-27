.libPaths("/data/scratch2/yaochung41/RLib")
source("functions.R")

library(TEKRABber)
library(dplyr)

load("primateBrainData.RData")

# genes
df_hm_gene <- hmGene[,c(-1)]
rownames(df_hm_gene) <- hmGene$geneID

# transposable elements
hsTEexp <- hmTE %>% select(-c(1,2,3))
rownames(hsTEexp) <- hmTE$name  #908 TEs

# genes convert to tpm
sample_counts <- colSums(df_hm_gene)

scaling_factor <- sample_counts / 1e6

df_hm_gene_tpm <- df_hm_gene
df_hm_gene_tpm <- t(t(df_hm_gene_tpm)/ scaling_factor + 1) * 1e6
df_hm_gene_tpm <- as.data.frame(df_hm_gene_tpm)

# tes convert to tpm
te_count <- colSums(hsTEexp)
te_scale <- te_count / 1e6
hsTE_tpm <- hsTEexp
hsTE_tpm <- t(t(hsTE_tpm)/ te_scale + 1) * 1e6
hsTE_tpm <- as.data.frame(hsTE_tpm)

cluster_meta <- read.csv("data/cluster_meta.csv")

set.seed(47)
gene_list <- rownames(df_hm_gene_tpm)
selected_genes <- replicate(1000,
                            sample(gene_list, size=337, replace=FALSE),
                            simplify = FALSE)



cluster_id_c1 <- cluster_meta %>%
    filter(cluster == "cluster1") %>%
    select(1)

cluster_id_c2 <- cluster_meta %>%
    filter(cluster == "cluster2") %>%
    select(1)

for (i in 1:1000){

    gene_set <- selected_genes[[i]]
    gene <- df_hm_gene_tpm
    te <- hsTE_tpm

    cluster_gene_c1 <- gene %>% select(cluster_id_c1$Run)
    cluster_gene_c2 <- gene %>% select(cluster_id_c2$Run)

    cluster_gene_c1 <- cluster_gene_c1[rownames(cluster_gene_c1) %in% gene_set, ]
    cluster_gene_c2 <- cluster_gene_c2[rownames(cluster_gene_c2) %in% gene_set, ]

    cluster_te_c1 <- te %>% select(cluster_id_c1$Run)
    cluster_te_c2 <- te %>% select(cluster_id_c2$Run)

    df_c1 <- corrOrthologTE(
        geneInput = cluster_gene_c1,
        teInput = cluster_te_c1,
        numCore = 5,
        fileDir = "./results_c1",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )

    df_c2 <- corrOrthologTE(
        geneInput = cluster_gene_c2,
        teInput = cluster_te_c2,
        numCore = 5,
        fileDir = "./results_c2",
        fileName = paste0("gene_", i, "_vs_TE_corr.csv")
    )

}

