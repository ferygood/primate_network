library(dplyr)
library(twice)

hsC1_corr_sig_select <- hsC1_corr_sig %>%
    filter(pair %in% c1_confirm_correlation$x) %>%
    mutate(kznf_age = ifelse(geneName %in% y_kznf, "young", "old")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(link_age = ifelse(kznf_age=="young" | te_age=="young", "young", "old")) %>%
    left_join(hg19rmsk_info, join_by(teName==gene_id))

hsC2_corr_sig_select <- hsC2_corr_sig %>%
    filter(pair %in% c2_confirm_correlation$x) %>%
    mutate(kznf_age = ifelse(geneName %in% y_kznf, "young", "old")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(link_age = ifelse(kznf_age=="young" | te_age=="young", "young", "old")) %>%
    left_join(hg19rmsk_info, join_by(teName==gene_id))

# add kznf age
y_kznf <- kznf_infer %>%
    filter(age=="young")


data("hg19rmsk_info")
hg19rmsk_info[hg19rmsk_info$family_id %in%
                  c("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"),]$family_id <- "SVA"



write.csv(hsC1_corr_sig_select,
          file = "~/Desktop/ferygood_github/primate_network/correlationTables/hmC1_335_sig.csv", row.names = F)

write.csv(hsC2_corr_sig_select,
          file = "~/Desktop/ferygood_github/primate_network/correlationTables/hmC2_64_sig.csv", row.names = F)


cbe_control_specific_select <- cbe_control_specific %>%
    mutate(kznf_age = ifelse(geneName %in% y_kznf, "young", "old")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(link_age = ifelse(kznf_age=="young" | te_age=="young", "young", "old")) %>%
    left_join(hg19rmsk_info, join_by(teName==gene_id))

tcx_control_specific_select <- tcx_control_specific %>%
    mutate(kznf_age = ifelse(geneName %in% y_kznf, "young", "old")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(link_age = ifelse(kznf_age=="young" | te_age=="young", "young", "old")) %>%
    left_join(hg19rmsk_info, join_by(teName==gene_id))

write.csv(cbe_control_specific_select,
          file="~/Desktop/ferygood_github/primate_network/correlationTables/cbe_81_control_sig.csv", row.names=F)

write.csv(tcx_control_specific_select,
          file="~/Desktop/ferygood_github/primate_network/correlationTables/tcx_55_control_sig.csv", row.names=F)
