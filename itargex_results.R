#############################
#### Data Rearrangement  ####
#############################

subdirs <- dir(path = "./output", pattern = "test_*", full.names = TRUE)
num_of_replicate = length(subdirs)

df <- data.frame()
for (i in 1: num_of_replicate) {
    next_df <- read.csv(file.path(subdirs[i], "cor_sig_local.csv"))
    df = rbind(df, next_df)
}
colnames(df)[1] <- c("gene")

observed_genes <- unique(df$gene)

mean_df <- c()
sd_df <- c()
for (i in 1: length(observed_genes)) {
    current_df <- subset(df, df$gene == observed_genes[i], select = c(-1))
    observed_prob <- nrow(current_df) / num_of_replicate
    names(observed_prob) <- c("prob")
    
    new_df_mean <- apply(current_df, 2, mean)
    new_df_sd <- apply(current_df, 2, sd)
    
    mean_df <- rbind(mean_df, c(observed_prob, new_df_mean))
    sd_df <- rbind(sd_df, new_df_sd)
}
mean_df <- as.data.frame(mean_df)
rownames(mean_df) <- observed_genes
sd_df <- as.data.frame(sd_df)
rownames(sd_df) <- observed_genes
pval_order <- order(mean_df$pval_beta)
mean_df <- mean_df[pval_order, ]
sd_df <- sd_df[pval_order, ]

write.csv(mean_df, file = "./output/10_times_mean.csv", quote = FALSE)
write.csv(sd_df, file = "./output/10_times_sd.csv", quote = FALSE)
