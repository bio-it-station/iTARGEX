#############################
#### Data Rearrangement  ####
#############################

args <- commandArgs(trailingOnly = TRUE)
subdirs <- dir(path = args[1], pattern = "output_*", full.names = TRUE)
num_of_replicate = length(subdirs)

df <- data.frame()
for (i in 1: num_of_replicate) {
    next_df <- read.csv(file.path(subdirs[i], "cor_sig_local.csv"))
    df = rbind(df, next_df)
}
colnames(df)[1] <- c("gene")

observed_genes <- unique(df$gene)

summary_df <- c()
for (i in 1: length(observed_genes)) {
    current_df <- subset(df, df$gene == observed_genes[i], select = c(-1))
    observed_prob <- nrow(current_df) / num_of_replicate
    names(observed_prob) <- c("prob")
    
    new_df_mean <- apply(current_df, 2, mean)
    new_df_sd <- apply(current_df, 2, sd)
    
    new_df_comb <- do.call(cbind, Map(cbind, new_df_mean, new_df_sd))
    names(new_df_comb) <- do.call(c, Map(cbind, names(new_df_mean), paste0(names(new_df_sd), "_sd")))
    new_df_comb <- c(observed_prob, new_df_comb)
    summary_df <- rbind(summary_df, new_df_comb)
}
summary_df <- as.data.frame(summary_df)
rownames(summary_df) <- observed_genes
summary_df <- summary_df[order(summary_df$pval_beta), ]

write.csv(format(summary_df, digits = 4), file = file.path(args[1], "summary.csv"), quote = FALSE)
