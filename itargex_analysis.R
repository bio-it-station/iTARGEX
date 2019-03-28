################################
#### Analysis of EM results ####
################################

library(data.table)
library(foreach)
library(doParallel)
library(scales)

args <- commandArgs(trailingOnly = TRUE) # Input 2 arguments. return vector with (input, output)
# args <- c('./input/rep_time.txt', './output/test') # test input (comment out after test)

del_array <- fread(file = "./input/del.txt", sep = "\t", nThread = 2) # deletion array
# del_array <- read.csv(file = "./sim/sim_y.csv") # simulated deletion array

input <- read.csv(file = args[1], sep = "\t") # trait data
dat.use <- merge(input, del_array, by = "ID") # merge data (deletion array and trait data)
dat.use <- na.omit(dat.use)
dat.use[, 2] <- scale(dat.use[, 2])
x <- dat.use[, 2]
y <- dat.use[, -c(1, 2)]
output_dir <- args[2]

param <- read.csv(file = file.path(output_dir, "local_param.csv"), row.names = 1)
weight_c1 <- fread(file = file.path(output_dir, "weight_comp1.csv"), nThread = 2)
weight_c2 <- fread(file = file.path(output_dir, "weight_comp2.csv"), nThread = 2)

#### correlation computation ####
# p-value computation based on t-test
t_test_pval <- function (n, r) {
	pval_tail_1 <- pt(r * sqrt(n - 2) / sqrt(1 - r ^ 2), n - 2)
    pval_tail_2 <- 1 - pval_tail_1
    2 * min(pval_tail_1, pval_tail_2)
}

test_beta <- function(x, y) {
    a <- lm(y ~ x)
    summary(a)$coefficients[2, c(1, 4)]
}

#### function for combine the results from parallel processes ####
comb <- function(x, ...) {  
    mapply(rbind, x, ..., SIMPLIFY=FALSE)
}

cl <- makeCluster(4)
registerDoParallel(cl)

cor_df <- foreach(i = 1: ncol(y), .combine = 'comb', .multicombine = TRUE) %dopar% {
# for (i in 1: ncol(y)) {
    
    # # Soft assign
    # cor_comp1 <- cor(x * weight_c1[i, ], y[, i] * weight_c1[i, ]) # correlation of component_1 computation
    # cor_comp2 <- cor(x * weight_c2[i, ], y[, i] * weight_c2[i, ]) # correlation of component_2 computation
    # pval_cor_comp1 <- t_test_pval(n = sample_num, r = cor_comp1)
    # pval_cor_comp2 <- t_test_pval(n = sample_num, r = cor_comp2)
    
    # Hard assign
    comp1_idx <- which(weight_c1[i, ] > 0.5)
    comp2_idx <- which(weight_c2[i, ] > 0.5)
    cor_comp1 <- cor(x[comp1_idx], y[, i][comp1_idx])
    cor_comp2 <- cor(x[comp2_idx], y[, i][comp2_idx])
    pval_cor_comp1 <- t_test_pval(n = length(comp1_idx), r = cor_comp1)
    pval_cor_comp2 <- t_test_pval(n = length(comp2_idx), r = cor_comp2)
    pval_beta_comp1 <- test_beta(x[comp1_idx], y[, i][comp1_idx])
    pval_beta_comp2 <- test_beta(x[comp2_idx], y[, i][comp2_idx])
    cor_local <- c(cor_comp1, pval_cor_comp1, pval_beta_comp1, cor_comp2, pval_cor_comp2, pval_beta_comp2)
    
    return(cor_local)
    
    # # Combine all data
    # cor_df <- rbind(cor_df, cor_local)
}

stopCluster(cl)

cor_df <- as.data.frame(cor_df, row.names = colnames(y)[1: ncol(y)])
colnames(cor_df) <- c("cor_c1", "pval_cor_c1", "beta_c1", "pval_beta_c1",
                      "cor_c2", "pval_cor_c2", "beta_c2", "pval_beta_c2")

# Apply Bonferroni correction to p-value
pval_adj_df <- sapply(cor_df[, c(2, 6)], p.adjust, method = "bonferroni")
colnames(pval_adj_df) <- c("adjp_c1", "adjp_c2")
cor_df <- cbind(cor_df, pval_adj_df)

# -log10(p-value) transformation
logp_df <- sapply(cor_df[, c(9, 10)], log10) * (-1) # apply -log10 transformation to dataframe
logp_df[is.infinite(logp_df)] <- 350 # Set the p-value with infinite value to 350
cor_df[, c(9, 10)] <- logp_df

# Concatenate more information to the dataframe 
cor_df <- cbind(cor_df, param$beta_1, param$beta_2, param$lambda_1, param$lambda_2)
colnames(cor_df)[c(11: 14)] <- c("em_beta_1", "em_beta_2", "ratio_1", "ratio_2")

# Write out the dataframe with all data
fwrite(format(cor_df, digits = 6), file = file.path(output_dir, "correlation.csv"), quote = F, row.names = T, nThread = 2)

# # GO annotation
# clus = read.csv('./1-s2.0-S0092867414003420-mmc1.csv', header = T, row.names = 1, stringsAsFactors = F)[2]
# cor_df <- merge(cor_df, clus, by = "row.names")
# colnames(cor_df)[1] <- c("gene")

# Select and output the significant cases
cor_sig_c1 <- subset(cor_df, adjp_c1 > -log10(0.01) & !(adjp_c2 > -log10(0.01)),
                     select = c(1: 4, 9, 11, 13))
cor_sig_c2 <- subset(cor_df, adjp_c2 > -log10(0.01) & !(adjp_c1 > -log10(0.01)),
                     select = c(5: 8, 10, 12, 14))
cor_sig_both <- subset(cor_df, adjp_c1 > -log10(0.01) & adjp_c2 > -log10(0.01))

cor_sig_local <- rbind(cor_sig_c1, setNames(cor_sig_c2, names(cor_sig_c1)))
colnames(cor_sig_local) <- c("cor", "pval_beta", "pval.adj", "em_beta", "ratio")

cor_sig_local <- cor_sig_local[order(cor_sig_local$pval_beta), ]

write.csv(format(cor_sig_both, digits = 6), file = file.path(output_dir, "cor_sig_both.csv"))
write.csv(format(cor_sig_local, digits = 6), file = file.path(output_dir, "cor_sig_local.csv"))

# Plot xy_scatter for significant cases
dir.create(file.path(output_dir, "figure_local"))
dir.create(file.path(output_dir, "figure_both"))
xy_plot <- function(sig_df, param. = param, weight_c1. = weight_c1, weight_c2. = weight_c2,
                    major = c("c1", "c2", "both")) {
    idx <- which(row.names(param.) %in% row.names(sig_df))
    for (i in 1: length(idx)) {
        # Set Output folder and select points for each distribution
        filename <- file.path(output_dir, "figure_local")
        if (major == "c1") {
            sig_points <- which(weight_c1.[idx[i], ] > 0.5)
            sig_param <- c(2, 3)
            nonsig_points <- which(weight_c2.[idx[i], ] > 0.5)
            nonsig_param <- c(6, 7)
        }
        else if (major == "c2"){
            sig_points <- which(weight_c2.[idx[i], ] > 0.5)
            sig_param <- c(6, 7)
            nonsig_points <- which(weight_c1.[idx[i], ] > 0.5)
            nonsig_param <- c(2, 3)
        }
        else {
            filename <- file.path(output_dir, "figure_both")
            if (length(which(weight_c1.[idx[i], ] > 0.5)) > length(which(weight_c2.[idx[i], ] > 0.5))) {
                sig_points <- which(weight_c1.[idx[i], ] > 0.5)
                sig_param <- c(2, 3)
                nonsig_points <- which(weight_c2.[idx[i], ] > 0.5)
                nonsig_param <- c(6, 7)
            }
            else {
                sig_points <- which(weight_c2.[idx[i], ] > 0.5)
                sig_param <- c(6, 7)
                nonsig_points <- which(weight_c1.[idx[i], ] > 0.5)
                nonsig_param <- c(2, 3)
            }
        }
        
        png(file.path(filename, paste0(row.names(param.)[idx[i]], ".png")),
            width = 16, height = 12, units = "cm", res = 300)
        # Nonsig_plot
        plot(x[nonsig_points], y[, idx[i]][nonsig_points],
             xlim = c(min(x), max(x)),
             ylim = c(min(y[, idx[i]]), max(y[, idx[i]])),
             main = row.names(param.)[idx[i]],
             xlab = "Traits", ylab = "log(Fold-change)",
             col = alpha("#5858FA", .5), pch = 19, cex = 0.5)
        
        # Sig_plot
        points(x[sig_points], y[, idx[i]][sig_points],
               col = alpha("#FA5858", .5), pch = 19, cex = 0.5)
        
        # Regression line
        abline(coef = param.[idx[i], nonsig_param], col = "#0000FF")
        # abline(lm(y[, idx[i]][nonsig_points] ~ x[nonsig_points]), col = "#0000FF")
        abline(coef = param.[idx[i], sig_param], col = "#FF0000")
        # abline(lm(y[, idx[i]][sig_points] ~ x[sig_points]), col = "#FF0000")
        dev.off()
    }
}

xy_plot(cor_sig_c1, major = "c1")
xy_plot(cor_sig_c2, major = "c2")
xy_plot(cor_sig_both, major = "both")

# #### manhatten like plot ####
# cor.global=read.delim("./cor_global_only.txt",header = T,stringsAsFactors = F)
# cor.local=read.delim("./cor_local_only.txt",header = T,stringsAsFactors = F)
# clus = read.csv(gene_description, header = T,stringsAsFactors = F) #functional clustering and functional description of the mutants
# atr=clus[clus$gene%in%union(cor.global[,1],cor.local[,1]),]
# dat.global=merge(cor.global,atr[c(1,2,4)],by= "gene")
# dat.global$logP=-log10(dat.global$pvalue)
# dat.global$logP[dat.global$logP==Inf]=350
# dat.local=merge(cor.local,atr[c(1,2,4)],by = "gene")
# dat.local$logP=-log10(dat.local$pvalue)
# dat.local$logP[dat.local$logP==Inf]=350
# png(paste0("./manhatten_like_plot_",colnames(dat.use)[2],".png"),width = 20,height = 12,units = "cm",pointsize = 10,res = 600)
# par(mar=c(14,4,4,2))
# plot(1,1,xlim = c(1,16),ylim = c(0,max(c(max(dat.global$logP),max(dat.local$logP)))+10),type = 'n',xlab = "",ylab = "-log10(p)",axes = F)
# rect(seq(0.5,15.5,2),rep(0,8),seq(1.5,16.5,2),rep(400,8),col = "grey70",border = NA)
# for(i in 1:length(unique(atr$functional.category))){
#     index1=dat.global$functional.category%in%unique(atr$functional.category)[i]
#     index2=dat.local$functional.category%in%unique(atr$functional.category)[i]
#     points(rnorm(i-0.15,0.05,n = sum(index1)),dat.global$logP[index1],col=ifelse(dat.global$logP[index1]<summary(dat.global$logP)[5],"black","orange"),cex=0.7,pch=19)
#     points(rnorm(i+0.15,0.05,n = sum(index2)),dat.local$logP[index2],col=ifelse(dat.local$logP[index2]<summary(dat.local$logP)[5],"black","green"),cex=0.7,pch=19)
# }
# axis(1,at = 1:16,labels = unique(atr$functional.category),las=2)
# axis(2)
# abline(h = summary(dat.global$logP)[5],col="red")
# abline(h = summary(dat.local$logP)[5],col="blue")
# dev.off()
