setwd("./iTARGEX/")
del_arr <- read.table("./del.txt", sep = '\t')

# # Test for one pair
# del_arr[,1]
# cor.test(del_arr[,3], del_arr[,2], method = "pearson")
# 
# library("ggpubr")
# ggscatter(del_arr, x = "opi3.del.vs..wt", y = "cho2.del.vs..wt",
#           add = "reg.line", conf.int = TRUE,
#           cor.coef = TRUE, cor.method = "pearson")

# Generate linkage list by pearson correlation
cor_mat <- 1 - cor(as.matrix(del_arr))
linkage <- hclust(as.dist(cor_mat), method = "average")

# Plot dendrogram with pcc and average method clustering
png(filename = "del_arr_dendrogram.png", res = 300)
plot(hclust(as.dist(cor_mat), method = "average"), hang = -1, labels = FALSE)
# abline(h = 0.6, col="red")
dev.off()

# Plot heatmap
library(pheatmap)
breaksList = seq(0, 2, by = 0.2)
col <- colorRampPalette(c("brown", "red", "yellow", "white"))(length(breaksList))

png(filename = "del_arr_heatmap.png", res = 300)
pheatmap(mat = cor_mat[linkage$order, linkage$order], col = col,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 0.8, cellheight = 0.8,
         breaks = breaksList)
dev.off()
