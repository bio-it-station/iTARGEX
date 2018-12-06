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
png(filename = "./output/del_arr_dendrogram.png", width = 30, height = 15, units = "cm", res = 300)
plot(hclust(as.dist(cor_mat), method = "average"), hang = -1, labels = FALSE)
# abline(h = 0.6, col="red")
dev.off()

# Plot heatmap
library(pheatmap)
breaksList = seq(0, 2, by = 0.2)
col <- colorRampPalette(c("red", "white", "white", "blue"))(length(breaksList) - 1)

png(filename = "./output/del_arr_heatmap.png", width = 25, height = 25, units = "cm", res = 300)
pheatmap(mat = cor_mat[linkage$order, linkage$order], col = col,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 0.9, cellheight = 0.9,
         breaks = breaksList,
         legend_breaks = c(0, 1, 2), legend_labels = c("Pos cor", "no cor", "Neg cor"))
dev.off()
