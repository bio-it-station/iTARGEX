
input_dir <- './input'
file_name <- 'phenotypes'
output_dir <- './outliers_detection_figures'

df <- read.csv(file = file.path(input_dir, paste0(file_name, '.csv')))
# png(file.path(output_dir, paste0(file_name, '.png')), width = 20, height = 10, units = "cm", res = 200)
hist(df[,2], breaks = 100, main = file_name, xlab = 'values')
# dev.off()

library(extremevalues)
dist_param = 'lognormal'
y <- getOutliers(df[,2], distribution=dist_param, method = "I")
# png(file.path(output_dir, paste0(file_name, dist_param, '.png')), width = 16, height = 12, units = "cm", res = 200)
outlierPlot(df[,2], y)
# dev.off()
outliers <- c(y$iRight, y$iLeft)
removed_df <- df[-outliers, ]
write.csv(removed_df, file = file.path(input_dir, paste0(file_name, '_rm_outliers.csv')), quote = FALSE, row.names = FALSE)
