library(scales)

output_dir <- "sim"
dir.create(file.path(output_dir, "sim_fig"), recursive = TRUE)

test.x <- rnorm(6000, 0, 0.7)
df.x <- data.frame(ID = c(1: 6000), test_X = test.x)
write.table(df.x, file = file.path(output_dir, "sim_x.csv"), quote = FALSE, row.names = FALSE ,sep = "\t")

beta_1 <- c(0.075)
beta_2 <- c(0, 0.001)
sigma_1 <- c(0.01, 0.025, 0.05, 0.075, 0.1)
sigma_2 <- c(0.01, 0.015, 0.02, 0.03, 0.04)
pi <- 0.15

df.y <- data.frame(ID = c(1:6000))
counter <- 1
for (idx_beta_2 in 1: length(beta_2)) {
    for (idx_beta_1 in 1:length(beta_1)) {
        for (idx_sigma_2 in 1: length(sigma_2)) {
            for (idx_sigma_1 in 1: length(sigma_1)) {
                name <- paste0("sim", counter)
                test.y <- vector()
                test.y[1: (pi * 6000)] <- 0 + test.x[1: (pi * 6000)] * beta_1[idx_beta_1] + rnorm((pi * 6000), 0, sigma_1[idx_sigma_1])
                test.y[(pi * 6000 + 1): 6000] <- 0 + test.x[(pi * 6000 + 1): 6000] * beta_2[idx_beta_2] + rnorm(((1 - pi) * 6000), 0, sigma_2[idx_sigma_2])
                
                df.y <- cbind(df.y, test.y)
                png(file.path(output_dir, "sim_fig", paste0(name, ".png")), width = 16, height = 12, units = "cm", res = 300)
                plot(test.x[1: (pi * 6000)], test.y[1: (pi * 6000)], main = name, col = alpha("gray", .5), pch = 19, cex = 0.5)
                points(test.x[(pi * 6000 + 1): 6000], test.y[(pi * 6000 + 1): 6000],
                       col = alpha("black", .5), pch = 19, cex = 0.5)
                mtext(paste("beta_1=", beta_1[idx_beta_1], "beta_2=", beta_2[idx_beta_2],
                            "sigma_1=", sigma_1[idx_sigma_1], "sigma_2=", sigma_2[idx_sigma_2], "pi=", pi), side = 3)
                dev.off()
                counter = counter + 1
                colnames(df.y)[counter] <- name
            }
        }
    }
}
write.csv(df.y, file = file.path(output_dir, "sim_y.csv"), quote = FALSE, row.names = FALSE)
