library(scales)

sim_dir = './simulation_190815'

test.x <- rnorm(6000, 0, 0.7)
df.x <- data.frame(ID = c(1: 6000), test_X = test.x)
write.csv(df.x, file = file.path(sim_dir, "sim_x.csv"), quote = FALSE, row.names = FALSE)
# max_x <- ceiling(max(test.x)) + 1
# min_x <- floor(min(test.x)) - 1

for (test in 1: 100) {
    output_dir <- paste0("sim_", test)
    # dir.create(file.path(output_dir, "sim_fig"), recursive = TRUE)
    dir.create(file.path(sim_dir, output_dir), recursive = TRUE)
    ############################################
    # Parameter space
    beta_1 <- c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1) # Slope of component 1
    sigma_1 <- 0.44 # Sigma of component 1
    beta_2 <- 0 # Slope of component 2
    sigma_2 <- 0.11 # Sigma of component 2
    pi <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0) # Ratio between component 1 and component 2
    ############################################
    
    df.y <- data.frame(ID = c(1: 6000))
    counter <- 1
    for (idx_pi in 1: length(pi)) {
        for (idx_beta_2 in 1: length(beta_2)) {
            for (idx_beta_1 in 1: length(beta_1)) {
                for (idx_sigma_2 in 1: length(sigma_2)) {
                    for (idx_sigma_1 in 1: length(sigma_1)) {
                        # name <- paste0("sim", counter)
                        name <- paste(pi[idx_pi], beta_1[idx_beta_1], sep = "_")
                        test.y <- vector()
                        comp1_num <- pi[idx_pi] * 6000
                        comp2_num <- 6000 - comp1_num
                        if (comp1_num != 0) {
                            test.y[1: comp1_num] <- test.x[1: comp1_num] * beta_1[idx_beta_1] + rnorm(comp1_num, 0, sigma_1[idx_sigma_1])
                        }
                        if (comp2_num != 0) {
                            test.y[(comp1_num + 1): 6000] <- test.x[(comp1_num + 1): 6000] * beta_2[idx_beta_2] + rnorm(comp2_num, 0, sigma_2[idx_sigma_2])
                        }
                        df.y <- cbind(df.y, test.y)
                        # png(file.path(output_dir, "sim_fig", paste0(name, ".png")), width = 16, height = 12, units = "cm", res = 200)
                        # plot(test.x[1: (pi * 6000)], test.y[1: (pi * 6000)], main = name, col = alpha("gray", .5), pch = 19, cex = 0.5, xlim = c(min_x, max_x), ylim = c(-0.5, 0.5))
                        # points(test.x[(pi * 6000 + 1): 6000], test.y[(pi * 6000 + 1): 6000],
                        #        col = alpha("black", .5), pch = 19, cex = 0.5)
                        # mtext(paste("beta_1=", beta_1[idx_beta_1], "beta_2=", beta_2[idx_beta_2],
                        #             "sigma_1=", sigma_1[idx_sigma_1], "sigma_2=", sigma_2[idx_sigma_2], "pi=", pi), side = 3)
                        # dev.off()
                        counter = counter + 1
                        colnames(df.y)[counter] <- name
                    }
                }
            }
        }
    }
    
    write.csv(df.y, file = file.path(sim_dir, output_dir, "sim_y.csv"), quote = FALSE, row.names = FALSE)
}
