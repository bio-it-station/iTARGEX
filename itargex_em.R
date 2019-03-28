#############################################################################################################
#### This program is used to analyze the Correlation of Yeast Gene feature and genetic perturbation data ####
#############################################################################################################

#### load package #####
library(mixtools)
library(Rcpp)
library(data.table)
library(foreach)
library(doParallel)

#### function of mixture regression(model II) ####
regmixEM <- function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2,
                      addintercept = TRUE, epsilon = 1e-08, maxit = 10000, verb = FALSE) {
    s <- sigma
    if (addintercept) {
        x <- cbind(1, x)
    }
    n <- length(y)
    p <- ncol(x)
    tmp <- regmix.init(y = y, x = x, lambda = lambda, beta = beta, s = s, k = k, addintercept = addintercept, arbmean = TRUE, arbvar = TRUE)
    lambda <- tmp$lambda
    beta <- tmp$beta
    s <- tmp$s
    k <- tmp$k
    diff <- 1
    iter <- 0
    xbeta <- x %*% beta
    res <- (y - xbeta) ^ 2
    comp <- t((lambda / sqrt(2 * pi * s ^ 2)) * t(exp(-t(t(res) / (2 * s ^ 2)))))
    obsloglik <- sum(log(apply(comp, 1, sum)))
    ll <- obsloglik
    z <- matrix(nrow = n, ncol = k)
    restarts <- 0
    
    ### EM iterative part ###
    while (diff > epsilon && iter < maxit) {
        ######################  original R code  ##################################
        # for (i in 1: n) {
        #     for (j in 1: k) {
        #         z.denom <- c()
        #         for (h in 1: k) {
        #             z.denom <- c(z.denom, (lambda[h] / lambda[j]) *
        #                              (s[j] / s[h]) *
        #                              exp(-0.5 * ((1 / s[h] ^ 2) * res[i, h] - (1 / s[j] ^ 2) * res[i, j]))
        #                         )
        #         }
        #         z[i, j] <- 1 / sum(z.denom)
        #     }
        # }
        ###########################################################################
        
        ######################  revise C++ code   #################################
        cppFunction('
                    NumericMatrix cpp_for(int n, int k, NumericVector lamda, NumericVector s, NumericMatrix res, NumericMatrix z) {
                    int i, j, h;
                    double shit;
                    for (i = 0; i < n; ++i) {
                    for (j = 0; j < k; ++j) {
                    shit = 0.0;
                    for (h = 0; h < k; ++h) {
                    shit = shit + (lamda[h] / lamda[j]) * (s[j] / s[h]) * exp(-0.5 * ((1 / pow(s[h], 2)) * res(i, h) - (1 / pow(s[j], 2)) * res(i, j)));
                    }
                    z(i, j) = 1 / shit;
                    }
                    }
                    return z;
                    }
                    ')
        
        z <- cpp_for(n, k, lambda, s, res, z)
        ###########################################################################
        z <- z / apply(z, 1, sum)
        lambda.new <- apply(z, 2, mean)
        if (sum(lambda.new < 1e-08) > 0 || is.na(sum(lambda.new))) {
            sing <- 1
        }
        else {
            if (addintercept) {
                lm.out <- lapply(1: k, function(i) lm(y ~ x[, -1], weights = z[, i]))
            }
            else {
                lm.out <- lapply(1: k, function(i) lm(y ~ x - 1, weights = z[, i]))
            }
            beta.new <- sapply(lm.out, coef)
            
            xbeta.new <- x %*% beta.new
            res <- (y - xbeta.new) ^ 2
            s.new <- sqrt(sapply(1: k, function(i) sum(z[, i] * (res[, i]))/sum(z[, i])))
            lambda <- lambda.new
            beta <- beta.new
            xbeta <- x %*% beta
            s <- s.new
            sing <- sum(s < 1e-08)
            comp <- lapply(1: k, function(i) lambda[i] * dnorm(y, xbeta[, i], s[i]))
            comp <- sapply(comp, cbind)
            compsum <- apply(comp, 1, sum)
            newobsloglik <- sum(log(compsum))
        }
        if (sing > 0 || is.na(newobsloglik) || newobsloglik < obsloglik || abs(newobsloglik) == Inf) {
            cat("Need new starting values due to singularity...", "\n")
            restarts <- restarts + 1
            if (restarts > 15) {
                stop("Too many tries!")
            }
            tmp <- regmix.init(y = y, x = x, k = k, addintercept = addintercept, arbmean = TRUE, arbvar = TRUE)
            lambda <- tmp$lambda
            beta <- tmp$beta
            s <- tmp$s
            k <- tmp$k
            diff <- 1
            iter <- 0
            xbeta <- x %*% beta
            res <- (y - xbeta) ^ 2
            comp <- t((lambda / sqrt(2 * pi * s ^ 2)) * t(exp(-t(t(res) / (2 * s ^ 2)))))
            obsloglik <- sum(log(apply(comp, 1, sum)))
            ll <- obsloglik
        }
        else {
            diff <- newobsloglik - obsloglik
            obsloglik <- newobsloglik
            ll <- c(ll, obsloglik)
            iter <- iter + 1
            if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", obsloglik, "\n")
            }
        }
    }
    ### EM iterative part ###
    
    scale.order <- order(s)
    sigma.min <- min(s)
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if (addintercept == FALSE) {
        beta <- rbind(0, beta)
        p = p + 1
    }
    rownames(beta) <- c(paste("beta", 0: (p - 1), sep = ""))
    colnames(beta) <- c(paste("comp", 1: k, sep = ""))
    colnames(z) <- c(paste("comp", 1: k, sep = ""))
    a <- list(x = x,
              y = y,
              lambda = lambda,
              beta = beta,
              sigma = s,
              loglik = obsloglik,
              posterior = z,
              all.loglik = ll,
              restarts = restarts,
              ft = "regmixEM")
    class(a) = "mixEM"
    a
    }

#### function of model selection ####
regmixmodel.sel.wy <- function (x, y, k = 2) {
    p <- function(emout) {
        length(emout$beta) + length(emout$sigma) + length(emout$lambda)
    }
    
    bic <- NULL
    global_param <- c()
    n <- length(y)
    for (i in 1: k) {
        if (i == 1) {
            a <- glm(y ~ x)
            beta <- matrix(a$coef, ncol = 1)
            beta_pval <- summary(a)$coef[2, 4]
            loglik <- sum(log(dnorm(y, mean = beta[1, ] + as.matrix(x) %*% beta[2: nrow(beta), ], sd = sd(a$res))))
            emout <- list(beta = beta, sigma = sd(a$res), lambda = 1, loglik = loglik)
            global_param <- data.frame(cbind(emout$beta[1,], emout$beta[2,], beta_pval, emout$sigma))
            colnames(global_param) <- c("int", "beta", "beta_pval", "sigma")
        }
        else {
            emout <- regmixEM(y, x, k = i)
        }
        bic[i] <- emout$loglik - log(n) * (p(emout) - 1) / 2
    }
    
    bic_df <- t(data.frame(bic, row.names = c("global", "local")))
    local_param <- data.frame(do.call(cbind, Map(cbind, emout$lambda, emout$beta[1,], emout$beta[2,], emout$sigma)))
    colnames(local_param) <- c("lambda_1", "int_1", "beta_1", "sigma_1", "lambda_2", "int_2", "beta_2", "sigma_2")
    weight_comp1 <- t(as.data.frame(emout$posterior[, 1]))
    weight_comp2 <- t(as.data.frame(emout$posterior[, 2]))
    
    out <- list("bic_df" = bic_df,
                "global_param" = global_param, "local_param" = local_param,
                "weight_comp1" = weight_comp1, "weight_comp2" = weight_comp2)
    out
}

#### function for combine the results from parallel EM ####
comb <- function(x, ...) {  
    mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

#### input data ####
args <- commandArgs(trailingOnly = TRUE) # Input 2 arguments. return vector with (input, output)
# args <- c('./input/rep_time.txt', './output/no_intercept') # test input (comment out after test)

del_array <- fread(file = "./input/del.txt", sep = "\t", nThread = 2) # deletion array
# del_array <- read.csv(file = "./sim/sim_y.csv") # simulated deletion array

input <- read.csv(file = args[1], sep = "\t") # trait data
dat.use <- merge(input, del_array, by = "ID") # merge data (deletion array and trait data)
dat.use <- na.omit(dat.use)
dat.use[, 2] <- scale(dat.use[, 2])

#### model selection ####
output_dir <- args[2]
dir.create(file.path(output_dir)) # create the directory for saving raw outputs

cl <- makeCluster(4)
registerDoParallel(cl)

results <- foreach(i = 3: length(dat.use),
                   .packages = c("mixtools", "Rcpp", "data.table"),
                   .combine = 'comb',
                   .multicombine = TRUE) %dopar% {
# for(i in 3: length(dat.use)){
    # set.seed(10)
    em <- regmixmodel.sel.wy(x = dat.use[, 2], y = dat.use[, i], k = 2) # EM and Model selection
    # if (abs(model_select$EM_param[3, 1]) < abs(model_select$EM_param[3, 2])) { # Switch comp1 with larger slope one
    #     model_select$EM_param <- model_select$EM_param[c(2, 1)]
    #     colnames(model_select$EM_param) <- c("comp1", "comp2")
    # }
    
    # Output files
    rownames(em$bic_df) <- colnames(dat.use)[i]
    rownames(em$global_param) <- colnames(dat.use)[i]
    rownames(em$local_param) <- colnames(dat.use)[i]
    rownames(em$weight_comp1) <- colnames(dat.use)[i]
    rownames(em$weight_comp2) <- colnames(dat.use)[i]
    
    result <- list("bic_df" = em$bic_df,
                   "global_param" = em$global_param,
                   "local_param" = em$local_param,
                   "weight_comp1" = em$weight_comp1,
                   "weight_comp2" = em$weight_comp2)
    
    # filename <- file.path(output_dir, "raw_output", colnames(dat.use)[i])
    # write.csv(model_select$model_select, file = paste0(filename, "_model_select.csv"), quote = F) # model selection for each mutant
    # write.csv(model_select$EM_param, file = paste0(filename, "_param.csv"), quote = F) # model II parameters
    # fwrite(model_select$weight, file = paste0(filename, "_weight.csv"), quote = F, nThread = 1) # posterior probabilities from model II
    return(result)
}

stopCluster(cl)

# Apply Bonferroni correction to global beta p-value
global_param <- results[["global_param"]]
pval_adj_df <- p.adjust(global_param[, 3], method = "bonferroni")
global_param['adj_pval'] <- pval_adj_df
global_param <- global_param[order(global_param$beta_pval), ]

write.csv(results[["bic_df"]], file = file.path(output_dir, "bic_df.csv"), quote = F) # model selection for each mutant
write.csv(global_param, file = file.path(output_dir, "global_param.csv"), quote = F) # model I parameters
write.csv(results[["local_param"]], file = file.path(output_dir, "local_param.csv"), quote = F) # model II parameters
fwrite(results[["weight_comp1"]], file = file.path(output_dir, "weight_comp1.csv"), quote = F, nThread = 2) # posterior probabilities from model II
fwrite(results[["weight_comp2"]], file = file.path(output_dir, "weight_comp2.csv"), quote = F, nThread = 2) # posterior probabilities from model II
