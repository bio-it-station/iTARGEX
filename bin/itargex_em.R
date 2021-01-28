#!/usr/bin/env Rscript
'
iTARGEX EM Algorithm.

Usage:
  itargex_em.R [options] <input-trait> <deletome-data> <output-folder>
  itargex_em.R -h | --help

Options:
  -h --help        Show this screen.
  --version        Show version.
  --iter=<int>     Max iteration number of EM algorithm [default: 1000].
  --ncpus=<int>    Number of CPU for parallel execuation [default: 4].
' -> doc

library(docopt)
args <- docopt(doc)
tryCatch(
    {
        args$iter <- as.numeric(args$iter)
        args$ncpus <- as.integer(args$ncpus)
    },
    warning = function(c) {
        message("Error: option should be numeric value.")
        quit(save = "no", status = 1)
    }
)

if (args$ncpus < 1) {
    message("Error: --ncpus should larger than 1.")
    quit(save = "no", status = 1)
} else if (args$iter < 500) {
    message("Error: --iter should larger than 500.")
    quit(save = "no", status = 1)
}

#### load package #####
library(mixtools)
library(Rcpp)
library(data.table)
library(foreach)
library(doParallel)

#### function of mixture regression(model II) ####
regmixEM <- function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2,
                      addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, maxit = 1000, verb = FALSE) {
    if (arbmean == FALSE && arbvar == FALSE) {
        stop(paste("Must change constraints on beta and/or sigma!","\n"))
    }
    
    #### TCH save initial value
        beta.init <- beta
        lambda.init <- lambda
    ####
    
    s <- sigma
    if (addintercept) {
        x <- cbind(1, x)
    }
    n <- length(y)
    p <- ncol(x)
    tmp <- regmix.init(y = y, x = x, lambda = lambda, beta = beta, s = s, k = k,
                       addintercept = addintercept, arbmean = arbmean, arbvar = arbvar)
    lambda <- tmp$lambda
    beta <- tmp$beta
    s <- tmp$s
    k <- tmp$k
    diff <- 1
    iter <- 0
    xbeta <- x %*% beta
    res <- (y - xbeta) ^ 2
    if (arbmean == FALSE) {
        res <- sapply(1: k, function(i) res)
    }
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
        cppFunction(
            '
            NumericMatrix cpp_for(int n, int k, NumericVector lambda, NumericVector s, NumericMatrix res, NumericMatrix z, int arbvar) {
            int i, j, h;
            double shit;
            for (i = 0; i < n; ++i) {
                for (j = 0; j < k; ++j) {
                    shit = 0.0;
                    for (h = 0; h < k; ++h) {
                        shit += (lambda[h] / lambda[j]) * (s[j * arbvar] / s[h * arbvar]) * exp(-0.5 * ((1 / pow(s[h * arbvar], 2)) * res(i, h) - (1 / pow(s[j * arbvar], 2)) * res(i, j)));
                    }
                    z(i, j) = 1 / shit;
                }
            }
            return z;
            }
            '
        )
        
        z <- cpp_for(n, k, lambda, s, res, z, arbvar)
        ###########################################################################
        z <- z / apply(z, 1, sum)
        lambda.new <- apply(z, 2, mean)
        
        if (sum(lambda.new < 1e-08) > 0 || is.na(sum(lambda.new))) {
            sing <- 1
        }
        else {
            if (arbmean == FALSE) {
                if (addintercept) {
                    beta.new <- lm(y ~ x[, -1], weights = apply(t(t(z) / (s ^ 2)), 1, sum))$coef
                }
                else {
                    beta.new <- lm(y ~ x - 1, weights = apply(t(t(z) / (s ^ 2)), 1, sum))$coef
                }
            }
            else {
                if (addintercept) {
                    lm.out <- lapply(1: k, function(i) lm(y ~ x[, -1], weights = z[, i]))
                }
                else {
                    lm.out <- lapply(1: k, function(i) lm(y ~ x - 1, weights = z[, i]))
                }
                
                beta.new <- sapply(lm.out, coef)
                #### wylai set beta_2
                # beta.new[, 2] <- c(0, 0)
                ####
                
            }
            xbeta.new <- x %*% beta.new
            res <- (y - xbeta.new) ^ 2
            if (arbmean == FALSE){
                res <- sapply(1: k, function(i) res)
            }
            if (arbvar) {
                s.new <- sqrt(sapply(1: k, function(i) sum(z[, i] * (res[, i])) / sum(z[, i])))
            }
            else {
                s.new <- sqrt(sum(z * res) / n)
            }
            lambda <- lambda.new
            beta <- beta.new
            xbeta <- x %*% beta
            s <- s.new
            sing <- sum(s < 1e-08)
            comp <- lapply(1: k, function(i) lambda[i] * dnorm(y, xbeta[, i * arbmean + (1 - arbmean)], s[i * arbvar + (1 - arbvar)]))
            comp <- sapply(comp, cbind)
            compsum <- apply(comp, 1, sum)
            newobsloglik <- sum(log(compsum))
        }
        if (sing > 0 || is.na(newobsloglik) || newobsloglik < obsloglik || abs(newobsloglik) == Inf) {
            cat("Need new starting values due to singularity...", "\n")
            restarts <- restarts + 1
            if (restarts > 5) {
                stop("Too many tries!")
            }
            
            #### TCH Using previosly saved initial value to start a new round
            tmp <- regmix.init(y = y, x = x, k = k, lambda = lambda, beta = beta, addintercept = addintercept, arbmean = arbmean, arbvar = arbvar)
            ####
            
            lambda <- tmp$lambda
            beta <- tmp$beta
            s <- tmp$s
            k <- tmp$k
            diff <- 1
            iter <- 0
            xbeta <- x %*% beta
            res <- (y - xbeta) ^ 2
            if (arbmean == FALSE){
                res <- sapply(1: k, function(i) res)
            }
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
    conv <- "Y"
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
        conv <- "N"
    }
    cat("number of iterations=", iter, "\n")
    if (addintercept == FALSE) {
        beta <- rbind(0, beta)
        p = p + 1
    }
    if (arbmean == FALSE) {
        z <- z[, scale.order]
        names(beta) <- c(paste0("beta", 0: (p - 1)))
        colnames(z) <- c(paste0("comp", 1: k))
        a <- list(x = x,
                  y = y,
                  lambda = lambda[scale.order],
                  beta = beta,
                  sigma = sigma.min,
                  scale = s[scale.order] / sigma.min,
                  loglik = obsloglik,
                  posterior = z[, scale.order],
                  all.loglik = ll,
                  restarts = restarts,
                  conv = conv,
                  iter = iter,
                  ft = "regmixEM")
    }
    else {
        rownames(beta) <- c(paste0("beta", 0: (p - 1)))
        colnames(beta) <- c(paste0("comp", 1: k))
        colnames(z) <- c(paste0("comp", 1: k))
        a <- list(x = x,
                  y = y,
                  lambda = lambda,
                  beta = beta,
                  sigma = s,
                  loglik = obsloglik,
                  posterior = z,
                  all.loglik = ll,
                  restarts = restarts,
                  conv = conv,
                  iter = iter,
                  ft = "regmixEM")
    }
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
            r <- cor(x, y)
            global_param <- data.frame(cbind(emout$beta[1,], emout$beta[2,], beta_pval, emout$sigma, r))
            colnames(global_param) <- c("int", "beta", "pval", "sigma", "cor")
        }
        else {
            beta <- matrix(as.numeric(global_param[c(1, 2)]), nrow = 2, ncol = 1)
            beta <- cbind(beta, c(0, 0))
            emout <- regmixEM(y, x, beta = beta, k = i, maxit = args$iter)
        }
        bic[i] <- emout$loglik - log(n) * (p(emout) - 1) / 2
    }
    
    bic_df <- t(data.frame(bic, row.names = c("global", "local")))
    local_param <- data.frame(do.call(cbind, Map(cbind, emout$lambda, emout$beta[1,], emout$beta[2,], emout$sigma)))
    local_param <- cbind(local_param, emout$conv, emout$iter)
    colnames(local_param) <- c("lambda_1", "int_1", "beta_1", "sigma_1", "lambda_2", "int_2", "beta_2", "sigma_2", "conv", "iter")
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

input <- read.csv(file = args$input_trait) # trait data
del_array <- fread(file = args$deletome_data, nThread = args$ncpus)
dat.use <- merge(input, del_array, by = "ID") # merge data (deletion array and trait data)
dat.use <- na.omit(dat.use)
dat.use[, 2] <- scale(dat.use[, 2])

#### model selection ####
output_dir <- args$output_folder
dir.create(file.path(output_dir)) # create the directory for saving raw outputs

cl <- makeCluster(args$ncpus)
registerDoParallel(cl)

results <- foreach(i = 3: length(dat.use),
                   .packages = c("mixtools", "Rcpp", "data.table"),
                   .combine = 'comb',
                   .multicombine = TRUE) %dopar% {
# for(i in 3: length(dat.use)){
    # set.seed(10)
    em <- regmixmodel.sel.wy(x = dat.use[, 2], y = dat.use[, i], k = 2) # EM and Model selection
    
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
    
    return(result)
}

stopCluster(cl)

# Apply Bonferroni correction to global beta p-value
global_param <- results[["global_param"]]
pval_adj_df <- p.adjust(global_param[, "pval"], method = "bonferroni")
global_param['adj_pval'] <- pval_adj_df
logp_df <- sapply(global_param[, 'adj_pval'], log10) * (-1) # apply -log10 transformation to dataframe
logp_df[is.infinite(logp_df)] <- 350 # Set the p-value with infinite value to 350
global_param[, 'adj_pval'] <- logp_df

write.csv(results[["bic_df"]], file = file.path(output_dir, "bic_df.csv"), quote = F) # model selection for each mutant
write.csv(global_param, file = file.path(output_dir, "global_param.csv"), quote = F) # model I parameters
write.csv(results[["local_param"]], file = file.path(output_dir, "local_param.csv"), quote = F) # model II parameters
fwrite(results[["weight_comp1"]], file = file.path(output_dir, "weight_comp1.csv"), quote = F, nThread = args$ncpus) # posterior probabilities from model II
fwrite(results[["weight_comp2"]], file = file.path(output_dir, "weight_comp2.csv"), quote = F, nThread = args$ncpus) # posterior probabilities from model II
