n <- 500
x <- rnorm(n, sd = 4)
y <- -0.1 + 0.01 * x + rnorm(n, sd = 0.05)
plot(x, y)
a <- lm(y ~ x)
beta <- matrix(a$coef, ncol = 1)
loglik <- sum(log(dnorm(y, mean = beta[1, ] + as.matrix(x) %*% beta[2: nrow(beta), ], sd = sd(a$res))))
emout <- list(beta = beta, sigma = sd(a$res), lambda = 1, loglik = loglik)
summary(a)$coefficients[,4]
BIC(a)

bic_log = c(BIC(a))
lowest_bic = BIC(a)
best_model <- emout
# for (i in 2:5) {
emout1 <- regmixEM(y , x, k = 2)
current_bic <- -2 * emout1$loglik + log(n) * (p(emout)) / 2
    # if (current_bic < lowest_bic) {
    #     lowest_bic <- current_bic
    #     best_model <- emout1
    # }
    # bic_log <- c(bic_log, current_bic)
# }
plot(bic_log, type = 'b')
