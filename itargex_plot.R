library(data.table)

args <- c('./input/rep_time.txt', './output/test_1') # test input (comment out after test)
param <- read.csv(file = file.path(args[2], "param.csv"), row.names = 1)

del_array <- fread(file = "./input/del.txt", sep = "\t", nThread = 2) # deletion array
input <- read.csv(file = args[1], sep = "\t") # trait data
dat.use <- merge(input, del_array, by = "ID") # merge data (deletion array and trait data)
dat.use <- na.omit(dat.use)
dat.use[, 2] <- scale(dat.use[, 2])
x <- dat.use[, 2]
y <- dat.use[, -c(1, 2)]
input_2 <- args[2]

df <- read.csv(file.path(input_2, "cor_sig_local.csv"))
colnames(df)[1] <- c("gene")

param <- param[which(rownames(param) %in% df$gene), ]

for (i in 1: nrow(df)) {
    idx <- which(colnames(dat.use) == df$gene[i])
    plot(x, y[, idx], main = df$gene[i], pch = 19, cex = 0.5)
    abline(coef = param[i, c(2, 3)], col = 'green')
    abline(coef = param[i, c(6, 7)], col = 'red')
}