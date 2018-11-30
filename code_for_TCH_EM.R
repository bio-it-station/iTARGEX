#######################################################################################################################################################################################################
############################################ This program is used to analyze the Correlation of Yeast Gene feature and genetic perturbation data. #####################################################
#######################################################################################################################################################################################################

#### load package #####
library(mixtools)
library(Rcpp)
#### input data ####
original_dir <- getwd()
args=commandArgs(trailingOnly = TRUE)
del.use=read.table("del.txt",header = T,stringsAsFactors = F,sep = "\t")
input=read.csv(file = args[1],stringsAsFactors = F,sep = "\t")[,c(1,as.integer(args[2]))]
del.use=cbind(row.names(del.use),del.use)
colnames(del.use)[1]='name'
colnames(input)[1]='name'
dat.use1=merge(input,del.use,by='name',all=T) #merge data (deletion array and the input gene feature data)
#### function of mixture regression(model II)####
regmixEM1 = function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2, 
                      addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, maxit = 10000, 
                      verb = FALSE) 
{
  if(arbmean == FALSE && arbvar == FALSE){
    stop(paste("Must change constraints on beta and/or sigma!","\n"))
  }
  s = sigma
  if (addintercept) {
    x = cbind(1, x)
  }
  n <- length(y)
  p <- ncol(x)
  tmp <- regmix.init(y = y, x = x, lambda = lambda, beta = beta, 
                     s = s, k = k, addintercept = addintercept, arbmean = arbmean, arbvar = arbvar)
  lambda <- tmp$lambda
  beta <- tmp$beta
  s <- tmp$s
  k <- tmp$k
  diff <- 1
  iter <- 0
  xbeta <- x %*% beta
  res <- (y - xbeta)^2
  if(arbmean == FALSE){
    res <- sapply(1:k,function(i) res)
  }
  #    comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[,i * arbmean + (1 - arbmean)], 
  #        s[i * arbvar + (1 - arbvar)]))
  #    comp <- sapply(comp, cbind)
  #    compsum <- apply(comp, 1, sum)
  #    obsloglik <- sum(log(compsum))
  comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
                                                             s^2)))))
  obsloglik <- sum(log(apply(comp, 1, sum)))
  ll <- obsloglik
  z = matrix(nrow = n, ncol = k)
  zz = matrix(nrow = n, ncol = k)
  restarts <- 0
  while (diff > epsilon && iter < maxit) {
    
    
    
    
    
    
    
    
    
    
    
    
    
    ######################  original R code  ################################
    #   for (i in 1:n) {
    #      for (j in 1:k) {
    #        z.denom = c()
    #        for (h in 1:k) {
    #          z.denom = c(z.denom, (lambda[h]/lambda[j]) * (s[j]/s[h]) * exp(-0.5 * ((1/s[h]^2) * res[i, h] - (1/s[j]^2) * res[i, j])))
    #        }
    #        z[i, j] = 1/sum(z.denom)
    #      }
    #    }
    ##########################################################################
    
    
    
    
    ######################  revise C++ code   ###################################
    cppFunction('NumericMatrix cpp_for(int n, int k, NumericVector lamda,NumericVector s, NumericMatrix res, NumericMatrix zz){
                int i, j, h;
                double shit;
                for(i=0;i<n;++i){
                for(j=0;j<k;++j){
                shit=0.0;
                for(h=0;h<k;++h){
                shit=shit+(lamda[h]/lamda[j]) * (s[j]/s[h]) * exp(-0.5 * ((1/pow(s[h],2)) * res(i, h) - (1/pow(s[j],2)) * res(i, j)));
                }
                zz(i,j)=1/shit;
                }
                }
                return zz;
  }')

    z <- cpp_for(n,k,lambda,s,res,zz)
    #############################################################################
    
    
    
    
    
    #   z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)
    z = z/apply(z,1,sum)
    lambda.new <- apply(z, 2, mean)
    if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
      sing <- 1
    }
    else {
      
      if (arbmean == FALSE) {
        if (addintercept) {
          beta.new <- lm(y~x[,-1],weights=apply(t(t(z)/(s^2)),1,sum))$coef
        }
        else beta.new <- lm(y~x-1,weights=apply(t(t(z)/(s^2)),1,sum))$coef
        #            beta.new <- sapply(lm.out, coef)
        #        beta.new1 <- apply(t(apply(z,2,sum)*t(beta.new)),1,sum)/n
        # beta.new2 <- lm(y~x[,-1],weights=apply(t(t(z)/(s^2)),1,sum))$coef
        # beta.new<-as.vector(solve(t(x) %*% sweep(x, 1, t(t(z)/(s^2)), "*")) %*% apply(t(t(z)/(s^2))*y*x,2,sum) )
      } else {
        if (addintercept) {
          lm.out <- lapply(1:k, function(i) lm(y ~ x[,-1], weights = z[, i]))
        }
        beta.new <- sapply(lm.out, coef)
        beta.new[,2] <- c(0,0)
      }
      xbeta.new <- x %*% beta.new
      res <- (y - xbeta.new)^2
      if(arbmean == FALSE){
        res <- sapply(1:k,function(i) res)
      }
      if (arbvar) {
        s.new <- sqrt(sapply(1:k, function(i) sum(z[, 
                                                    i] * (res[, i]))/sum(z[, i])))
      }
      else s.new <- sqrt(sum(z * res)/n)
      lambda <- lambda.new
      beta <- beta.new
      xbeta <- x%*%beta
      s <- s.new
      sing <- sum(s < 1e-08)
      comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[,i * arbmean + (1 - arbmean)], 
                                                        s[i * arbvar + (1 - arbvar)]))
      comp <- sapply(comp, cbind)
      compsum <- apply(comp, 1, sum)
      newobsloglik <- sum(log(compsum))
      #            comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
      #                s^2)))))
      #            newobsloglik <- sum(log(apply(comp, 1, sum)))
    }
    if (sing > 0 || is.na(newobsloglik) || newobsloglik < obsloglik || abs(newobsloglik) == 
        Inf){# || sum(z) != n) {
      cat("Need new starting values due to singularity...", 
          "\n")
      restarts <- restarts + 1
      if(restarts>15) stop("Too many tries!")
      tmp <- regmix.init(y = y, x = x, k = k, addintercept = addintercept, 
                         arbmean = arbmean, arbvar = arbvar)
      lambda <- tmp$lambda
      beta <- tmp$beta
      s <- tmp$s
      k <- tmp$k
      diff <- 1
      iter <- 0
      xbeta <- x %*% beta
      res <- (y - xbeta)^2
      if(arbmean == FALSE){
        res <- sapply(1:k,function(i) res)
      }
      #    comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[,i * arbmean + (1 - arbmean)], 
      #        s[i * arbvar + (1 - arbvar)]))
      #    comp <- sapply(comp, cbind)
      #    compsum <- apply(comp, 1, sum)
      #    obsloglik <- sum(log(compsum))
      comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
                                                                 s^2)))))
      obsloglik <- sum(log(apply(comp, 1, sum)))
      ll <- obsloglik
      
    }
    else {
      diff <- newobsloglik - obsloglik
      obsloglik <- newobsloglik
      ll <- c(ll, obsloglik)
      iter <- iter + 1
      if (verb) {
        cat("iteration=", iter, "diff=", diff, "log-likelihood", 
            obsloglik, "\n")
      }
    }
}
  scale.order = order(s)
  sigma.min = min(s)
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  if(arbmean == FALSE){
    z=z[,scale.order]
    names(beta) <- c(paste("beta", ".", 0:(p-1), sep = ""))
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, y=y, lambda = lambda[scale.order], beta = beta, sigma = sigma.min, scale = s[scale.order]/sigma.min, loglik = obsloglik, 
           posterior = z[,scale.order], all.loglik=ll, restarts = restarts, ft="regmixEM")
    class(a) = "mixEM"
    a
  } else {
    rownames(beta) <- c(paste("beta", ".", 0:(p-1), sep = ""))
    colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, y=y, lambda = lambda, beta = beta, sigma = s, loglik = obsloglik, 
           posterior = z, all.loglik=ll, restarts = restarts, ft="regmixEM")
    class(a) = "mixEM"
    a
  }
  }

#### function of model selection ####
regmixmodel.sel.wy = function (x, y, w=NULL, k = 2, ...) 
{
  aic <- NULL
  bic <- NULL
  caic <- NULL
  icl <- NULL
  AIC <- function(emout, p) {
    emout$loglik - (p - 1)
  }
  BIC <- function(emout, p, n) {
    emout$loglik - log(n) * (p - 1)/2
  }
  CAIC <- function(emout, p, n) {
    emout$loglik - (log(n) + 1) * (p - 1)/2
  }
  ICL <- function(emout, p, n) {
    BIC(emout, p, n) - sum(emout$lambda * log(emout$lambda))
  }
  
  p <- function(emout) {length(emout$beta) + length(emout$sigma) + (length(emout$scale)-1)*is.null(emout$scale) +length(emout$lambda)}
  n <- length(y)
  for (i in 1:k) {
    if (i == 1) {
      a <- glm(y ~ x)
      beta <- matrix(a$coef, ncol = 1)
      loglik <- log(prod(dnorm(y, mean = beta[1, ] + 
                                 as.matrix(x) %*% beta[2:nrow(beta), ], sd = sd(a$res))))
      emout <- list(beta = beta, sigma = sd(a$res), 
                    lambda = 1, loglik = loglik)
    }
    else emout <- regmixEM1(y, x, k = i, ...)
    P = p(emout)  
    aic[i] <- AIC(emout, p=P)
    bic[i] <- BIC(emout, p=P, n=n)
    caic[i] <- CAIC(emout, p=P, n=n)
    icl[i] <- ICL(emout, p=P, n=n)
  }
  
  out = rbind(aic, bic, caic, icl)
  Winner = apply(out, 1, function(x) (1:length(x))[x == max(x)])
  colnames(out) = 1:k
  rownames(out) = c("AIC", "BIC", "CAIC", "ICL")
  out1=cbind(out, Winner)
  out2=list("model_selection"=out1,"regmixEM_res"=emout)
  out2
}

#### model selection ####
#setwd("~/yeast/")
output_dir <- args[3]
dir.create(output_dir) #create the directory for analysis
setwd(output_dir)
for(i in 3:length(dat.use1)){
  set.seed(10)
  test = na.omit(dat.use1[,c(2,i)])
  test[,1]=scale(test[,1])
  test[,2]=as.numeric(paste(test[,2]))
  model_select = regmixmodel.sel.wy(x=test[,1],y=test[,2],k=2) #model selection
  model_II_parameter=data.frame(comp1=c(model_select$regmixEM_res$lambda[1],model_select$regmixEM_res$beta[,1],model_select$regmixEM_res$sigma[1]),comp2=c(model_select$regmixEM_res$lambda[2],model_select$regmixEM_res$beta[,2],model_select$regmixEM_res$sigma[2]),row.names = c("lambda","beta0","beta1","sigma"))  #result of model II parameters
  write.table(model_II_parameter,file = paste0("./",colnames(dat.use1)[i],".txt"),quote = F,row.names = T) #result of model II parameters
  write.table(model_select$model_selection,file = paste0("./",colnames(dat.use1)[i],"_model_selection",".txt"),quote = F) #result of model selection for each mutant 
  write.table(model_select$regmixEM_res$posterior,file = paste0("./",colnames(dat.use1)[i],"_weight.txt"),quote = F,row.names = F) #result of posterior probabilities of 6000 genes for model II
}
write.table(na.omit(dat.use1)[,1], "gene_order.txt", quote = F, col.names = F, row.names = F)

#### correlation computation ####
correlation_computation=function(dat.use1){
  #Global effect correlation computation
  cor_global=c()
  for(i in 3:length(dat.use1)){
    test = na.omit(dat.use1[,c(2,i)])
    test[,1]=scale(test[,1])
    test[,2]=as.numeric(paste(test[,2]))
    cor_glo=cor(test[,1],test[,2]) #correlation computation
    cor_glo_p.val=2*min(pt(cor_glo*sqrt(nrow(test)-2)/sqrt(1-cor_glo^2),nrow(test)-2,lower.tail = F),pt(cor_glo*sqrt(nrow(test)-2)/sqrt(1-cor_glo^2),nrow(test)-2)) #p-value computation
    cor_pre=c(cor_glo,cor_glo_p.val)
    cor_global=rbind(cor_global,cor_pre) #result of global correlation
    row.names(cor_global)[i-2]=colnames(dat.use1)[i]
  }
  #Local effect correlation computation
  cor_local=c()
  for(i in 3:length(dat.use1)){
    test = na.omit(dat.use1[,c(2,i)])
    test[,1]=scale(test[,1])
    test[,2]=as.numeric(paste(test[,2]))
    weight=read.table(file = paste0("./",colnames(dat.use1)[i],"_weight.txt"),header = T)
    cor_comp1=cor(test[which(weight[,1]>=0.5),1]*weight[which(weight[,1]>=0.5),1],test[which(weight[,1]>=0.5),2]*weight[which(weight[,1]>=0.5),1]) #correlation_component1 computation
    cor_comp2=cor(test[which(weight[,2]>=0.5),1]*weight[which(weight[,2]>=0.5),2],test[which(weight[,2]>=0.5),2]*weight[which(weight[,2]>=0.5),2]) #correlation_component2 computation
    ratio=read.table(file = paste0("./",colnames(dat.use1)[i],".txt"),header = T) #loading the ratio of responsive gene for each mutant
    cor_comp1_p.val=2*min(pt(cor_comp1*sqrt(nrow(test[which(weight[,1]>=0.5),1])-2)/sqrt(1-cor_comp1^2),nrow(test[which(weight[,1]>=0.5),1])-2,lower.tail = F),pt(cor_comp1*sqrt(nrow(test[which(weight[,1]>=0.5),1])-2)/sqrt(1-cor_comp1^2),nrow(test[which(weight[,1]>=0.5),1])-2)) #p-value_component1 computation
    cor_comp2_p.val=2*min(pt(cor_comp2*sqrt(nrow(test[which(weight[,2]>=0.5),1])-2)/sqrt(1-cor_comp2^2),nrow(test[which(weight[,2]>=0.5),1])-2,lower.tail = F),pt(cor_comp2*sqrt(nrow(test[which(weight[,2]>=0.5),1])-2)/sqrt(1-cor_comp2^2),nrow(test[which(weight[,2]>=0.5),1])-2)) #p-value_component2 computation
    corr_pre=cbind(cor_comp1,cor_comp1_p.val,cor_comp2,cor_comp2_p.val,ratio[1,])
    cor_local=rbind(cor_local,corr_pre) #result of local correlation
    row.names(cor_local)[i-2]=colnames(dat.use1)[i]
  }
  #correlation_merge
  cor_all=cbind(cor_local[,c(1,2,5)],cor_global) #combine the correlation of global and local effect for each mutant
  colnames(cor_all)=c("cor_local","pval_local","ratio of responsive genes","cor_global","pval_global")
  cor_all[,6]=NA
  colnames(cor_all)[6]="components"
  for(i in 3:length(dat.use1)){
    aa = read.table(file = paste0("./",colnames(dat.use1)[i],"_model_selection",".txt"),sep = " ",stringsAsFactors = F)
    cor_all[(i-2),6]=aa[2,3]
  } #add the column for the result of model selection
  cor_final=matrix(NA,nrow(cor_all),3)
  for(i in 1:nrow(cor_all)){
    if(cor_all[i,6]==1){
      cor_final[i,1]=cor_all[i,4]
      cor_final[i,2]=cor_all[i,5]
      cor_final[i,3]=1
    }
    if(cor_all[i,6]==2){
      cor_final[i,1]=cor_all[i,1]
      cor_final[i,2]=cor_all[i,2]
      cor_final[i,3]=cor_all[i,3]
    }
  } #classification of the mutants belong to model1 or model2
  row.names(cor_final)=row.names(cor_all)
  colnames(cor_final)=c("correlation","pvalue","ratio of responsive genes")
  list(cor_final=cor_final,cor_all=cor_all)
}
cor_final_all=correlation_computation(dat.use1)
#### correlation significance ####
gene_description <- file.path(original_dir, '1-s2.0-S0092867414003420-mmc1.csv')
correlation_significance=function(cor_final_all){
  clus = read.csv(gene_description, header = T,stringsAsFactors = F) #functional clustering and functional description of the mutants
  cor_final=cor_final_all[[1]]
  for (i in 1:nrow(cor_final)) {
    row.names(cor_final)[i] = strsplit(row.names(cor_final)[i],"[.]")[[1]][1]
  }
  for (i in 1:nrow(cor_final)) {
    row.names(cor_final)[i] = toupper(row.names(cor_final)[i])
  }#change the lowercase to uppercase
  rownames(cor_final)[rownames(cor_final)=="HPA1"]="ELP3"
  rownames(cor_final)[rownames(cor_final)=="CAF17"]="IBA57"
  rownames(cor_final)[rownames(cor_final)=="CAP"]="SRV2"
  rownames(cor_final)[rownames(cor_final)=="CDK8"]="SSN3"
  rownames(cor_final)[rownames(cor_final)=="CYCC"]="SSN8"
  rownames(cor_final)[rownames(cor_final)=="FMP38"]="GEP3"
  rownames(cor_final)[rownames(cor_final)=="MED12"]="SRB8"
  rownames(cor_final)[rownames(cor_final)=="MED13"]="SSN2"
  rownames(cor_final)[rownames(cor_final)=="KRH2"]="GPB1"
  rownames(cor_final)[rownames(cor_final)=="KRH1"]="GPB2"
  rownames(cor_final)[rownames(cor_final)=="KEM1"]="XRN1"
  rownames(cor_final)[rownames(cor_final)=="YNR004W"]="SWM2"
  rownames(cor_final)[rownames(cor_final)=="YNL213C"]="RRG9"
  rownames(cor_final)[rownames(cor_final)=="YMR293C"]="HER2"
  rownames(cor_final)[rownames(cor_final)=="YER064C"]="VHR2"
  rownames(cor_final)[rownames(cor_final)=="YOL138C"]="RTC1"
  rownames(cor_final)[rownames(cor_final)=="CHF1"]="DMA1"
  rownames(cor_final)[rownames(cor_final)=="RIS1"]="ULS1"
  rownames(cor_final)[rownames(cor_final)=="YHL010C"]="ETP1"
  rownames(cor_final)[rownames(cor_final)=="SSN6"]="CYC8"
  rownames(cor_final)[rownames(cor_final)=="MED31"]="SOH1"
  rownames(cor_final)[rownames(cor_final)=="MED5"]="NUT1"
  rownames(cor_final)[rownames(cor_final)=="MED9"]="CSE2"
  rownames(cor_final)[rownames(cor_final)=="MED3"]="PGD1"
  rownames(cor_final)[rownames(cor_final)=="MED18"]="SRB5"
  rownames(cor_final)[rownames(cor_final)=="MED20"]="SRB2"
  rownames(cor_final)[rownames(cor_final)=="MED15"]="GAL11"
  rownames(cor_final)[rownames(cor_final)=="MED16"]="SIN4"
  rownames(cor_final)[rownames(cor_final)=="CAC1"]="RLF2"
  rownames(cor_final)[rownames(cor_final)=="YGR071C"]="ENV11"
  rownames(cor_final)[rownames(cor_final)=="YDR266C"]="HEL2"
  rownames(cor_final)[rownames(cor_final)=="NOT4"]="MOT2"
  rownames(cor_final)[rownames(cor_final)=="YIL014C"]="YIL014C-A"
  rownames(cor_final)[rownames(cor_final)=="YBR028C"]="YPK3"
  rownames(cor_final)[rownames(cor_final)=="ATG24"]="SNX4"
  rownames(cor_final)[rownames(cor_final)=="YNR047W"]="FPK1"
  rownames(cor_final)[rownames(cor_final)=="YKL161C"]="KDX1"
  cor_final=cbind(row.names(cor_final),cor_final)
  colnames(cor_final)[1]="gene"
  cor_final=merge(cor_final,clus[,c(1,3)],by="gene")
  cor_global_only = cor_final[cor_final[,4]==1,] #global_correlation
  cor_local_only = cor_final[cor_final[,4]!=1,] #local_correlation
  write.table(cor_global_only[,-4],file = "./cor_global_only.txt",sep = "\t",quote = F,row.names = F) #gene list recording the genetic perturbation which has global effect on the input gene feature
  write.table(cor_local_only,file = "./cor_local_only.txt",sep = "\t",quote = F,row.names = F) #gene list recording the genetic perturbation which has local effect on the input gene feature
  
  cor_final_sig=cor_final[as.numeric(paste(cor_final[,3]))<=0.01/nrow(cor_final_all[[2]]),] #bonferroni correction at alpha=0.01
  cor_final_sig[,3]=-log10(as.numeric(paste(cor_final_sig[,3]))) #-log10(p)
  cor_final_sig_global = subset(cor_final_sig,cor_final_sig[,4]==1) #global effect significant gene
  cor_final_sig_local = subset(cor_final_sig,cor_final_sig[,4]!=1) #local effect significant gene
  if (nrow(cor_final_sig_global)) {
    cor_final_sig_global[!is.finite(as.numeric(paste(cor_final_sig_global[,3]))),3]=350 #change Inf p-value to 350 (global)
  }
  if (nrow(cor_final_sig_local)) {
    cor_final_sig_local[!is.finite(as.numeric(paste(cor_final_sig_local[,3]))),3]=350 #change Inf p-value to 350 (local)
  }
  write.table(cor_final_sig_global[,-4],file = "./cor_sig_global_only.txt",sep = "\t",quote = F,row.names = F) #gene list recording the genetic perturbation which has global effect on the input gene feature
  write.table(cor_final_sig_local,file = "./cor_sig_local_only.txt",sep = "\t",quote = F,row.names = F) #gene list recording the genetic perturbation which has local effect on the input gene feature
}
correlation_significance(cor_final_all)

#### manhatten like plot ####
cor.global=read.delim("./cor_global_only.txt",header = T,stringsAsFactors = F)
cor.local=read.delim("./cor_local_only.txt",header = T,stringsAsFactors = F)
clus = read.csv(gene_description, header = T,stringsAsFactors = F) #functional clustering and functional description of the mutants
atr=clus[clus$gene%in%union(cor.global[,1],cor.local[,1]),]
dat.global=merge(cor.global,atr[c(1,2,4)],by= "gene")
dat.global$logP=-log10(dat.global$pvalue)
dat.global$logP[dat.global$logP==Inf]=350
dat.local=merge(cor.local,atr[c(1,2,4)],by = "gene")
dat.local$logP=-log10(dat.local$pvalue)
dat.local$logP[dat.local$logP==Inf]=350
png(paste0("./manhatten_like_plot_",colnames(dat.use1)[2],".png"),width = 20,height = 12,units = "cm",pointsize = 10,res = 600)
par(mar=c(14,4,4,2))
plot(1,1,xlim = c(1,16),ylim = c(0,max(c(max(dat.global$logP),max(dat.local$logP)))+10),type = 'n',xlab = "",ylab = "-log10(p)",axes = F)
rect(seq(0.5,15.5,2),rep(0,8),seq(1.5,16.5,2),rep(400,8),col = "grey70",border = NA)
for(i in 1:length(unique(atr$functional.category))){
  index1=dat.global$functional.category%in%unique(atr$functional.category)[i]
  index2=dat.local$functional.category%in%unique(atr$functional.category)[i]
  points(rnorm(i-0.15,0.05,n = sum(index1)),dat.global$logP[index1],col=ifelse(dat.global$logP[index1]<summary(dat.global$logP)[5],"black","orange"),cex=0.7,pch=19)
  points(rnorm(i+0.15,0.05,n = sum(index2)),dat.local$logP[index2],col=ifelse(dat.local$logP[index2]<summary(dat.local$logP)[5],"black","green"),cex=0.7,pch=19)
}
axis(1,at = 1:16,labels = unique(atr$functional.category),las=2)
axis(2)
abline(h = summary(dat.global$logP)[5],col="red")
abline(h = summary(dat.local$logP)[5],col="blue")
dev.off()

#enrichment test-global (top20)
GO_result_global=c()
for (i in 1:length(goterm)) {
  m=length(unique(sgd[sgd$V5==goterm[i],3]))
  k=length(cor.global[order(cor.global[,3])[1:20],1])
  x=length(intersect(k,m))
  n=length(union(cor.global[,1],cor.local[,1]))-m
  p.val = phyper(x,m,n,k,lower.tail = F)
  pre=c(goterm[i],p.val)
  GO_result_global=rbind(GO_result_global,pre)
}
colnames(GO_result_global)=c("GO_term","p.value")
write.table(GO_result_global,file = "./GO_result_global.txt",quote = F,row.names = F)

#enrichment test-local (top20)
GO_result_local=c()
for (i in 1:length(goterm)) {
  m=length(unique(sgd[sgd$V5==goterm[i],3]))
  k=length(cor.local[order(cor.local[,3])[1:20],1])
  x=length(intersect(k,m))
  n=length(union(cor.global[,1],cor.local[,1]))-m
  p.val = phyper(x,m,n,k,lower.tail = F)
  pre=c(goterm[i],p.val)
  GO_result_local=rbind(GO_result_local,pre)
}
colnames(GO_result_local)=c("GO_term","p.value")
write.table(GO_result_local,file = "./GO_result_local.txt",quote = F,row.names = F)
