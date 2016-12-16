#################################################################
## Method of Divide and Combine
## last updated on 2016.12.16
#################################################################

##########################################
## dependent packages
##########################################
# install first if you haven't, use install.packages()
library(foreach) # package for function foreach
library(doSNOW)  # package for parallelization

library(glmnet)  #
library(MASS)    #

##########################################
## data generator, data stored as .RData
##########################################
# description: generate K data sets and store in separate .RData files (each .RData contains variables X and y)
datagenerator <- function(n, p, K, nsig, type, corstruct, rho, tempdatadir, categorical=FALSE, seed=NA){
  n           # sample size in each data set
  p           # dimension of beta
  K           # number of data sets
  nsig        # number of nonzero beta
  type        # c("gaussian", "binomial", "poisson", "cox")
  corstruct   # c("ind", "cs", "ar1") correlation structure of X
  rho         # rho is associated with corstruct
  tempdatadir # name of the director to store simulated data
  categorical # set to TRUE to let some X be dichotomized
  seed        # set random seed
  
  if(length(n)!=1 & length(n)!=K){ stop("n must be a number of a vector of length K.") }
  if(!is.na(seed)){ set.seed(seed) }
  seed.list <- sample(1:1e8, K, replace=FALSE)
  if(length(n)==1){ n <- rep(n,K) }
  
  beta <- rep(0, p)
  ssig <- sort(sample(c(1:p), nsig, replace=FALSE))
  if(type=="gaussian"){ beta[ssig] <- runif(nsig, 0.3, 0.3) }
  if(type=="binomial"){ beta[ssig] <- runif(nsig, 0.3, 0.3) }
  if(type=="poisson") { beta[ssig] <- runif(nsig, 0.1, 0.1) }
  
  dir.create(tempdatadir)
  for(k in 1:K){
    set.seed(seed.list[k])
    if(corstruct=="ind"){ Sigma <- diag(rep(1,p)) }
    if(corstruct=="cs") { Sigma <- matrix(rho,p,p); diag(Sigma) <- 1 }
    if(corstruct=="ar1"){ Sigma <- rho^(abs(outer(1:p, 1:p, "-"))) }
    X <- mvrnorm(n=n[k], mu=rep(0,p), Sigma=Sigma)
    if(categorical==TRUE){ some <- sample(1:p, round(p/2), replace=FALSE); X[,some] <- (X[,some]>0) }
    
    set.seed(seed.list[k] + 1234567)
    if(type=="gaussian"){ y <- X[,ssig]%*%beta[ssig] + rnorm(n[k]) }
    if(type=="binomial"){ expit <- function(a){exp(drop(a))/(1+exp(drop(a)))}; y <- rbinom(n = n[k], size = 1, prob = expit( X[,ssig]%*%beta[ssig] )) }
    if(type=="poisson") { y <- rpois(n = n[k], lambda = exp( X[,ssig]%*%beta[ssig] )) }
    y <- drop(y)
    save(X, y, file = paste(tempdatadir, "/", k, ".RData", sep=""))
  }
  beta
}

# description: read K .RData files and combine
fulldata <- function(K, tempdatadir){
  X.full <- c(); y.full <- c()
  for(k in 1:K){
    load(paste(tempdatadir, "/", k, ".RData", sep=""))
    X.full <- rbind(X.full, X)
    if(length(dim(y))==2){ 
      y.full <- rbind(y.full, y)
    } else{ 
      y.full <- c(y.full, y) 
    }
  }
  list(X=X.full, y=y.full)
}


# description: modac main function
MODAC <- function(K, tempdatadir, type, tau, seed=NA){
  
  if(type=="gaussian"){ invlink <- function(x){x}; invlinkdiv <- function(x){1}; aphi <- function(beta, X, y){out <- sum((y-X%*%beta)^2)/(length(y)-sum(beta!=0))} }
  if(type=="binomial"){ invlink <- function(x){exp(x)/(1 + exp(x))}; invlinkdiv <- function(x){exp(x)/(1 + exp(x))^2}; aphi <- function(beta, X, y){out <- 1} }
  if(type=="poisson") { invlink <- function(x){exp(x)};              invlinkdiv <- function(x){exp(x)}; aphi <- function(beta, X, y){out <- 1} }
  evec <- function(beta, X, y){ y - sapply(drop(X%*%beta), invlink) }
  pmat <- function(beta, X, y){ out <- sapply(drop(X%*%beta), invlinkdiv); diag(out) }

  list_out <- foreach(k = 1:K, .packages = "glmnet") %dopar% {
    load(paste(tempdatadir, "/", k, ".RData", sep=""))
    tau <- 0
    eta1 <- tau
    X <- X[complete.cases(cbind(X,y)),]
    y <- y[complete.cases(cbind(X,y))]
    
    intercept <- FALSE
    fit <- cv.glmnet(X, y, family=type, standardize=FALSE, intercept=intercept)
    if(intercept==FALSE){ betahat <- drop(coef(fit, s="lambda.min"))[-1] }
    if(intercept==TRUE) { betahat <- drop(coef(fit, s="lambda.min")); X <- cbind(1,X) }
    
    n <- length(y); p <- length(betahat)
    if(n <= p){stop("n has to be larger than p")}
    
    a <- aphi(betahat, X, y)
    P <- pmat(betahat, X, y)
    e <- evec(betahat, X, y)
    eta <- rep(eta1,p)
    if(type=="gaussian"){ 
      A <- solve( t(X)%*%X + diag(n*eta) ) %*% t(X)
      betahat.c  <- drop(betahat + A%*%e)
      Sigmahatinv <- (1/a) * t(X)%*%X + diag(eta)
    } else{ 
      A <- solve( t(X)%*%P%*%X + diag(n*eta) ) %*% t(X)
      betahat.c  <- drop(betahat + A%*%e)
      Sigmahatinv <- (1/a) * t(X)%*%P%*%X + diag(eta)
    }
    list(betahat.c=betahat.c, Sigmahatinv=Sigmahatinv)
  }
  length(list_out)
  
  sum1 <- 0; sum2 <- 0
  for(i in 1:K){
    sum1 <- sum1 + list_out[[i]]$Sigmahatinv 
    sum2 <- sum2 + list_out[[i]]$Sigmahatinv %*% list_out[[i]]$betahat.c
  }
  
  varb <- solve(sum1)
  betahat <- drop(varb %*% sum2)
  pvalue <- 2*pnorm(-abs(betahat / sqrt(diag(varb))))
  winner <- (pvalue < 0.05)
  
  cbind(betahat=betahat, stdev=sqrt(diag(varb)), pvalue=pvalue)
}


##########################################
## Test functions
##########################################
N <- 10000  # total sample size
K <- 20     # number of sub-datasets
p <- 200    # number of covariates
nsig <- 10  # number of significant covariates
type <- "gaussian" # response type c("gaussian", "binomial", "poisson")
corstruct <- "cs"  # correlation structure of the generated X
rho <- 0.8         # correlation parameter
n <- round(N/K)    # sample size in sub-dataset
tempdatadir <- 'tempdatafile'   # temporary data storage folder


# tempdatadir <- paste("Temp_", outputfilename, sep="")
beta <- datagenerator(n=n, p=p, K=K, nsig=nsig, type=type, corstruct=corstruct, rho=rho, tempdatadir=tempdatadir)
print(beta[beta!=0])  # print non-zero beta

# GLM ON FULL DATA
time.all.read <- system.time(data.full <- fulldata(K=K, tempdatadir=tempdatadir))[3]  # read in the full data
time.all.glm <- system.time(result.A.glm <- summary(glm(data.full$y ~ data.full$X-1, family=type)))[3] + time.all.read
print(paste("GLM on full data finishes in ", round(time.all.glm,2), " seconds.", sep=""))
rm(data.full)

# MODAC ON DIVIDED DATA
cl = makeCluster(rep("localhost", 4), type = "SOCK") # here 4 is number of processes used on the local machine
registerDoSNOW(cl)
time.MODAC.local.parallel <- system.time(result.P.MODAC <- MODAC(K=K, tempdatadir=tempdatadir, type=type, tau=0))[3]
stopCluster(cl)
print(paste("Under local parallelization (4 cores), MODAC on divided data finishes in ", round(time.MODAC.local.parallel,2), " seconds.", sep=""))

# remove temp data folder
unlink(tempdatadir, recursive = TRUE)  # remove the data












