#' A test script, to simulate some geno pheno data to test the EM algorithms. 
#' 
#' Date: Oct 24, 2018

set.seed(100)

N = 1000  
M = 500       # The number of causal SNPs
M_alt = 0  # The number of null SNPs
M_total = M + M_alt
h2 = 0.3

# the MAF 
freq <- runif(n = M_total, min = 0.01, max = 0.5)

Z <- t(replicate(n = N, expr = rbinom(n = M_total, size = 2, prob = freq)))
Z <- apply(Z, MARGIN = 2, scale)


idx <- sort(sample(1:(M+M_alt), M))
beta <- rnorm(n = M)
Z_causal <- Z[, idx]

g <- Z_causal %*% beta
g <- scale(g)

X <- cbind(1, runif(N, min=-1, max=1))
beta_X <- c(1.0, pi)

# c. residual. 
residual <- rnorm(n = nrow(Z))
residual <- scale(residual)

# combine phenotype
pheno <- scale(sqrt(h2)*g + sqrt(1 - h2)*residual) + scale(X%*%beta_X)


res_ML <- EM_LMM(X, Z, pheno, type="ML")

res_REML <- EM_LMM(X, Z, pheno, type = "REML")

res_E <- E_LMM(X, Z, pheno)

## It seems taht E_LMM converge faster! 



