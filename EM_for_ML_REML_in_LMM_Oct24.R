#' This script will implement the EM algorithm for Linear Mixed Models. 
#' Two different EM algorithms are listed here:
#'     1. Maximum likelihood (ML) EM;
#'     2. Restricted Maximum likelihood (REML) EM;
#'     
#'  I also implement an alternative, old-fashioned "E[u|y]" algorithm to estimate 
#'    variance component and beta. 
#' 
#'  Note: 
#'  1. I only implement the one-factor scenario to the above algorithms, i.e. 
#'    One random effect only. 
#'  2. The way we calculate the convergence value: 
#'    sqrt( (sigma_u - sigma_u_old)^2 + (sigma_e - sigma_e_old)^2 ) 
#' 
#' Date: Oct 24, 2018
#' 

EM_LMM <- function(X, Z, y, theta = c(0.5, 0.5), max_iter = 100, thres = 1e-6, type="ML"){
    # Parameters: 
    # X: the incidence matrix for fixed effect;
    # Z: the incidence matrix for one random effect;
    # y: the vector of the realized/observed outcomes;
    # theta: 3-element vector, the starting values for sigma_u and sigma_e, and beta;
    # max_iter: the maximum iteration number before reaching convergence;
    # thres: the convergence threshold;
    # type: a string, "ML" for ML algorithm or "REML" for REML algorithm
    
    sigma_u = theta[1]
    sigma_e = theta[2]
    num_level = ncol(Z)
    num_y = nrow(y)
    
    ZZt = tcrossprod(Z)
    
    i_iter = 1
    diff = 1
    
    if ( type == "ML"){
        while(i_iter <= max_iter & diff > thres){
            
            V = ZZt*sigma_u + diag(sigma_e, ncol = num_y, nrow  = num_y)
            V_inv = solve(V)
            XtVX_inv = solve(t(X)%*%V_inv%*%X)
            
            P = V_inv - V_inv %*% X %*% XtVX_inv %*% t(X) %*% V_inv
            ytP = t(y) %*% P
            sigma_u_new = sigma_u + (sigma_u^2/num_level) * (ytP %*% ZZt %*% t(ytP) - sum(diag(t(Z) %*% V_inv %*% Z)) )
            sigma_e_new = sigma_e + (sigma_e^2/num_y) * (ytP %*% t(ytP) - sum(diag(V_inv)) )
            
            diff = sqrt( (sigma_u - sigma_u_new)^2 + (sigma_e - sigma_e_new)^2 )
            sigma_u = as.numeric(sigma_u_new)
            sigma_e = as.numeric(sigma_e_new)
            
            cat("The", i_iter, "iteration produces sigma_u =", sigma_u, "and sigma_e =", sigma_e, " \n")
            i_iter = i_iter + 1
        }
    }else if (type == "REML") {
        while(i_iter <= max_iter & diff > thres){
            
            V = ZZt*sigma_u + diag(sigma_e, ncol = num_y, nrow  = num_y)
            V_inv = solve(V)
            XtVX_inv = solve(t(X)%*%V_inv%*%X)
            
            P = V_inv - V_inv %*% X %*% XtVX_inv %*% t(X) %*% V_inv
            ytP = t(y) %*% P
            sigma_u_new = sigma_u + (sigma_u^2/num_level) * (ytP %*% ZZt %*% t(ytP) - sum(diag(t(Z) %*% P %*% Z)) )
            sigma_e_new = sigma_e + (sigma_e^2/num_y) * (ytP %*% t(ytP) - sum(diag(P)) )
            
            diff = sqrt( (sigma_u - sigma_u_new)^2 + (sigma_e - sigma_e_new)^2 )
            sigma_u = as.numeric(sigma_u_new)
            sigma_e = as.numeric(sigma_e_new)
            
            cat("The", i_iter, "iteration produces sigma_u =", sigma_u, "and sigma_e =", sigma_e, " \n")
            i_iter = i_iter + 1
        }
    }
    
    V = ZZt*sigma_u + diag(sigma_e, ncol = num_y, nrow  = num_y)
    V_inv = solve(V)
    
    beta = solve(t(X) %*% V_inv %*%X) %*% t(X) %*% V_inv %*% y
    
    u = sigma_u * t(Z) %*% V_inv %*% y
    
    res = list()
    res$u = u
    res$sigma = cbind(sigma_u, sigma_e )
    res$beta = c(beta)
    # names(res) = c("sigma_u", "sigma_e", "beta")
    return(res)
    
}


E_LMM <- function(X, Z, y, theta = c(0.5, 0.5), max_iter = 100, thres = 1e-6){
    sigma_u = theta[1]
    sigma_e = theta[2]
    num_y = nrow(y)
    num_level = ncol(Z)
    ZZt = tcrossprod(Z)
    
    diff = 1
    i_iter = 1
    while (i_iter <= max_iter & diff >= thres){
        V = ZZt*sigma_u + diag(sigma_e, nrow = num_y, ncol = num_y)
        V_inv = solve(V)
        XtVX_inv = solve(t(X)%*%V_inv%*%X)
        P = V_inv - V_inv %*% X %*% XtVX_inv %*% t(X) %*% V_inv
        
        u = sigma_u * t(Z) %*% P %*% y
        e = sigma_e * P %*% y
        
        sigma_u_new = crossprod(u)/(sigma_u * sum(diag(P %*% ZZt)))
        sigma_e_new = crossprod(e)/(sigma_e * sum(diag(P)))
        
        diff = sqrt( (sigma_u - sigma_u_new)^2 + (sigma_e - sigma_e_new)^2 )
        sigma_u = as.numeric(sigma_u_new)
        sigma_e = as.numeric(sigma_e_new)
        
        cat("The", i_iter, "iteration produces sigma_u =", sigma_u, "and sigma_e =", sigma_e, " \n")
        i_iter = i_iter + 1
    }
    
    V = ZZt*sigma_u + diag(sigma_e, ncol = num_y, nrow  = num_y)
    V_inv = solve(V)
    beta = solve(t(X) %*% V_inv %*%X) %*% t(X) %*% V_inv %*% y
    
    res = list()
    res$u = u
    res$sigma = cbind(sigma_u, sigma_e)
    res$beta = c(beta)
    
    return(res)
}













