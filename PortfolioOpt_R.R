library(nloptr)
library(data.table)

run_opt <- function(returns, cov_matrix, boundDat, minTot){
  spp <- names(cov_matrix)
  sppUse <- spp
  mean_returns <- colMeans(returns)
  bounds <- boundDat[Spp %chin% spp,] 
  bndNew <- bounds
  
  while(length(mean_returns) > 1){
    ##maybe temp
    target <- set_target(mean_returns,cov_matrix)
  }
}


set_target <- function(mean_returns, cov_matrix){
  num_ass <- length(mean_returns)
  eq_constr <- function(x) sum(x) - 1
  min_var_w <- slsqp(x0 = rep(1/num_ass,num_ass),fn = portfolio_volatility,
                     lower = rep(0,num_ass), upper = rep(1,num_ass),
                     heq = eq_constr, mean_returns = mean_returns, cov_matrix = cov_matrix)
  min_var <- portfolio_return(min_var_w$par,mean_returns)
  target <- seq(min_var,max(mean_returns),length.out = 20)
  return(target)
}

portfolio_volatility <- function(weights, mean_returns, cov_matrix){
  #weights <- as.matrix(weights)
  return(t(weights) %*% cov_matrix %*% weights)
}

portfolio_return <- function(weights, mean_returns){
  return(sum(weights*mean_returns))
}

efficient_return <- function(mean_returns, cov_matrix, target, bounds){
  
}

portfolio_stdev <- function(weights){
  
}