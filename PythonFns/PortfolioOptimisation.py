# -*- coding: utf-8 -*-
"""
Portfolio Optimisation for Tree Species

@author: Kiri Daust
Many functions based on Ricky Kim's script:
https://towardsdatascience.com/efficient-frontier-portfolio-optimisation-in-python-e7844051e7f
"""

import pandas as pd  
import numpy as np
import scipy.optimize as sco
#import math

#df = pd.read_csv("TestReturns.csv")
#cov_matrix = pd.read_csv("TestSigma.csv",index_col = 0)
#spp = cov_matrix.columns
#returns = df
#mean_returns = returns.mean()

#risk_free_rate = 0

###setup basic functions##########
def portfolio_volatility(weights, mean_returns, cov_matrix): ##for minimising stdev
    return weights.T @ cov_matrix @ weights
    
def portfolio_return(weights, mean_returns, cov_matrix): ##for maximising return
    returns = np.sum(mean_returns*weights) 
    return returns
    
def neg_port_return(weights, mean_returns, cov_matrix): ##for maximising return
    #retAdj = [x - 2 if x < 3 else x for x in mean_returns]
    retAdj = mean_returns
    returns = np.sum(retAdj*weights) 
    return -(returns)

def set_bounds(boundDat, spp):
    tempInd = boundDat['Spp'].isin(spp)
    newBnd = boundDat[tempInd]
    minV = newBnd['minWt']
    maxV = newBnd['maxWt']
    return [(mi,ma) for mi,ma in zip(minV,maxV)]
    

###calculate efficient frontier points
def efficient_return(mean_returns, cov_matrix, target, bounds): ###efficient frontier for stdev
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix)

    constraints = ({'type': 'eq', 'fun': lambda x: np.sum(mean_returns*x) - target},
                   {'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
    result = sco.minimize(portfolio_volatility, num_assets*[1./num_assets,], args=args, method='SLSQP', bounds=bounds, constraints=constraints)
    return result
    
def efficient_stdev(mean_returns, cov_matrix, target, bounds): ##efficient frontier for return
   num_assets = len(mean_returns)
   args = (mean_returns, cov_matrix)

   def portfolio_stdev(weights): ###stdev
       return np.sqrt(np.dot(weights.T, np.dot(cov_matrix, weights)))

   constraints = ({'type': 'eq', 'fun': lambda x: portfolio_stdev(x) - target},##stdev == target
                  {'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
   result = sco.minimize(neg_port_return, num_assets*[1./num_assets,], args=args, method='SLSQP', bounds=bounds, constraints=constraints, jac = '2-point')
   return result
    
def set_target(mean_returns, cov_matrix):
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix)
    constraints = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1}) ##have to sum to 1
    bound = (0.0,1.0)
    bounds = tuple(bound for asset in range(num_assets))
    min_var_w = sco.minimize(portfolio_volatility, num_assets*[1./num_assets,], args=args,
                        method='SLSQP', bounds=bounds, constraints=constraints)
    min_var = portfolio_return(min_var_w['x'], mean_returns, cov_matrix)
    min_round = min_var
    max_round = max(mean_returns)
    target = np.linspace(min_round, max_round, 20)
    return target

def ef_weights_v2(returns, cov_matrix, boundDat, minTot): ##optType either mean or stdev - main R function
    spp = cov_matrix.columns
    sppUse = spp
    mean_returns = returns.mean()
    bounds = set_bounds(boundDat, spp)
    bndNew = bounds
    
    while(len(mean_returns) > 1):
        temp = []
        target = set_target(mean_returns,cov_matrix)
        testRet = [target[0]+0.5,target[-1]-0.5]
        for ret in testRet:
            temp.append(efficient_return(mean_returns, cov_matrix, ret, bndNew))
        w = np.array([x['x'] for x in temp])
        if(any(np.max(w, axis = 0) < minTot)):
            sppUse = np.array(list(sppUse[np.max(w, axis = 0) > minTot]))###remove spcies from covmat, returns, and bounds
            cov_matrix = cov_matrix.loc[sppUse, sppUse]
            mean_returns = mean_returns[sppUse]
            indUse = [b for a,b in zip(spp, range(len(spp))) if a in sppUse]
            bndNew = [bounds[i] for i in indUse]
            
        else:
            break
    
    efficients = []
    min_risk = target[0]
    for ret in target:
        efficients.append(efficient_return(mean_returns, cov_matrix, ret, bndNew))
    ep_w = [x['x'] for x in efficients]
    w_df = np.array(ep_w)
    w_df = pd.DataFrame(w_df)
    w_df.columns = sppUse
    ep_stdev = np.array([x['fun'] for x in efficients])
    ep_ret = np.array(target)
    sharpe = (ep_ret-min_risk)/ep_stdev
    return(w_df,ep_stdev,ep_ret,sharpe)


####efficient frontier - main function in R
#def ef_weights(returns, cov_matrix, target, minWt, maxWt, minTot): ##optType either mean or stdev - main R function
#    spp = cov_matrix.columns
#    sppUse = spp
#    mean_returns = returns.mean()
#    bounds = set_bounds(minWt, maxWt, mean_returns)
#    bndNew = bounds
#    while(True): ###if any < minTot, loop through and remove
#        efficients = []
#        for ret in target:
#            efficients.append(efficient_return(mean_returns, cov_matrix, ret, bndNew))
#        ep_w = [x['x'] for x in efficients]
#        w_df = np.array(ep_w)
#        if(any(np.max(w_df, axis = 0) < minTot)):
#            sppUse = list(spp[np.max(w_df, axis = 0) > minTot])###remove spcies from covmat, returns, and bounds
#            cov_matrix = cov_matrix.loc[sppUse, sppUse]
#            mean_returns = mean_returns[sppUse]
#            indUse = [b for a,b in zip(spp, range(len(spp))) if a in sppUse]
#            bndNew = [bounds[i] for i in indUse]
#            
#        else:
#            break
#        
#    w_df = pd.DataFrame(w_df)
#    w_df.columns = sppUse
#    ep_optVal = [portfolio_volatility(x['x'], mean_returns, cov_matrix) for x in efficients]        
#    ep_optVal = np.array(ep_optVal)
#    
#    return(w_df,ep_optVal)

