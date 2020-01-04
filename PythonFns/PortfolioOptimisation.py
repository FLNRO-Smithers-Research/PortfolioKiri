# -*- coding: utf-8 -*-
"""
Portfolio Optimisation for Tree Species

@author: Kiri Daust
Many functions based on Ricky Kim's script:
https://towardsdatascience.com/efficient-frontier-portfolio-optimisation-in-python-e7844051e7f
"""

import pandas as pd  
import numpy as np
#import matplotlib.pyplot as plt
import scipy.optimize as sco
import math

#df = pd.read_csv("TestReturns.csv")
#spp = cov_matrix.columns
#returns = df
#mean_returns = returns.mean()
#cov_matrix = pd.read_csv("TestSigma.csv",index_col = 0)
#risk_free_rate = 0

###setup basic functions##########
def portfolio_volatility(weights, mean_returns, cov_matrix): ##for minimising stdev
    return np.sqrt(np.dot(weights.T, np.dot(cov_matrix, weights)))
    
def portfolio_return(weights, mean_returns, cov_matrix): ##for maximising return
    returns = np.sum(mean_returns*weights) 
    return returns
    
def neg_port_return(weights, mean_returns, cov_matrix): ##for maximising return
    #retAdj = [x - 2 if x < 3 else x for x in mean_returns]
    retAdj = mean_returns
    returns = np.sum(retAdj*weights) 
    return -(returns)

def set_bounds(minV, maxV, mean_returns):
    if(isinstance(minV, (float,int)) and isinstance(maxV, (float,int))):
        return [(minV,maxV) for asset in range(len(mean_returns))]
    if(len(minV) != len(maxV) or len(minV) !=  len(mean_returns)):
        return -1
    return [(mi,ma) for mi,ma in zip(minV,maxV)]
    

###calculate efficient frontier points
def efficient_return(mean_returns, cov_matrix, target, bounds): ###efficient frontier for stdev
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix)

    def portfolio_return(weights):
        return np.sum(mean_returns*weights)

    constraints = ({'type': 'eq', 'fun': lambda x: portfolio_return(x) - target},
        # {'type': 'ineq', 'fun': lambda x: portfolio_return(x) - (target+0.01)},
        #             {'type': 'ineq', 'fun': lambda x: (target-0.01) - portfolio_return(x)},
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
    
def set_target(returns, cov_matrix):
    mean_returns = returns.mean()
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix)
    constraints = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
    bound = (0.0,1.0)
    bounds = tuple(bound for asset in range(num_assets))
    min_var_w = sco.minimize(portfolio_volatility, num_assets*[1./num_assets,], args=args,
                        method='SLSQP', bounds=bounds, constraints=constraints)
    min_var = portfolio_return(min_var_w['x'], mean_returns, cov_matrix)
    min_round = round(min_var,0)
    max_round = math.ceil(max(mean_returns))
    ##min_round = round(min_var*2,1)/2
    target = np.linspace(min_round, max_round+1, 20)
    return target

def ef_weights_cbst(returns, cov_matrix, target, minWt, maxWt, minTot): ##optType either mean or stdev - main R function
    spp = cov_matrix.columns
    sppUse = spp
    mean_returns = returns.mean()
    bounds = set_bounds(minWt, maxWt, mean_returns)
    bndNew = bounds
    while(len(mean_returns) > 1): ###if any < minTot, loop through and remove
        efficients = []
        for ret in target:
            efficients.append(efficient_stdev(mean_returns, cov_matrix, ret, bndNew))
        ep_w = [x['x'] for x in efficients]
        w_df = np.array(ep_w)
        if(any(np.max(w_df, axis = 0) < minTot)):
            sppUse = np.array(list(sppUse[np.max(w_df, axis = 0) > minTot]))###remove spcies from covmat, returns, and bounds
            cov_matrix = cov_matrix.loc[sppUse, sppUse]
            mean_returns = mean_returns[sppUse]
            indUse = [b for a,b in zip(spp, range(len(spp))) if a in sppUse]
            bndNew = [bounds[i] for i in indUse]
            
        else:
            break
        
    w_df = pd.DataFrame(w_df)
    w_df.columns = sppUse
    ep_optVal = [portfolio_volatility(x['x'], mean_returns, cov_matrix) for x in efficients]        
    ep_optVal = np.array(ep_optVal)
    ep_retVal = [portfolio_return(x['x'], mean_returns, cov_matrix) for x in efficients]        
    ep_retVal = np.array(ep_retVal)
    
    return(w_df,ep_optVal,ep_retVal)


###efficient frontier - main function in R
def ef_weights(returns, cov_matrix, target, minWt, maxWt, minTot): ##optType either mean or stdev - main R function
    spp = cov_matrix.columns
    sppUse = spp
    mean_returns = returns.mean()
    bounds = set_bounds(minWt, maxWt, mean_returns)
    bndNew = bounds
    while(True): ###if any < minTot, loop through and remove
        efficients = []
        for ret in target:
            efficients.append(efficient_return(mean_returns, cov_matrix, ret, bndNew))
        ep_w = [x['x'] for x in efficients]
        w_df = np.array(ep_w)
        if(any(np.max(w_df, axis = 0) < minTot)):
            sppUse = list(spp[np.max(w_df, axis = 0) > minTot])###remove spcies from covmat, returns, and bounds
            cov_matrix = cov_matrix.loc[sppUse, sppUse]
            mean_returns = mean_returns[sppUse]
            indUse = [b for a,b in zip(spp, range(len(spp))) if a in sppUse]
            bndNew = [bounds[i] for i in indUse]
            
        else:
            break
        
    w_df = pd.DataFrame(w_df)
    w_df.columns = sppUse
    ep_optVal = [portfolio_volatility(x['x'], mean_returns, cov_matrix) for x in efficients]        
    ep_optVal = np.array(ep_optVal)
    
    return(w_df,ep_optVal)


#max_sharpe = max_sharpe_ratio(mean_returns, cov_matrix, risk_free_rate)
#sdp, rp = portfolio_annualised_performance(max_sharpe['x'], mean_returns, cov_matrix)
#max_sharpe_allocation = pd.DataFrame(max_sharpe.x,index=returns.columns,columns=['allocation'])
#max_sharpe_allocation.allocation = [round(i*100,2)for i in max_sharpe_allocation.allocation]
#max_sharpe_allocation = max_sharpe_allocation.T
#
#min_vol = min_variance(mean_returns, cov_matrix)
#sdp_min, rp_min = portfolio_annualised_performance(min_vol['x'], mean_returns, cov_matrix)
#min_vol_allocation = pd.DataFrame(min_vol.x,index=returns.columns,columns=['allocation'])
#min_vol_allocation.allocation = [round(i*100,2)for i in min_vol_allocation.allocation]
#min_vol_allocation = min_vol_allocation.T
#
#target = np.linspace(25, 5, 25)
#efficient_portfolios = efficient_frontier(mean_returns, cov_matrix, target)
#ep_w = [x['x'] for x in efficient_portfolios]
#w_df = pd.DataFrame(ep_w)
#w_df.columns = spp
#w_df = w_df.set_index(target)
#ep = [x['fun'] for x in efficient_portfolios]
#plt.figure(figsize=(10, 7))
#plt.plot(ep, target, linestyle='-.', color='black', label='efficient frontier')

# def neg_sharpe_ratio(weights, mean_returns, cov_matrix, risk_free_rate): ###for maximising sharpe ration
#     p_ret = np.sum(mean_returns*weights ) 
#     p_var = np.sqrt(np.dot(weights.T, np.dot(cov_matrix, weights)))
#     return -(p_ret - risk_free_rate) / p_var
# 
# def max_sharpe_ratio(returns, cov_matrix, risk_free_rate): ##weights for max sharpe ratio
#     mean_returns = returns.mean()
#     num_assets = len(mean_returns)
#     args = (mean_returns, cov_matrix, risk_free_rate)
#     constraints = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
#     bound = (0.0,1.0)
#     bounds = tuple(bound for asset in range(num_assets))
#     result = sco.minimize(neg_sharpe_ratio, num_assets*[1./num_assets,], args=args,
#                         method='SLSQP', bounds=bounds, constraints=constraints)
#     r2 = pd.DataFrame(result.x,index=cov_matrix.columns,columns=['Weight'])
#     return r2.T
# 
#else:
#        efficients = []
#        for ret in target:
#            efficients.append(efficient_stdev(mean_returns, cov_matrix, ret, bounds))
#        ep_optVal = [portfolio_return(x['x'], mean_returns, cov_matrix) for x in efficients]
