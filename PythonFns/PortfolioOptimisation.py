# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:43:29 2019

@author: Kiri Daust
"""

import pandas as pd  
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sco

#df = pd.read_csv("TestReturns.csv")
#spp = cov_matrix.columns
#returns = df
#mean_returns = returns.mean()
#cov_matrix = pd.read_csv("CovMatSmal.csv",index_col = 0)
#risk_free_rate = 0

def portfolio_annualised_performance(weights, mean_returns, cov_matrix):
    returns = np.sum(mean_returns*weights ) 
    std = np.sqrt(np.dot(weights.T, np.dot(cov_matrix, weights)))
    return std, returns

def portfolio_volatility(weights, mean_returns, cov_matrix):
    return portfolio_annualised_performance(weights, mean_returns, cov_matrix)[0]

def neg_sharpe_ratio(weights, mean_returns, cov_matrix, risk_free_rate):
    p_var, p_ret = portfolio_annualised_performance(weights, mean_returns, cov_matrix)
    return -(p_ret - risk_free_rate) / p_var

def efficient_return(mean_returns, cov_matrix, target):
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix)

    def portfolio_return(weights):
        return portfolio_annualised_performance(weights, mean_returns, cov_matrix)[1]

    constraints = ({'type': 'eq', 'fun': lambda x: portfolio_return(x) - target},
                   {'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
    bounds = tuple((0,1) for asset in range(num_assets))
    result = sco.minimize(portfolio_volatility, num_assets*[1./num_assets,], args=args, method='SLSQP', bounds=bounds, constraints=constraints)
    return result


def efficient_frontier(mean_returns, cov_matrix, returns_range):
    efficients = []
    for ret in returns_range:
        efficients.append(efficient_return(mean_returns, cov_matrix, ret))
    return efficients

def min_variance(mean_returns, cov_matrix):
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix)
    constraints = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
    bound = (0.0,1.0)
    bounds = tuple(bound for asset in range(num_assets))

    result = sco.minimize(portfolio_volatility, num_assets*[1./num_assets,], args=args,
                        method='SLSQP', bounds=bounds, constraints=constraints)

    return result

def ef_weights(returns, cov_matrix, target, risk_free_rate):
    spp = cov_matrix.columns
    mean_returns = returns.mean()
    efficient_portfolios = efficient_frontier(mean_returns, cov_matrix, target)
    ep_w = [x['x'] for x in efficient_portfolios]
    w_df = pd.DataFrame(ep_w)
    w_df.columns = spp
    return(w_df)
    
def max_sharpe_ratio(mean_returns, cov_matrix, risk_free_rate):
    num_assets = len(mean_returns)
    args = (mean_returns, cov_matrix, risk_free_rate)
    constraints = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
    bound = (0.0,1.0)
    bounds = tuple(bound for asset in range(num_assets))
    result = sco.minimize(neg_sharpe_ratio, num_assets*[1./num_assets,], args=args,
                        method='SLSQP', bounds=bounds, constraints=constraints)
    r2 = pd.DataFrame(result.x,index=cov_matrix.columns,columns=['Weight'])
    return r2

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
