'''
Created on 9 Feb 2015

@author: daniel
'''

import numpy as np

def Assign_Rules_PDFS(dist, stoch_dim):
    
    # return rules, PDF
    if dist=='normal':
        return [10 for dummy in range(stoch_dim)], [1.0/np.sqrt(np.pi) for dummy in range(stoch_dim)]
    elif dist=='lognormal':
        return [10 for dummy in range(stoch_dim)], [1.0/np.sqrt(np.pi) for dummy in range(stoch_dim)]
    elif dist=='uniform':
        return [1 for dummy in range(stoch_dim)], [0.5 for dummy in range(stoch_dim)]
    else:
        print 'Distribution not recognised'
        raise ValueError

def Assign_Sampling_Routine(dist, depend):
    
    if dist=='normal':
        return Norm_Sample
    elif dist=='lognormal':
        return LogNorm_Sample
    elif dist=='uniform':
        return Uni_Sample
    else:
        print 'Distribution not recognised'
        raise ValueError
            
            
            
def LogNorm_Sample(means, C, rvs):
    return np.exp(means + C.dot(rvs))

def Norm_Sample(means, C, rvs):
    return means + C.dot(rvs)
    
def Uni_Sample(means, C, rvs):
    return means + np.sqrt(3.0)*C.dot(rvs)         