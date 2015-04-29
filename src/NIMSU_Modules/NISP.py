'''
Created on 24 Apr 2015

@author: daiel
'''

import numpy as np
from SparseGrids.SparseGridBaseClass import sparse_grid, CalcSparseSet


def NISP(control_info, sampler, engine):

    # identify which dimensions are to be split.
    
    # split elements
    
    # list of elements
    # each element needs to know domain, pdf, etc.
    
    

    me_dims=[1]
    

    multi=[False for i in range(sampler.Total_stoch_dim+1)]
    multi[2]=True
    
    print multi

    rules=np.ones(sampler.Total_stoch_dim,int)
    sg=sparse_grid(rules, \
                sampler.Total_stoch_dim, \
                2, \
                engine.Nominal, \
                Adaptive=False, \
                AdaptTol=control_info.hdmr_quad_tol, \
                MaxVals=40, \
                ValueType=control_info.ValueType)

    #for i in range(sg.Num_Points):
    #    print sg.Points[:,i] 
    pass    
