'''
Created on 4 Feb 2015

@author: daniel
'''
from UserInterface.interface import ReadControlFile, ReadParamFile
from Sampling.samples import SampleUncertainties
from SparseGrids.LibSG import sparse_grid, CalcSparseSet
from HDMR.HDMR_Base_Class import HDMR_Base
import os
import sys
import numpy as np
from collections import namedtuple

from HDMR.Test.Test_HDMR_Base_Class import engine_test, control_info, sampler_test_base
from DataTypes.Results import listData, singleData


if __name__ == "__main__":
    
#     filename="../1D.control"
#     
#     control_info = ReadControlFile(open(filename,'r'), sys.stdout)
#     
#     input_params=ReadParamFile(open("../"+control_info.param_infile,'r'), sys.stdout)
#     
#     sampler = SampleUncertainties(input_params, control_info, sys.stdout)
#     
        
    Total_stoch_dim=3
        
    # Setup the sampling object        
    mean = np.zeros(Total_stoch_dim,float)
    sd = (1.0/np.sqrt(3.0))*np.ones(Total_stoch_dim,float)
    sampler = sampler_test_base(Total_stoch_dim, mean, sd)
        
        
    # Setup the engine object
    val=sum(mean**2)        
    if control_info.ValueType=='single':
        nominal = singleData(val)
    elif control_info.ValueType=='list':
        nominal = listData(val,val*np.ones(control_info.Results_listlen,float))
    engine = engine_test(control_info, nominal)
        
        
    # Create new HDMR object
    newHDMR = HDMR_Base(control_info, engine, sampler) 
    
    # Create some test data
    if control_info.ValueType=='single':
        a = singleData(1.0)
    elif control_info.ValueType=='list':
        a = listData(1.0,np.ones(control_info.Results_listlen,float))
    newHDMR.data_struct.append([ [ [0,1], 0.0, 4.0*a, 16.0*a, None, False ] ,\
                                 [ [0,2], 0.0, 5.0*a, 25.0*a, None, False ] ,\
                                 [ [1,2], 0.0, 6.0*a, 36.0*a, None, False ]  ] )
    for data in newHDMR.data_struct[0]:
        if data[4] == None:
            rules = [1]
            data[4] = sparse_grid(rules, \
                                1, \
                                2, \
                                newHDMR.R0, \
                                Adaptive=False, \
                                AdaptTol=newHDMR.quadTol, \
                                MaxVals=40, \
                                ValueType=newHDMR.ValueType)

            data[4].Values=[]
            # set the quadrature values to x**2
            for p in range(data[4].Num_Points):
                data[4].Values.append( a * sum(data[4].Points[:,p]**2) )
    
    dimension=2
    Samples=newHDMR.createSamples(dimension, sampler, newHDMR.data_struct[dimension-1])
    print Samples
    
    
