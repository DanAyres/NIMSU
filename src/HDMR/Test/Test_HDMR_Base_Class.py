'''
Created on 14 Apr 2015

@author: daiel
'''
import unittest
from HDMR.HDMR_Base_Class import HDMR_Base
import numpy as np
from collections import namedtuple
from DataTypes.Results import singleData,listData
from SparseGrids.LibSG import sparse_grid

class sampler_test_base():
    
    def __init__(self, Total_stoch_dim, mean, sd):
        
        self.Total_stoch_dim=Total_stoch_dim
        self.rules = np.ones(self.Total_stoch_dim,int)
        self.PDF = 0.5*np.ones(self.Total_stoch_dim,float)
        
        self.mean = mean
        self.sd = sd
        
    def sample(self,rvs):
        """
            assume uniformly distributed, uncorrelated uncertainties.
        """
        samples={}
        
        for pos in range(self.Total_stoch_dim):
            key='val'+str(pos)
            samples[key] = self.mean[pos] + np.sqrt(3.0) * self.sd[pos] * rvs[pos] 
        
        return samples
    

class engine_test():
    
    def __init__(self, control_info, Nominal):
        """
        
        """
        
        self.control_info = control_info
        self.Nominal = Nominal
        
        pass

    def interface(self, Samples):
        
        Results={}
        if self.control_info.ValueType == 'single':
            for key in Samples:
                val=0.0
                for name in Samples[key]:
                    val += Samples[key][name]**2
                Results[key] = singleData(val)
        elif self.control_info.ValueType=='list':
            for key in Samples:
                val=0.0
                for name in Samples[key]:
                    val += Samples[key][name]**2
                Results[key] = listData(val,listt=val*np.ones(self.control_info.Results_listlen,float))
        
        return Results 

control = namedtuple('control', 'ValueType Results_listlen hdmr_quad_adapt hdmr_quad_tol hdmr_sg_level hdmr_metric')
control_info = control(ValueType = 'single'      ,\
                        Results_listlen = 2      ,\
                        hdmr_quad_tol = 1.0E-7   ,\
                        hdmr_metric = 1.0E-4     ,\
                        hdmr_sg_level = [2,2,2,2],\
                        hdmr_quad_adapt = False    
                        )
 
class Test(unittest.TestCase):


    def testsubsetContributions(self):
        
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
        newHDMR.data_struct[0]=[ [ [0], 0.0, 1.0*a, 1.0*a, None, False ] ,\
                                 [ [1], 0.0, 2.0*a, 4.0*a, None, False ] ,\
                                 [ [2], 0.0, 3.0*a, 9.0*a, None, False ]  ]
        newHDMR.data_struct.append([ [ [0,1], 0.0, 4.0*a, 16.0*a, None, False ] ,\
                                     [ [0,2], 0.0, 5.0*a, 25.0*a, None, False ] ,\
                                     [ [1,2], 0.0, 6.0*a, 36.0*a, None, False ]  ] )

        newHDMR.data_struct.append([ [ [0,1,2], 0.0, 7.0*a, 49.0*a, None, False ] ])
        
        # check the one dimensional sets have zero subset contributions
        for i in range(3):
            mean,var= newHDMR.subsetContributions([i])
            self.assertEqual(mean.val, 0.0); self.assertEqual(var.val, 0.0)
            
        # check the two dimensional sets
        mean,var= newHDMR.subsetContributions([0,1])
        self.assertEqual(mean.val, 3.0); self.assertEqual(var.val, 5.0)
        mean,var= newHDMR.subsetContributions([0,2])
        self.assertEqual(mean.val, 4.0); self.assertEqual(var.val, 10.0)
        mean,var= newHDMR.subsetContributions([1,2])
        self.assertEqual(mean.val, 5.0); self.assertEqual(var.val, 13.0)
        
        # check the three dimensional set
        mean,var= newHDMR.subsetContributions([0,1,2])
        self.assertEqual(mean.val, 21.0); self.assertEqual(var.val, 91.0)
        
        # remove some sets
        del newHDMR.data_struct[1][0] 
        mean,var= newHDMR.subsetContributions([0,1,2])
        self.assertEqual(mean.val, 17.0); self.assertEqual(var.val, 75.0)
        
    def testHDMR_functions(self):
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
        
    def testcreateSamples(self):
        
        Total_stoch_dim=4
        
        
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
        

    def testfindPreviousPoints(self):
        Total_stoch_dim=3
        quad_level=2
        
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
        

        # Init the sparse grid objects
        rules = [1]
        for data in newHDMR.data_struct[0]:  
            
            data[4] = sparse_grid(rules, \
                                        1, \
                                        quad_level, \
                                        newHDMR.R0, \
                                        Adaptive=False, \
                                        AdaptTol=newHDMR.quadTol, \
                                        MaxVals=40, \
                                        ValueType=newHDMR.ValueType)    
    
            data[4].Values=[]
            for p in range(data[4].Num_Points):
                data[4].Values.append( a * float(p) )
            
        newHDMR.data_struct.append([ [ [0,1], 0.0, 4.0*a, 16.0*a, None, False ] ] )
#                                      [ [0,2], 0.0, 5.0*a, 25.0*a, None, False ] ,\
#                                      [ [1,2], 0.0, 6.0*a, 36.0*a, None, False ]  ] )

        rules=[1,1]
        for data in newHDMR.data_struct[1]:  
            data[4] = sparse_grid(rules, \
                                        2, \
                                        quad_level, \
                                        newHDMR.R0, \
                                        Adaptive=False, \
                                        AdaptTol=newHDMR.quadTol, \
                                        MaxVals=40, \
                                        ValueType=newHDMR.ValueType)    

        dims=newHDMR.data_struct[1][0][0]
        rvs=np.zeros(Total_stoch_dim,float)
        
        for p in range(newHDMR.data_struct[1][0][4].Num_Points):
            
            points= newHDMR.data_struct[1][0][4].Points[:,p]
            if all(points==0.0): continue
            
            try:
                print newHDMR.findPreviousPoints(dims, points, sampler, rvs).val, points
            except AttributeError:
                print newHDMR.findPreviousPoints(dims, points, sampler, rvs), points
            

        #
        # Test the 3 dimensional terms
        #

        newHDMR.data_struct[1] = [ [ [0,1], 0.0, 4.0*a, 16.0*a, None, False ]  ,\
                                      [ [0,2], 0.0, 5.0*a, 25.0*a, None, False ] ,\
                                      [ [1,2], 0.0, 6.0*a, 36.0*a, None, False ]  ] 
        
            
        rules=[1,1]
        for data in newHDMR.data_struct[1]:  
            data[4] = sparse_grid(rules, \
                                        2, \
                                        quad_level, \
                                        newHDMR.R0, \
                                        Adaptive=False, \
                                        AdaptTol=newHDMR.quadTol, \
                                        MaxVals=40, \
                                        ValueType=newHDMR.ValueType)
            data[4].Values=[]
            for p in range(data[4].Num_Points):
                data[4].Values.append( a * float(p) )
            
        newHDMR.data_struct.append( [ [ [0,1,2], 0.0, 4.0*a, a, None, False ] ])

        rules=[1,1,1]
        for data in newHDMR.data_struct[2]:  
            data[4] = sparse_grid(rules, \
                                        3, \
                                        quad_level, \
                                        newHDMR.R0, \
                                        Adaptive=False, \
                                        AdaptTol=newHDMR.quadTol, \
                                        MaxVals=40, \
                                        ValueType=newHDMR.ValueType)
            
            
        print '\n\n'
        dims=newHDMR.data_struct[2][0][0]
        for p in range(newHDMR.data_struct[2][0][4].Num_Points):
            
            points= newHDMR.data_struct[2][0][4].Points[:,p]
            if all(points==0.0): continue
            
            
            try:
                print newHDMR.findPreviousPoints(dims, points, sampler, rvs).val, points
            except AttributeError:
                print newHDMR.findPreviousPoints(dims, points, sampler, rvs), points