'''
Created on 6 May 2015

@author: daiel
'''
import numpy as np
from NIMSU_Modules.DataType_Results import singleData, listData
from collections import namedtuple

class sampler_test_base():
    
    def __init__(self, Total_stoch_dim, mean, sd, rules):
        
        self.Total_stoch_dim=Total_stoch_dim
        self.rules = rules
        self.PDF = []
        
        self.mean = mean
        self.sd = sd
        
        for r in self.rules:
            
            if r in [1,2,3,4]:
                self.PDF.append(0.5)
            elif r in [5,6,10]:
                self.PDF.append(np.sqrt(np.pi))
            elif r in [7,8]:
                print 'Not implemented'
                raise
            elif r==9:
                print 'Not implemented'
                raise    
            else:
                raise ValueError("Quadrature rule not recognised")     
        
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
                lst=np.zeros(self.control_info.Results_listlen,float)
                for name in Samples[key]:
                    val += Samples[key][name]**2
                    lst+=np.power(Samples[key][name], range(1,self.control_info.Results_listlen+1) )
                Results[key] = listData(val,listt=lst)
#                 Results[key] = listData(val,listt=val*np.ones(self.control_info.Results_listlen,float))
        
        return Results 

control = namedtuple('control', 'ValueType Results_listlen hdmr_quad_adapt hdmr_quad_tol hdmr_sg_level hdmr_metric')
control_info = control(ValueType = 'single'      ,\
                        Results_listlen = 2      ,\
                        hdmr_quad_tol = 1.0E-7   ,\
                        hdmr_metric = 1.0E-4     ,\
                        hdmr_sg_level = [2,2,2,2],\
                        hdmr_quad_adapt = False    
                        )
