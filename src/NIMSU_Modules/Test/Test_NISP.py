'''
Created on 6 May 2015

@author: daiel
'''
import unittest
from NIMSU_Modules.NISP import NISP
from NIMSU_Modules.Test.DummyClasses import engine_test, sampler_test_base
from collections import namedtuple
import numpy as np
from NIMSU_Modules.DataType_Results import singleData, listData

control = namedtuple('control', 'ValueType ,\
                                 use_scm ,\
                                 scm_quad_adapt ,\
                                 scm_use_me,\
                                 scm_me_dims,\
                                 scm_quad_tol,\
                                 scm_sg_level,\
                                 scm_pce ,\
                                 scm_poly_order ,\
                                 Results_listlen')
control_info = control(ValueType = 'list'      ,\
                       use_scm=False,\
                       scm_quad_adapt=False,\
                       scm_use_me=False,\
                       scm_me_dims=[0],\
                       scm_quad_tol=1E-10,\
                       scm_sg_level=[6],\
                       scm_pce=True,\
                       scm_poly_order=4,\
                       Results_listlen=4)

class Test(unittest.TestCase):


    def testName(self):
        
        Total_Stoch_Dim=2
        
        # Setup the sampling object        
        mean = np.zeros(Total_Stoch_Dim,float)
        sd = (1.0/np.sqrt(3.0))*np.ones(Total_Stoch_Dim,float)
        rules=np.ones(Total_Stoch_Dim,int)
        sampler = sampler_test_base(Total_Stoch_Dim, mean, sd, rules)
        
        
        
        # Setup the engine object
        val=sum(mean**2)        
        if control_info.ValueType=='single':
            nominal = singleData(val)
        elif control_info.ValueType=='list':
            nominal = listData(val,val*np.ones(control_info.Results_listlen,float))
        engine = engine_test(control_info, nominal)
        
        
        NISP(control_info, sampler, engine)
        


