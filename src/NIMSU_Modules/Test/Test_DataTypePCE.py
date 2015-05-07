'''
Created on 29 Apr 2015

@author: daiel
'''
import unittest

from NIMSU_Modules.DataTypePCE import PCE_Element, IndextoBasis, PCE, PceBasis
from NIMSU_Modules.Test.DummyClasses import engine_test, sampler_test_base
from NIMSU_Modules.DataTypePCE import Hermite_Coefficients, Legendre_Coefficients
from NIMSU_Modules.DataType_Results import listData, singleData

from NIMSU_Modules.SparseGridBaseClass import sparse_grid

import numpy as np
from collections import namedtuple

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
                       scm_me_dims=None,\
                       scm_quad_tol=1E-10,\
                       scm_sg_level=[6],\
                       scm_pce=True,\
                       scm_poly_order=7,\
                       Results_listlen=3)

class Test(unittest.TestCase):

    def testPolyCoefficients(self):
        """
            Check that the Coefficients for the
             
                Hermite
                Legendre
                
            polynomials are correct.
        
        """
        
        Pmax=7
        val=0.0
        rv=np.power(val,range(Pmax+1)  )
        L=Legendre_Coefficients(Pmax)
        H=Hermite_Coefficients(Pmax)
        
        # The Hermite polynomials evaluated at x=0 correspond to the Hermite numbers
        Hermite_Numbers=np.asarray([(1-2*np.mod(i/2,2))*float(semifact(i-1)) if np.mod(i,2)==0 else 0.0 for i in range(Pmax+1) ])
        testi=np.allclose(Hermite_Numbers, H.dot(rv), rtol=1e-10, atol=1e-10)
        self.assertEqual(testi,True)

        # Legendre tests        
        val=1.0
        rv=np.power(val,range(Pmax+1)  )
        testi=np.allclose(np.ones(Pmax+1,float), L.dot(rv), rtol=1e-10, atol=1e-10)
        self.assertEqual(testi,True)
        

    def TtestPolyProd(self):
        """
            Test the PCE_Element method <PolyProd>.
            
            The integral over the polynomial product
            divided by the product of the normalisations 
            coefficients should equal unity.
            
            Test for uniform, gaussian and a mixture of the two.
        
        """
        
        
        Total_Stoch_Dim=3
        SG_LevelMax=4
        ValueType='single'
        Pmax=6
        
        
        # Uniform with Clensahw-Curtis quadrature (Tensor Product)        
        Rules=np.ones(Total_Stoch_Dim,int)
        # Create sparse grid instance
        sg=sparse_grid(Rules, \
                    Total_Stoch_Dim, \
                    SG_LevelMax,\
                    ValueType,
                    TensorProd=True)
        
        # create polynomial instance
        newPCE=PCE(Total_Stoch_Dim, 
                 Pmax,
                 ValueType,
                 Rules,
                 SG_LevelMax, 
                 Pmin=0, 
                 Dims=None, 
                 ME_Adapt_Dims=None, 
                 QuadIdx=None,
                 CalcQuadrature=True)

        for j in range(Pmax+1):
            orders=j*np.ones(Total_Stoch_Dim,int)
            val=0.0
            for i in range(sg.Num_Points):    
                
                x=sg.Points[:,i]
                val+= newPCE.Elements[0].PolyProd(orders, x, range(Total_Stoch_Dim))**2*sg.Weights[i]
            
            
            testi=np.allclose(newPCE.PDF*val, newPCE.Elements[0].PolyNorms(orders,range(Total_Stoch_Dim)), rtol=1e-10, atol=1e-10)
            self.assertEqual(testi,True)
        

        SG_LevelMax=3
        Rules=5*np.ones(Total_Stoch_Dim,int)
        # Create sparse grid instance
        sg=sparse_grid(Rules, \
                    Total_Stoch_Dim, \
                    SG_LevelMax,\
                    ValueType,
                    TensorProd=True)
        
        # create polynomial instance
        newPCE=PCE(Total_Stoch_Dim, 
                 Pmax,
                 ValueType,
                 Rules,
                 SG_LevelMax, 
                 Pmin=0, 
                 Dims=None, 
                 ME_Adapt_Dims=None, 
                 QuadIdx=None,
                 CalcQuadrature=True)

        for j in range(Pmax+1):
            orders=j*np.ones(Total_Stoch_Dim,int)
            val=0.0
            for i in range(sg.Num_Points):    
                
                x=sg.Points[:,i]
                val+= newPCE.Elements[0].PolyProd(orders, x, range(Total_Stoch_Dim))**2*sg.Weights[i]
            
            
            testi=np.allclose(newPCE.PDF*val, newPCE.Elements[0].PolyNorms(orders,range(Total_Stoch_Dim)), rtol=1e-6, atol=1e-10)
            self.assertEqual(testi,True)
            #print (newPCE.PDF*val- newPCE.Elements[0].PolyNorms(orders,range(Total_Stoch_Dim)))/newPCE.Elements[0].PolyNorms(orders,range(Total_Stoch_Dim)), testi
            
            
    def testName(self):
        
        Total_Stoch_Dim=1
        
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
        
        newPCE=PCE(control_info, sampler, engine.Nominal)                 
        
        adapt_Int_Complete=False
        while not adapt_Int_Complete:
        
            Samples = newPCE.CalculateSamples()
            
            Results=engine.interface(Samples)
            
            adapt_Int_Complete = newPCE.Update(Results)
            
            print adapt_Int_Complete
         
        # Project onto the polynomial basis   
        for ele in newPCE.Elements:
            for poly in ele.pce:
                for pos,p in enumerate(ele.SparseGrid.Index):
                    pts=ele.SparseGrid.Points[poly[2],p]
                    poly[0] += newPCE.PDF * ele.SparseGrid.Weights[pos] * ele.PolyProd(poly[1], pts, poly[2]) * ele.SparseGrid.Values[p]
                poly[0] /= ele.PolyNorms(poly[1], poly[2])
        
                print  poly[0].val, poly[0].list
            #print 
def semifact(n):
    """
        Compute the double factorial or semi-factorial: n!!
        
        Args:    
            n:    int
            
        
    """
        
    if n<= 0:
        return 1
    else:
        return n * semifact(n-2)