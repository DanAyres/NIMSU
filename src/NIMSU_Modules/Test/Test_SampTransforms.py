'''
Created on 9 Feb 2015

@author: daniel
'''

from NIMSU_Modules.SampTransforms import Transform_Statistics, FormCorrelation, LN2N_Covar, LN2N_Mean, Transform_LogNormal
import numpy as np

import unittest

class TestTransformStatistics(unittest.TestCase):
    
    def testNormal(self):
        
        NParams=4
        
        Covariance_Matrix = np.random.rand(NParams**2).reshape(NParams,NParams)
        Sigma=np.sqrt(Covariance_Matrix.diagonal())
        Means=np.random.rand(NParams)
        dist='normal'
        depend='corr'
        Corr_Matrix=FormCorrelation(Covariance_Matrix, Sigma, NParams)
        
        test_Cov, test_Corr, test_Sigma, test_Means = Transform_Statistics(Covariance_Matrix, 
                                                                           Corr_Matrix, 
                                                                           Sigma, 
                                                                           Means, 
                                                                           NParams, 
                                                                           dist, 
                                                                           depend)
        
        
        
        self.assertTrue(np.allclose(test_Cov, Covariance_Matrix, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(test_Corr, Corr_Matrix, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(test_Sigma, Sigma, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(test_Means, Means, rtol=1e-05, atol=1e-08))
    
    
    def testLN2N_Covar(self):
        NParams=3
        
        # Test covariance etc
        Covariance_Matrix=np.zeros([NParams,NParams],float)
        for i in range(NParams):
            for j in range(NParams):
                Covariance_Matrix[i,j] = (i+1)*(j+1)
            Covariance_Matrix[i,i] = np.exp(1) - 1
        
        Means=np.ones(NParams, float)
        
        test_Cov = np.zeros([NParams,NParams],float)
        for i in range(NParams):
            for j in range(NParams):
                test_Cov[i,j] = np.log( (i+1)*(j+1) + 1.0 )
            test_Cov[i,i] = 1.0    
        
        N_Cov,N_Corr,N_SD = LN2N_Covar(Means, Covariance_Matrix, NParams)
        
        self.assertTrue(np.allclose(N_Cov, test_Cov, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Corr, test_Cov, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_SD, np.ones(NParams,float), rtol=1e-05, atol=1e-08))
        
        
        # Standard deviation should equal 3
        Covariance_Matrix=np.zeros([NParams,NParams],float)
        for i in range(NParams):
            for j in range(NParams):
                Covariance_Matrix[i,j] = (i+1)*(j+1)
            Covariance_Matrix[i,i] = np.exp(9) - 1
        N_Cov,N_Corr,N_SD = LN2N_Covar(Means, Covariance_Matrix, NParams)
        
        self.assertTrue(np.allclose(N_SD, np.ones(NParams,float)*3, rtol=1e-05, atol=1e-08))

    def testLN2N_Mean(self):
        
        NParams=3

        # Test mean values
        Sigma = np.zeros(NParams,float)
        Means = np.asarray([float(i) for i in range(1,NParams+1)] )
        N_Means=LN2N_Mean(Means, Sigma)
        test_Means = np.asarray( [np.log(i) for i in range(1,NParams+1)] )
        self.assertTrue(np.allclose(N_Means, test_Means, rtol=1e-05, atol=1e-08))
        
        
        Means = np.ones(NParams,float)
        Sigma = np.asarray([float(i) for i in range(1,NParams+1)] )
        N_Means=LN2N_Mean(Means, Sigma)
        test_Means = np.asarray( [np.log( 1.0/np.sqrt(i**2 + 1.0) ) for i in range(1,NParams+1)] )
        self.assertTrue(np.allclose(N_Means, test_Means, rtol=1e-05, atol=1e-08))
        
        
    def testTransformLogNormal_corr(self):

        NParams=3
        
        depend='corr'
        Covariance_Matrix=np.zeros([NParams,NParams],float)
        for i in range(NParams):
            for j in range(NParams):
                Covariance_Matrix[i,j] = (i+1)*(j+1)
            Covariance_Matrix[i,i] = np.exp(1) - 1
        
        Means=np.ones(NParams, float)
        Sigma = np.asarray([float(i) for i in range(1,NParams+1)] )
        
        test_Cov = np.zeros([NParams,NParams],float)
        for i in range(NParams):
            for j in range(NParams):
                test_Cov[i,j] = np.log( (i+1)*(j+1) + 1.0 )
            test_Cov[i,i] = 1.0
        
        test_Means = np.asarray( [np.log( 1.0/np.sqrt(i**2 + 1.0) ) for i in range(1,NParams+1)] )
        
        
        N_Cov, N_Corr, N_Sigma, N_Means = Transform_LogNormal(Covariance_Matrix, 
                                                              Sigma, 
                                                              Means, 
                                                              NParams, 
                                                              depend)
        
        
        
        self.assertTrue(np.allclose(N_Cov, test_Cov, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Corr, test_Cov, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Sigma, np.ones(NParams,float), rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Means, test_Means, rtol=1e-05, atol=1e-08))
        

    def testTransformLogNormal_no(self):
        
        NParams=3
        Covariance_Matrix=None

        depend='no'
        Sigma=np.asarray([ np.sqrt(np.exp(1) - 1) for i in range(NParams) ])
        Means=np.ones(NParams, float)
        test_Cov=np.zeros([NParams,NParams],float)
        for i in range(NParams):
            test_Cov[i,i] = 1.0
        test_Means=np.asarray([ -0.5 for i in range(NParams) ])
        N_Cov, N_Corr, N_Sigma, N_Means = Transform_LogNormal(Covariance_Matrix, 
                                                              Sigma, 
                                                              Means, 
                                                              NParams, 
                                                              depend)
        
        self.assertEqual(N_Corr, None, "Correlation not equal to zero")
        self.assertTrue(np.allclose(N_Cov, test_Cov, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Means, test_Means, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Sigma, np.ones(NParams,float), rtol=1e-05, atol=1e-08))
     
     
    def testTransformLogNormal_yes(self):
        
        NParams=3
        Covariance_Matrix=None
           
        depend='yes'
        Sigma=np.asarray([ np.sqrt(np.exp(1) - 1) for i in range(NParams) ])
        Means=np.ones(NParams, float)
        test_Cov=np.zeros([NParams,1],float)
        for i in range(NParams):
            test_Cov[i,0] = 1.0
        test_Means=np.asarray([ -0.5 for i in range(NParams) ])
        N_Cov, N_Corr, N_Sigma, N_Means = Transform_LogNormal(Covariance_Matrix, 
                                                              Sigma, 
                                                              Means, 
                                                              NParams, 
                                                              depend)
        
        self.assertEqual(N_Corr, None, "Correlation not equal to zero")
        self.assertTrue(np.allclose(N_Cov, test_Cov, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Means, test_Means, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(N_Sigma, np.ones(NParams,float), rtol=1e-05, atol=1e-08))
        