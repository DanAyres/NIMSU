'''
Created on 10 Feb 2015

@author: daniel
'''

import unittest
from Sampling.Decompose import Karhunen_Loeve, Cholesky, Cumulative_Eigenvalues, Decompose_Covariance_Matrix
from Sampling.Form_Covariance import Exponential_kernel
from StringIO import StringIO
import numpy as np

class TestDecomposeCovariance(unittest.TestCase):
    
    def testDecompose_Covariance_Matrix(self):
        
        NParams = 3
        Cov_Mat = np.zeros([NParams,NParams],float)
        depend = 'yes'
        decomp = 'chol'
        kl_theta = 1.0
        output=StringIO()
        
        self.assertRaises(ValueError, Decompose_Covariance_Matrix, Cov_Mat, depend, decomp, NParams, kl_theta, output)
        
        Cov_Mat = 4.0*np.ones([NParams,1],float)
        Sampling_Matrix, stoch_dim=Decompose_Covariance_Matrix(Cov_Mat, depend, decomp, NParams, kl_theta, output)
        self.assertEqual(stoch_dim, 1) 
        self.assertTrue(np.allclose(Sampling_Matrix, np.ones([NParams,1],float)*2, rtol=1e-05, atol=1e-08))     
        
    
    def testCumulative_Eigenvalues(self):
        
        N=100
        eig_vals = np.linspace(0, 1, N, endpoint=True)
        
        tol=0.98
        stoch_dim = Cumulative_Eigenvalues(eig_vals, tol)
        self.assertEqual(stoch_dim, 99)
        
        tol=0.90
        stoch_dim = Cumulative_Eigenvalues(eig_vals, tol)
        self.assertEqual(stoch_dim, 95)
        
    
    def testKarhunen_Loeve(self):
        
        NParams = 5
        cor_len = 1.0
        
        kl_theta = 1.0 
        
        Cov_Mat = Exponential_kernel(NParams, cor_len)
        
        Sampling_Matrix, stoch_dim = Karhunen_Loeve(Cov_Mat, NParams, kl_theta)
        
        self.assertTrue(np.allclose(Sampling_Matrix.dot(Sampling_Matrix.T), Cov_Mat,  rtol=1e-05, atol=1e-08), "Error with KL")
        self.assertEqual(stoch_dim, NParams, "Wrong stoch_dim")
        
        
        
    def testCholesky(self):
        
        NParams=5
        cor_len = 1.0
        
        Cov_Mat = Exponential_kernel(NParams, cor_len)
        Sampling_Matrix = Cholesky(Cov_Mat)
        
        self.assertTrue(np.allclose(Sampling_Matrix.dot(Sampling_Matrix.T), Cov_Mat), "Error with Cholesky")
        
        
        # From the wiki page
        NParams = 3
        Cov_Mat = np.zeros([NParams,NParams],float)
        Cov_Mat[0,0] = 4;   Cov_Mat[0,1] = 12; Cov_Mat[0,2] = 16
        Cov_Mat[1,0] = 12;  Cov_Mat[1,1] = 37; Cov_Mat[1,2] = 43
        Cov_Mat[2,0] = -16; Cov_Mat[2,1] = -43; Cov_Mat[2,2] = 98
        
        test_Mat = np.zeros([NParams,NParams],float)
        test_Mat[0,0] = 2.0
        test_Mat[1,0] = 6.0;  test_Mat[1,1] = 1.0
        test_Mat[2,0] = -8.0; test_Mat[2,1] = 5.0; test_Mat[2,2] = 3.0
        self.assertTrue(np.allclose(Cholesky(Cov_Mat), test_Mat))
