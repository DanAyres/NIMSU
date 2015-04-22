'''
Created on 9 Feb 2015

@author: daniel
'''

import unittest
import numpy as np
from Sampling.Form_Covariance import Exponential_kernel, Calculate_Covariance_Matrix, Form_Covariance_Matrix, PosDef

from StringIO import StringIO
#import sys

class TestFormCovariance(unittest.TestCase):
    
    def testExpoKernel(self):
        
        N=3
        l=0.1
        a=np.zeros([N,N], float)
        for i in range(N):
            for j in range(i+1):
                a[i,j] = np.exp( -np.fabs(i-j)/l )
                if j!=i: a[j,i] = a[i,j]
        
        b=Exponential_kernel(N, l)
        
        self.assertTrue(np.allclose(a,b, rtol=1e-05, atol=1e-08))
        
    def testCalcCov(self):
        
        N=3
        a=np.ones([N,N], float)
        Sigma=np.zeros(N,float)
    
        for i in range(N):
            Sigma[i] = i+1
            
        b=Calculate_Covariance_Matrix(a, Sigma, N)
        for i in range(N):
            for j in range(N):
                a[i,j] = (i+1)*(j+1)
        
        self.assertTrue(np.allclose(a,b, rtol=1e-05, atol=1e-08))
    
    
    def testPosDef(self):

        output=StringIO()
        N=3
        Sigma=np.zeros(N,float)
    
        for i in range(N):
            Sigma[i] = i+1

        # Test reconstruction of covariance
        kernel_args=0.1
        Corr_Matrix=Exponential_kernel(N, kernel_args)
        Covariance_Matrix=Calculate_Covariance_Matrix(Corr_Matrix, Sigma, N)
        test_Matrix,flag=PosDef(Covariance_Matrix, N, output)
        self.assertTrue(np.allclose(Covariance_Matrix, test_Matrix, rtol=1e-05, atol=1e-08))
        self.assertEqual(flag, False, "POSDEF negative eigenvalues")
        
        # Test catch negative eigenvalues
        N=2
        Covariance_Matrix=np.zeros([N,N],float)
        Covariance_Matrix[0,0] = 1;Covariance_Matrix[1,1]=1
        Covariance_Matrix[0,1] = 2;Covariance_Matrix[1,0]=2
        test_Matrix1,flag=PosDef(Covariance_Matrix, N, output)
        test_Matrix=np.ones([N,N], float)
        test_Matrix *=1.5
        self.assertTrue(np.allclose(test_Matrix, test_Matrix1, rtol=1e-05, atol=1e-08))
        self.assertEqual(flag, True, "POSDEF NO negative eigenvalues")
        

    def testFormCovMat(self):
        
        output=StringIO()
        
        # depend==no
        N=3
        Corr_Matrix=np.zeros([N,N], float)
        test_Matrix=np.zeros([N,N], float)
        Sigma=np.zeros(N,float)
    
        for i in range(N):
            Sigma[i] = i+1
            
        depend='no'
        kernel='exponential'
        kernel_args=0.1
        
        Covariance_Matrix=Form_Covariance_Matrix(Corr_Matrix, Sigma, N, depend, kernel, kernel_args, output)
        
        for i in range(N):
            test_Matrix[i,i] = Sigma[i]**2
                
        self.assertTrue(np.allclose(Covariance_Matrix, test_Matrix, rtol=1e-05, atol=1e-08))
        
        # depend==yes
        depend='yes'
        Covariance_Matrix=Form_Covariance_Matrix(Corr_Matrix, Sigma, N, depend, kernel, kernel_args, output)
        self.assertTrue(np.allclose(Covariance_Matrix[:,0], Sigma**2, rtol=1e-05, atol=1e-08))
        
        
        # Test bad kernel
        Corr_Matrix=None
        depend='corr'
        kernel='dave'
        self.assertRaises(ValueError, Form_Covariance_Matrix, Corr_Matrix, Sigma, N, depend, kernel, kernel_args, output)
        
        
        # Test expo kernel
        test_Matrix = Exponential_kernel(N, kernel_args)
        for i in range(N):
            for j in range(N):
                test_Matrix[i,j] *= Sigma[i]*Sigma[j]
        kernel='exponential'
        Covariance_Matrix=Form_Covariance_Matrix(Corr_Matrix, Sigma, N, depend, kernel, kernel_args, output)
        self.assertTrue(np.allclose(Covariance_Matrix, test_Matrix, rtol=1e-05, atol=1e-08))
        
        
        
        
        
        
        