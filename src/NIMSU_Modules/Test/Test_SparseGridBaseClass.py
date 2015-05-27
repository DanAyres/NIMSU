import unittest
from NIMSU_Modules.SparseGridBaseClass import sparse_grid, Compute_Grid, CalcSparseSet
from FORTRAN_Modules.sandiainterface import reduced_sg as Sandia
import numpy as np
import sys
from NIMSU_Modules.DataType_Results import singleData, listData, hdf5Data

from collections import namedtuple
from operator import mul


class TestLibSG(unittest.TestCase):


    def TtestErrorMeasure(self):
        
        """
         Test the sparse_grid.ErrorMeasure method 
        """
        
        num_dim=1
        LevelMax=4
        rules=3*np.ones(num_dim,int)
        exponent=2
        nom=intfunc_uni(np.zeros(num_dim,float),exponent)
        
        vt='list'
        listlen=3
        new_int=[listData(1.0,np.zeros(listlen,float)),listData(3.0,np.zeros(listlen,float))]
        old_int=[listData(0.5,np.zeros(listlen,float)),listData(1.0,np.zeros(listlen,float))]
        
        
        
        wmix = 0.0
        Nominal=listData(nom,listt=nom*np.ones(listlen,float))
        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix, Listlen=listlen)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 2.0)        

        wmix = 1.0
        Nominal=listData(nom,listt=nom*np.ones(listlen,float))
        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix, Listlen=listlen)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 1.0)


        vt='single'
        wmix = 0.0
        Nominal=singleData(nom)
        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 2.0)        

        wmix = 1.0
        Nominal=listData(nom)
        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 1.0)


    def TtestSparseGrid_uni(self):
        """
            test the sparse grid construction via the evaluation of some test integrals.
            
            A series of integrals are performed on [-1,1] using a constant weight = 1.
              
             The integrand is
              
                 sum_{i=1}^N x_i^k prod_{j=1}^N dx_j
              
             with solution
              
                I  = 2^N * N/(k+1)    for k even
                I  = 0                for k odd
            
            
        """
        
        LevelMax=4
        r_list=[1,2,3,4]
        exponent=2
        
        for num_dim in range(1,5):
            for r in r_list:
                #
                # Check integration works with list data
                #    
                vt='list'
                listlen=3
                rules=r*np.ones(num_dim,int)
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Listlen=listlen )
                
                Idx=range(newgrid.Num_Points)
                newgrid.Values=[listData(0.0,np.zeros(listlen,float)) for i in range(newgrid.Num_Points)]
                
                for i in Idx:
                    val=intfunc_uni(newgrid.Points[:,i], exponent)
                    newgrid.Values[i].val = val
                    newgrid.Values[i].list = val*np.ones(3,float)
        
                integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                
                exact=intfunc_uni_analytic(exponent, num_dim)
                self.assertTrue(np.allclose(integral[0].val, exact, rtol=1e-10, atol=1e-10))
                self.assertTrue(np.allclose(integral[0].list, exact*np.ones(listlen,float), rtol=1e-10, atol=1e-10))
                
                #
                # Check integration works with single data
                #
                vt='single'
                exponent=2
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt)
                
                Idx=range(newgrid.Num_Points)
                newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
                for i in Idx:
                    newgrid.Values[i].val = intfunc_uni(newgrid.Points[:,i], exponent)
                    
                integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                exact=intfunc_uni_analytic(exponent, num_dim)
                self.assertTrue(np.allclose(integral[0].val, exact, rtol=1e-10, atol=1e-10))
     
    def TtestSparseGrid_Laguerre(self):
        """
            Test the construction and application of Gauss-Laguerre quadratures
            via the evaluation of test integrals.
            
            The weight function and support for the Laguerre polynomials
            are
            
                supp. [0,inf]
                w(x) = exp(-x)
            
            The weight and support for the Generalize Lag polys are
            
                supp. [0,inf]
                w(x) = x^a exp(-x)    a > -1
            
            The following test integrals have been defined for Gauss-Laguerre integrals:
            
            1)  int exp(-px) = 1/p  [0,inf]
                -->  int exp(-x) * exp( x[1-p] ) = 1/p
                
            2) int exp(-p^2 x^2) = sqrt(pi) / (2p)
               -->   int exp(-x) * exp( x[1-p^2 x]) = sqrt(pi) / (2p)
            
             
            Integrals for the Generalized GL integrals:
             
            3) int (1/sqrt(x)) * exp (-px) = sqrt(pi/p)
        
        
            4) int x^{n-1} exp (-px) = (1/p^n) * G(n)
            
            G: gamma function
        """
        
        from scipy.special import gamma
        
        vt='single'
        
        LevelMax=8
        order_max = 5
        num_dim=2
        
        # The basic Laguerre quadratures
        q=7
        for num_dim in range(1,3):
            for order in range(1,order_max):
                Nominal = singleData( np.exp(0.0) )
                        
                rules=q*np.ones(num_dim,int)
                                           
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal, TensorProd=False)
                         
                integral=0.0
                for p in range(newgrid.Num_Points):
                    temp = Laguerre_test(newgrid.Points[:,p], order)
                    integral += newgrid.Weights[p] * temp
                          
                int1=(1.0/order)**num_dim
                int2=(np.sqrt(np.pi)/(2.0*order))**num_dim
                self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))
                self.assertTrue(np.allclose(integral[1], int2, rtol=1e-4, atol=1e-3))   

        
            # The generalised Laguerre quadratures
            q=8
            Dist_args=[[-0.5] for dummy in range(num_dim)]  # alpha = -0.5
            for order in range(1,order_max):
                        
                rules=q*np.ones(num_dim,int)
                                           
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args)
                         
                integral=0.0
                for p in range(newgrid.Num_Points):
                    temp = Laguerre_test(newgrid.Points[:,p], order)
                    integral += newgrid.Weights[p] * temp
                 
                int1=np.sqrt(np.pi/order)**num_dim
                self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))

            for n in range(1,5):
                Dist_args=[[n-1] for dummy in range(num_dim)]  # alpha = -0.5
                for order in range(1,order_max):
                            
                    rules=q*np.ones(num_dim,int)
                                               
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args)
                    
                    
                    integral=0.0
                    for p in range(newgrid.Num_Points):
                        temp = Laguerre_test(newgrid.Points[:,p], order)
                        integral += newgrid.Weights[p] * temp
                     
                    int1 = ((1.0/order**n) * gamma(n))**num_dim
                    self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))

    def TtestSparseGrid_gauss(self):
        """
             Tests the ComputeGrid method which, given an index set, computes
             sparse grid weights and points.
             
             All of the Gaussian weight rules are tested:
             
             rule 5: Gauss Hermite
             rule 6: Generalized Gauss Hermite
             rule 10: Hermite Genz-Keister
             
             
             The integrand is
              
                 sum_{i=1}^N x_i^k prod_{j=1}^N dx_j
              
             with solution
              
                I  = 2^N * N/(k+1)    for k even
                I  = 0                for k odd
        """
  
        q_rules=[5,6,10]
        LevelMax=4
        vt='single'
  
        for q in q_rules:
            for exponent in range(1,2):
                for num_dim in range(1,7):
                    
                    if q==10: LevelMax=4
  
                    rules=q*np.ones(num_dim,int)
                    
                      
                    #exponent=0
                      
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt)
                
                    Idx=range(newgrid.Num_Points)
                    newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
                    for i in Idx:
                        newgrid.Values[i].val = intfunc_gauss(newgrid.Points[:,i], exponent)
                        
                    integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                    self.assertTrue(np.allclose(integral[0].val, num_dim*intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))
                      
                    newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
                    for i in Idx:
                        newgrid.Values[i].val = intfunc_gauss2(newgrid.Points[:,i],exponent)
                     
                    integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                    self.assertTrue(np.allclose(integral[0].val, intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))

    def TtestSparseGrid_Jacobi(self):
        """
        Test the Gauss-Jacobi quadrature rules.
        
        The weight function:
        
            w(x,a,b) = (1-x)^a * (1+x)^b
        
        
        
        """
        from scipy.special import gamma

        vt='single'
        LevelMax=5
        
        alpha=[i for i in range(1,3)]
        beta=alpha

        q=9
        for num_dim in range(1,5):
            for a in alpha:
                for b in beta:
                    if a==b: continue
                    Dist_args=[[b,a] for i in range(num_dim)]  # alpha = -0.5
                            
                    rules=q*np.ones(num_dim,int)
                    
                    pdf = (np.power(2,a+b+1)*gamma(a+1)*gamma(b+1) / gamma(a+b+2))**num_dim
                                               
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args)
                    self.assertTrue(np.allclose(sum(newgrid.Weights), pdf , rtol=1e-6, atol=1e-6))
                             
                    integral=0.0
                    for p in range(newgrid.Num_Points):
                        x=newgrid.Points[:,p]
                        temp = Jacobi_test(x)
                        integral += newgrid.Weights[p] * temp
                     
                    mean=integral[1]/pdf
                    var= integral[2]/pdf - mean**2
                    m2=(a-b)/float(a+b+2)
                    uv = 4*(a+1)*(b+1)/float(( (a+b+2)**2 * (3+a+b) )) + m2**2
                    v2= uv**num_dim - np.power(m2,2*num_dim)
                    #print "%5i %5i %5i %15.5e %15.5e %10.5f %10.5f %10.5f" %(a,b, num_dim, mean, m2**num_dim, mean/m2**num_dim, sum(newgrid.Weights), pdf)
                    self.assertTrue(np.allclose(m2**num_dim, mean , rtol=1e-6, atol=1e-6))
                    
                    self.assertTrue(np.allclose(v2, var , rtol=1e-6, atol=1e-6))

    def TtestSparseGrid_mixed(self):
        from scipy.special import gamma
        """
            Test the sparse grid construction of mixed univariate rules via
            the integration of test functions.
            
            
            The integrand is a follows
            
            [uniform, gauss, Laguerre, Jacobi]
            
            [x1,x2,x3,x4]
            
            f(x1,x2,x3,x4) = 
        
        """
        

        rules=[1,5,8,9]
        num_dim=4
        LevelMax=6
        vt='single'

        uni_pdf = 2.0
        gauss_pdf = np.sqrt(np.pi)
        alpha=[i for i in range(1,5)]
        beta=alpha
        
        for n in range(1,5):
            for order in range(1,2):
                for a in alpha:
                    for b in beta:
                
                        Dist_args=[None, None, [n-1] , [b,a] ]
                        
                        beta_pdf = np.power(2,a+b+1)*gamma(a+1)*gamma(b+1) / gamma(a+b+2)
                        lag_pdf = gamma(n)
                        
                        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args)
                        
                        self.assertTrue(np.allclose(sum(newgrid.Weights), uni_pdf * gauss_pdf * beta_pdf * lag_pdf, rtol=1e-6, atol=1e-6))
                        
                        
                        integral=0.0
                        for p in range(newgrid.Num_Points):
                            temp = intfunc_mixed(newgrid.Points[:,p], 1.0)
                            integral += newgrid.Weights[p] * temp
                         
                        int1=(2.0/3.0) * np.sqrt(np.pi) * (1.0/order**n) * gamma(n) * (a-b)/float(a+b+2)
                        self.assertTrue(np.allclose(integral/beta_pdf, int1, rtol=1e-3, atol=1e-3))



    def TtestAdapt_gauss(self):
        """
             Test the adaptive sparse grid method.
             
             All of the Gaussian weight rules are tested:
             
             rule 5: Gauss Hermite
             rule 6: Generalized Gauss Hermite
             rule 10: Hermite Genz-Keister
             
             
             The integrand is
              
                 sum_{i=1}^N x_i^k prod_{j=1}^N dx_j
              
             with solution
              
                I  = 2^N * N/(k+1)    for k even
                I  = 0                for k odd
        """
  
        q_rules=[5,6,10]
        LevelMax=4

        # Test for a single value
        vt='single'
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,7):
                    
                    Nominal = singleData(intfunc_gauss(np.zeros(num_dim,float),exponent))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=100)
                                
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            values.append( singleData( intfunc_gauss(newgrid.Points[:,p],exponent) )  )
                                
                        adapt_flag=newgrid.Adapt(values)
                    
                    integral = newgrid.Current_Integral[0].val 
                    self.assertTrue(np.allclose(integral, num_dim*intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))

        # Test for a list of values
        vt='list'
        listlen=2
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,7):
                    
                    val=intfunc_gauss(np.zeros(num_dim,float),exponent)
                    Nominal = listData(val,val*np.ones(listlen,float))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=100, Listlen=listlen)
                                
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            val=intfunc_gauss(newgrid.Points[:,p],exponent)
                            values.append( listData(val,val*np.ones(listlen,float)) )
                                
                        adapt_flag=newgrid.Adapt(values)
                    
                    integral = newgrid.Current_Integral[0].val 
                    self.assertTrue(np.allclose(integral, num_dim*intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))

    def TtestAdapt_uni(self):

  
        q_rules=[1,2,3,4]
        LevelMax=6

        # Test for a single value
        vt='single'
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,5):
                    
                    Nominal = singleData(intfunc_uni(np.zeros(num_dim,float),exponent))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=1000)
                                
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            values.append( singleData( intfunc_uni(newgrid.Points[:,p],exponent) )  )
                                
                        adapt_flag=newgrid.Adapt(values)
                    
                    integral = newgrid.Current_Integral[0].val 
                    self.assertTrue(np.allclose(integral, intfunc_uni_analytic(exponent, num_dim), rtol=1e-6, atol=1e-6))

        # Test for a list of values
        vt='list'
        listlen=2
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,5):
                    
                    val=intfunc_uni(np.zeros(num_dim,float),exponent)
                    Nominal = listData(val,val*np.ones(listlen,float))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=1000, Listlen=listlen)
                                
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            val=intfunc_uni(newgrid.Points[:,p],exponent)
                            values.append( listData(val,val*np.ones(listlen,float)) )
                                
                        adapt_flag=newgrid.Adapt(values)
                    
                    integral = newgrid.Current_Integral[0].val 
                    self.assertTrue(np.allclose(integral, intfunc_uni_analytic(exponent, num_dim), rtol=1e-6, atol=1e-6))

    def TtestAdapt_Jacobi(self):

        from scipy.special import gamma
  
        print '\n\n'
        q_rules=[9]
        LevelMax=6
        alpha=[i for i in range(1,3) ]
        beta=alpha
        
        # Test for a single value
        vt='single'
        for q in q_rules:
            for num_dim in range(1,8):
                for a in alpha:
                    for b in beta:
                        if a==b: continue
                        
                        Dist_args=[[b,a] for i in range(num_dim)]  
                        pdf = (np.power(2,a+b+1)*gamma(a+1)*gamma(b+1) / gamma(a+b+2))**num_dim
                        Nominal = singleData( reduce(mul,np.zeros(num_dim,float) ) )
                        rules=q*np.ones(num_dim,int)
                                  
                        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,  
                                              Adaptive=True, AdaptTol=1.0E-16, 
                                              MaxVals=1000,Distribution_Args=Dist_args)
                                
                        adapt_flag=True
                        while adapt_flag:
                            values=[]
                            for p in newgrid.ComputeIndex:
                                val=reduce(mul,newgrid.Points[:,p])
                                values.append( singleData( val )  )
                                    
                            adapt_flag=newgrid.Adapt(values)
                        
                        integral = newgrid.Current_Integral[0].val
                        mean=integral/pdf 
                        m2=((a-b)/float(a+b+2))**num_dim
                        print "%5i %5i %5i %15.5e %15.5e %10.5f %10.5f %10.5f" %(a,b, num_dim, mean, m2, mean/m2, sum(newgrid.Weights), pdf)
                        #self.assertTrue(np.allclose(mean,m2, rtol=1e-6, atol=1e-6))
    def TtestAdapt_Laguerre(self):

        from scipy.special import gamma
        
        vt='single'
        
        LevelMax=8
        order_max = 5
        num_dim=2
        
        # The basic Laguerre quadratures
        q=7
        for num_dim in range(1,2):
            print
            for order in range(1,order_max):
                Nominal = singleData( Laguerre_test(np.zeros(num_dim,float),order)[1] )
                        
                rules=q*np.ones(num_dim,int)
                                           
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal,
                                    Adaptive=True, AdaptTol=1.0E-16, 
                                    MaxVals=1000)
                         
                adapt_flag=True
                while adapt_flag:
                    values=[]
                    for p in newgrid.ComputeIndex:
                        val=Laguerre_test(newgrid.Points[:,p],order)[1]
                        values.append( singleData( val )  )
                    adapt_flag=newgrid.Adapt(values)
                integral = newgrid.Current_Integral[0].val
                          
                int1=(1.0/order)**num_dim
                int2=(np.sqrt(np.pi)/(2.0*order))**num_dim
                #self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))
                self.assertTrue(np.allclose(integral, int2, rtol=1e-4, atol=1e-3))   

            # The generalised Laguerre quadratures
            q=8
            Dist_args=[[-0.5] for dummy in range(num_dim)]  # alpha = -0.5
            for order in range(1,order_max):
                        
                rules=q*np.ones(num_dim,int)
                Nominal = singleData( Laguerre_test(np.zeros(num_dim,float),order)[0] )
                                           
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal, Distribution_Args=Dist_args,
                                      Adaptive=True, AdaptTol=1.0E-16,MaxVals=1000)
                         
                adapt_flag=True
                while adapt_flag:
                    values=[]
                    for p in newgrid.ComputeIndex:
                        val=Laguerre_test(newgrid.Points[:,p],order)[0]
                        values.append( singleData( val )  )
                    adapt_flag=newgrid.Adapt(values)
                integral = newgrid.Current_Integral[0].val
                 
                int1=np.sqrt(np.pi/order)**num_dim
                self.assertTrue(np.allclose(integral, int1, rtol=1e-6, atol=1e-6))

            for n in range(1,5):
                Dist_args=[[n-1] for dummy in range(num_dim)]  # alpha = -0.5
                for order in range(1,order_max):
                            
                    rules=q*np.ones(num_dim,int)
                    Nominal = singleData( Laguerre_test(np.zeros(num_dim,float),order)[0] )
                                               
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal, Distribution_Args=Dist_args,
                                      Adaptive=True, AdaptTol=1.0E-5,MaxVals=1000)
                    
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            val=Laguerre_test(newgrid.Points[:,p],order)[0]
                            values.append( singleData( val )  )
                        adapt_flag=newgrid.Adapt(values)
                    integral = newgrid.Current_Integral[0].val
                     
                    print newgrid.Errors[:10]
                    int1 = ((1.0/order**n) * gamma(n))**num_dim
                    self.assertTrue(np.allclose(integral, int1, rtol=1e-6, atol=1e-6))

    def TtestAdapt_mixed(self):
        from scipy.special import gamma
        """
            Test the sparse grid construction of mixed univariate rules via
            the integration of test functions.
            
            
            The integrand is a follows
            
            [uniform, gauss, Laguerre, Jacobi]
            
            [x1,x2,x3,x4]
            
            f(x1,x2,x3,x4) = 
        
        """
        

        rules=[1,5,8,9]
        num_dim=4
        LevelMax=6
        vt='single'

        uni_pdf = 2.0
        gauss_pdf = np.sqrt(np.pi)
        alpha=[i for i in range(1,5)]
        beta=alpha
        
        for n in range(1,5):
            for order in range(1,2):
                for a in alpha:
                    for b in beta:
                
                        Dist_args=[None, None, [n-1] , [b,a] ]
                        
                        beta_pdf = np.power(2,a+b+1)*gamma(a+1)*gamma(b+1) / gamma(a+b+2)
                        lag_pdf = gamma(n)
                        
                        Nominal = singleData( intfunc_mixed(np.zeros(num_dim,float), order)  )
                        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal, Distribution_Args=Dist_args,
                                      Adaptive=True, AdaptTol=1.0E-15,MaxVals=1000)
                        
                        self.assertTrue(np.allclose(sum(newgrid.Weights), uni_pdf * gauss_pdf * beta_pdf * lag_pdf, rtol=1e-6, atol=1e-6))
                        
                        
                        adapt_flag=True
                        while adapt_flag:
                            values=[]
                            for p in newgrid.ComputeIndex:
                                val=intfunc_mixed(newgrid.Points[:,p], order)
                                values.append( singleData( val )  )
                            adapt_flag=newgrid.Adapt(values)
                            integral = newgrid.Current_Integral[0].val
                        
                         
                        int1=(2.0/3.0) * np.sqrt(np.pi) * (1.0/order**n) * gamma(n) * (a-b)/float(a+b+2)
                        print integral/beta_pdf, int1, newgrid.Num_Points, newgrid.Errors[:10]
                        self.assertTrue(np.allclose(integral/beta_pdf, int1, rtol=1e-3, atol=1e-3))



    def TtestTensorProd_uni(self):
        """
            check the tensor product construction via the evaluation
            of some test integrals.
                    
        """
        
        LevelMax=2
        r_list=[1,2,3,4]
        exponent=2
        
        for num_dim in range(1,5):
            for r in r_list:
                #
                # Check integration works with list data
                #    
                vt='list'
                listlen=3
                rules=r*np.ones(num_dim,int)
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Listlen=listlen, TensorProd=True )
                
                
                Idx=range(newgrid.Num_Points)
                newgrid.Values=[listData(0.0,np.zeros(listlen,float)) for i in range(newgrid.Num_Points)]
                
                for i in Idx:
                    val=intfunc_uni(newgrid.Points[:,i], exponent)
                    newgrid.Values[i].val = val
                    newgrid.Values[i].list = val*np.ones(3,float)
        
                integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                
                exact=intfunc_uni_analytic(exponent, num_dim)
                self.assertTrue(np.allclose(integral[0].val, exact, rtol=1e-10, atol=1e-10))
                self.assertTrue(np.allclose(integral[0].list, exact*np.ones(listlen,float), rtol=1e-10, atol=1e-10))
                
                #
                # Check integration works with single data
                #
                vt='single'
                exponent=2
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, TensorProd=True)
                
                Idx=range(newgrid.Num_Points)
                newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
                for i in Idx:
                    newgrid.Values[i].val = intfunc_uni(newgrid.Points[:,i], exponent)
                    
                integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                exact=intfunc_uni_analytic(exponent, num_dim)
                self.assertTrue(np.allclose(integral[0].val, exact, rtol=1e-10, atol=1e-10))
        
    def TtestTensorProd_gauss(self):
        """
             Tests the ComputeGrid method which, given an index set, computes
             sparse grid weights and points.
             
             All of the Gaussian weight rules are tested:
             
             rule 5: Gauss Hermite
             rule 6: Generalized Gauss Hermite
             rule 10: Hermite Genz-Keister
             
             
             The integrand is
              
                 sum_{i=1}^N x_i^k prod_{j=1}^N dx_j
              
             with solution
              
                I  = 2^N * N/(k+1)    for k even
                I  = 0                for k odd
        """
  
        q_rules=[5,6,10]
        LevelMax=2
        vt='single'
  
        for q in q_rules:
            for exponent in range(1,4):
                for num_dim in range(1,4):
                    
                    if q==10: min(LevelMax,4)
  
                    rules=q*np.ones(num_dim,int)
                      
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, TensorProd=True)
                
                    Idx=range(newgrid.Num_Points)
                    newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
                    for i in Idx:
                        newgrid.Values[i].val = intfunc_gauss(newgrid.Points[:,i], exponent)
                        
                    integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                    self.assertTrue(np.allclose(integral[0].val, num_dim*intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))
                      
                    newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
                    for i in Idx:
                        newgrid.Values[i].val = intfunc_gauss2(newgrid.Points[:,i],exponent)
                     
                    integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
                    print integral[0].val, intfunc_gauss_analytic2(exponent, num_dim)
                    self.assertTrue(np.allclose(integral[0].val, intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))
   
    def TtestTensorProd_Laguerre(self):
        """
            Test the construction and application of Gauss-Laguerre quadratures
            via the evaluation of test integrals.
            
            The weight function and support for the Laguerre polynomials
            are
            
                supp. [0,inf]
                w(x) = exp(-x)
            
            The weight and support for the Generalize Lag polys are
            
                supp. [0,inf]
                w(x) = x^a exp(-x)    a > -1
            
            The following test integrals have been defined for Gauss-Laguerre integrals:
            
            1)  int exp(-px) = 1/p  [0,inf]
                -->  int exp(-x) * exp( x[1-p] ) = 1/p
                
            2) int exp(-p^2 x^2) = sqrt(pi) / (2p)
               -->   int exp(-x) * exp( x[1-p^2 x]) = sqrt(pi) / (2p)
            
             
            Integrals for the Generalized GL integrals:
             
            3) int (1/sqrt(x)) * exp (-px) = sqrt(pi/p)
        
        
            4) int x^{n-1} exp (-px) = (1/p^n) * G(n)
            
            G: gamma function
        """
        
        from scipy.special import gamma
        
        vt='single'
        
        LevelMax=5
        order_max = 5
        num_dim=2
        
        # The basic Laguerre quadratures
        q=7
        for num_dim in range(1,3):
            for order in range(1,order_max):
                Nominal = singleData( np.exp(0.0) )
                        
                rules=q*np.ones(num_dim,int)
                                           
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Nominal, TensorProd=True)
                         
                integral=0.0
                for p in range(newgrid.Num_Points):
                    temp = Laguerre_test(newgrid.Points[:,p], order)
                    integral += newgrid.Weights[p] * temp
                          
                int1=(1.0/order)**num_dim
                int2=(np.sqrt(np.pi)/(2.0*order))**num_dim
                self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))
                self.assertTrue(np.allclose(integral[1], int2, rtol=1e-4, atol=1e-3))   
                
        
            # The generalised Laguerre quadratures
            q=8
            Dist_args=[[-0.5] for dummy in range(num_dim)]  # alpha = -0.5
            for order in range(1,order_max):
                        
                rules=q*np.ones(num_dim,int)
                                           
                newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args, TensorProd=True)
                         
                integral=0.0
                for p in range(newgrid.Num_Points):
                    temp = Laguerre_test(newgrid.Points[:,p], order)
                    integral += newgrid.Weights[p] * temp
                 
                int1=np.sqrt(np.pi/order)**num_dim
                self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))

            for n in range(1,5):
                Dist_args=[[n-1] for dummy in range(num_dim)]  # alpha = -0.5
                for order in range(1,order_max):
                            
                    rules=q*np.ones(num_dim,int)
                                               
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args, TensorProd=True)
                    
                    integral=0.0
                    for p in range(newgrid.Num_Points):
                        temp = Laguerre_test(newgrid.Points[:,p], order)
                        integral += newgrid.Weights[p] * temp
                     
                    int1 = ((1.0/order**n) * gamma(n))**num_dim
                    self.assertTrue(np.allclose(integral[0], int1, rtol=1e-6, atol=1e-6))

    def TtestTensorProd_Jacobi(self):
        """
        Test the Gauss-Jacobi quadrature rules.
        
        The weight function:
        
            w(x,a,b) = (1-x)^a * (1+x)^b
        
        
        
        """
        from scipy.special import gamma

        vt='single'
        LevelMax=3
        
        alpha=[i for i in range(1,3)]
        beta=alpha

        q=9
        for num_dim in range(1,5):
            for a in alpha:
                for b in beta:
                    if a==b: continue
                    Dist_args=[[b,a] for i in range(num_dim)]  # alpha = -0.5
                            
                    rules=q*np.ones(num_dim,int)
                    
                    pdf = (np.power(2,a+b+1)*gamma(a+1)*gamma(b+1) / gamma(a+b+2))**num_dim
                                               
                    newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args, TensorProd=True)
                    self.assertTrue(np.allclose(sum(newgrid.Weights), pdf , rtol=1e-6, atol=1e-6))
                             
                    integral=0.0
                    for p in range(newgrid.Num_Points):
                        x=newgrid.Points[:,p]
                        temp = Jacobi_test(x)
                        integral += newgrid.Weights[p] * temp
                     
                    mean=integral[1]/pdf
                    var= integral[2]/pdf - mean**2
                    m2=(a-b)/float(a+b+2)
                    uv = 4*(a+1)*(b+1)/float(( (a+b+2)**2 * (3+a+b) )) + m2**2
                    v2= uv**num_dim - np.power(m2,2*num_dim)
                    #print "%5i %5i %5i %15.5e %15.5e %10.5f %10.5f %10.5f" %(a,b, num_dim, mean, m2**num_dim, mean/m2**num_dim, sum(newgrid.Weights), pdf)
                    self.assertTrue(np.allclose(m2**num_dim, mean , rtol=1e-6, atol=1e-6))
                    
                    self.assertTrue(np.allclose(v2, var , rtol=1e-6, atol=1e-6))
                    
    def TtestTensorProd_mixed(self):
        from scipy.special import gamma
        """
            Test the sparse grid construction of mixed univariate rules via
            the integration of test functions.
            
            
            The integrand is a follows
            
            [uniform, gauss, Laguerre, Jacobi]
            
            [x1,x2,x3,x4]
            
            f(x1,x2,x3,x4) = 
        
        """
        

        rules=[1,5,8,9]
        num_dim=4
        LevelMax=2
        vt='single'

        uni_pdf = 2.0
        gauss_pdf = np.sqrt(np.pi)
        alpha=[i for i in range(1,5)]
        beta=alpha
        
        for n in range(1,5):
            for order in range(1,2):
                for a in alpha:
                    for b in beta:
                
                        Dist_args=[None, None, [n-1] , [b,a] ]
                        
                        beta_pdf = np.power(2,a+b+1)*gamma(a+1)*gamma(b+1) / gamma(a+b+2)
                        lag_pdf = gamma(n)
                        
                        newgrid = sparse_grid(rules, num_dim, LevelMax, vt, Distribution_Args=Dist_args, TensorProd=True)
                        
                        self.assertTrue(np.allclose(sum(newgrid.Weights), uni_pdf * gauss_pdf * beta_pdf * lag_pdf, rtol=1e-6, atol=1e-6))
                        
                        
                        integral=0.0
                        for p in range(newgrid.Num_Points):
                            temp = intfunc_mixed(newgrid.Points[:,p], 1.0)
                            integral += newgrid.Weights[p] * temp
                         
                         
                        int1=(2.0/3.0) * np.sqrt(np.pi) * (1.0/order**n) * gamma(n) * (a-b)/float(a+b+2)
                        self.assertTrue(np.allclose(integral/beta_pdf, int1, rtol=1e-3, atol=1e-3))
                    
                    
    def TtestLargestErrorFowardNeighbours(self):
        """
             This method tests "LargestErrorFowardNeighbours" of the "sparse_grid"
             class from the LibSG module.
              
             "LargestErrorFowardNeighbours" finds the index with the largest error
             and computes its admissible forward neighbours.
              
        """
        num_dim=1
        Nominal = 1.0
        rules=3*np.ones(num_dim,int)
        vt='single'
  
        # check that the maximum level is not exceeded (one dimension only)
        for level in range(2,7):
            newgrid = sparse_grid(rules, num_dim, level, vt, Nominal,  Adaptive=True, AdaptTol=1.0E-10, MaxVals=20)
            LargestErrorFowardNeighbours_genidx(level, newgrid)
            newidx=newgrid.LargestErrorFowardNeighbours()
            self.assertEqual(len(newidx), 0)
          
        # check the full tensor grid is allowed and check levelMax is not exceeded (2D)
        num_dim=2
        rules=3*np.ones(num_dim,int)
        for level in range(2,7):
            newgrid = sparse_grid(rules, num_dim, level, vt, Nominal,  Adaptive=True, AdaptTol=1.0E-10, MaxVals=100)
            LargestErrorFowardNeighbours_genidx2D(level,newgrid)
            newidx=newgrid.LargestErrorFowardNeighbours()
              
            check= newgrid.Idx[:,newidx][:,0] == ([level,level])
            self.assertEqual(all(check),True)
            newgrid.Errors[newgrid.Active[:newgrid.N_Active][0]] = 10.0
            newidx=newgrid.LargestErrorFowardNeighbours()
            self.assertEqual(len(newidx), 0)

def intfunc_mixed(x,p):
    
    return x[0]**2 * x[1]**2 * np.exp(x[2]*(1.0-p)) * x[3]

def Jacobi_test(x):
    return np.asarray([reduce(mul,np.power(x,i)) for i in range(3)])

def Laguerre_test(x,p):
    return  np.asarray([ reduce(mul,np.exp(x*(1.0-p))) , reduce(mul,np.exp(x*(1.0-p**2*x))) ])

def LargestErrorFowardNeighbours_genidx2D(level, newgrid):
    counter = 1
    newgrid.N_Old=0; newgrid.N_Idx=0
    for i in range(level+1):
        for j in range(level+1):
            bnj= counter -1 
            bni= counter - (level+1)
            fnj= counter +1 
            fni= counter + (level+1)
            if i==0:bni=-1
            if j==0: bnj=-1
            
            if fnj>=(level+1)**2-1:fnj=-1
            if fni>=(level+1)**2-1:fni=-1
            
            newgrid.Idx[0,counter-1] = j; newgrid.Idx[1,counter-1] = i; newgrid.N_Idx += 1
            
            newgrid.N_Backward[0,counter-1] = bnj; newgrid.N_Backward[1,counter-1] = bni
            newgrid.Old[counter-1] = counter-1
            newgrid.N_Forward[0,counter-1] = fnj; newgrid.N_Forward[1,counter-1] = fni
            newgrid.N_Old +=1
            newgrid.Indicator[counter-1] = 0
            
            if counter >= (level+1)**2-1:
                newgrid.N_Old-=1
                newgrid.N_Forward[:,counter-1] = -1
                newgrid.Indicator[counter-1] = 1
                newgrid.Active[0] = counter-1
                newgrid.N_Active = 1 
                newgrid.heap_length=1
                newgrid.Errors[counter-1] = 5.0
                break
            counter+=1
            
def LargestErrorFowardNeighbours_genidx(level, newgrid):
    for i in range(2,level+1):
        newgrid.Idx[0,i] = i
        newgrid.N_Backward[:,i] = i
        newgrid.N_Forward[:,i-1] = i+1
        newgrid.N_Idx += 1
        
    for i in range(1,level):
        newgrid.Old[i] = i
        newgrid.N_Old += 1
        newgrid.Indicator[i] = 0
        
    newgrid.Active[newgrid.N_Active-1] = level
    newgrid.heap_length=1
    newgrid.Errors[level-1] = 5.0

def intfunc_gauss(x,exp):
    return sum(np.power(x,exp))   

def intfunc_gauss2(x,exp):
    return np.prod(np.power(x,exp))   

def intfunc_gauss_analytic2(exponent, num_dim):
    if np.mod(exponent,2) == 0:
        return np.power(np.pi, num_dim/2.0)
    else:
        return 0.0
    
def intfunc_uni_analytic(exponent, num_dim):
    if np.mod(exponent,2) == 0:
        return np.power(2,num_dim)*(num_dim/float(exponent+1))
    else: 
        return 0.0

def intfunc_uni(x,exp):
    return sum(np.power(x,exp))

