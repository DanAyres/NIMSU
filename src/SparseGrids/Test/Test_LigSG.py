import unittest
from SparseGrids.LibSG import sparse_grid, Compute_Grid, CalcSparseSet
from FORTRAN_Modules.sandiainterface import reduced_sg as Sandia
import numpy as np
import sys
from DataTypes.Results import singleData, listData, hdf5Data

class TestLibSG(unittest.TestCase):


    def testErrorMeasure(self):
        num_dim=1
        LevelMax=4
        rules=3*np.ones(num_dim,int)
        exponent=2
        nom=intfunc_uni(np.zeros(num_dim,float),exponent)
        
        new_int=[listData(1.0,np.zeros(3,float)),listData(3.0,np.zeros(3,float))]
        old_int=[listData(0.5,np.zeros(3,float)),listData(1.0,np.zeros(3,float))]
        
        vt='list'
        wmix = 0.0
        Nominal=listData(nom,listt=nom*np.ones(3,float))
        newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix, ValueType=vt)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 2.0)        

        wmix = 1.0
        Nominal=listData(nom,listt=nom*np.ones(3,float))
        newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix, ValueType=vt)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 1.0)


        vt='single'
        wmix = 0.0
        Nominal=singleData(nom)
        newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix, ValueType=vt)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 2.0)        

        wmix = 1.0
        Nominal=listData(nom)
        newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, wmix=wmix, ValueType=vt)
        Error=newgrid.ErrorMeasure(new_int, old_int)
        self.assertEqual(Error, 1.0)

    def testAdapt1_valtype(self):
        num_dim=1
        LevelMax=4
        rules=3*np.ones(num_dim,int)
        

        # Check integration works with list data
        vt='list'
        exponent=2
        nom=intfunc_uni(np.zeros(num_dim,float),exponent)
        Nominal=listData(nom,nom*np.ones(3,float))
        newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, ValueType=vt)
        
        Idx=range(newgrid.Num_Points)
        newgrid.Values=[listData(0.0,np.zeros(3,float)) for i in range(newgrid.Num_Points)]
        
        for i in Idx:
            val=intfunc_uni(newgrid.Points[:,i], exponent)
            newgrid.Values[i].val = val
            newgrid.Values[i].list = val*np.ones(3,float)

        integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
        
        exact=intfunc_uni_analytic(exponent, num_dim)
        self.assertTrue(np.allclose(integral[0].val, exact, rtol=1e-10, atol=1e-10))
        self.assertTrue(np.allclose(integral[0].list, exact*np.ones(3,float), rtol=1e-10, atol=1e-10))
    
        # Check integration works with single data
        vt='single'
        exponent=2
        Nominal=intfunc_uni(np.zeros(num_dim,float),exponent)
        newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal, Adaptive=False, AdaptTol=1.0E-10, MaxVals=20, ValueType=vt)
        
        Idx=range(newgrid.Num_Points)
        newgrid.Values=[singleData(0.0) for i in range(newgrid.Num_Points)]
        
        for i in Idx:
            newgrid.Values[i].val = intfunc_uni(newgrid.Points[:,i], exponent)
            
        integral=newgrid.CalcIntegral(newgrid.Weights, Idx)
        exact=intfunc_uni_analytic(exponent, num_dim)
        self.assertTrue(np.allclose(integral[0].val, exact, rtol=1e-10, atol=1e-10))
     
    def testComputeGrid_uni(self):
        """
             This tests the "Compute_Grid" grid method from the LibSG module.
              
             A series of integrals are performed on [-1,1] using a constant weight = 1.
              
             The integrand is
              
                 sum_{i=1}^N x_i^k prod_{j=1}^N dx_j
              
             with solution
              
                I  = 2^N * N/(k+1)    for k even
                I  = 0                for k odd
                 
              
        """
          
        q_rules=[1,2,3,4]
          
        LevelMax=4
        for q in q_rules:
            print 'rule = ', q
            for exponent in range(1,9):
                for num_dim in range(1,7):
                    rules=q*np.ones(num_dim,int)
                    growth=np.zeros(num_dim,int)
                    q_max=LevelMax*num_dim
                    sc=np.zeros(num_dim,int)
                    p=[]
                    tol=np.sqrt( sys.float_info.epsilon )
             
                    Idx=CalcSparseSet(0, LevelMax, num_dim)
                    Coeff= Sandia.calculate_coefficients(Idx, q_max)
                    points, weights =Compute_Grid(Idx, Coeff, q_max, rules, growth, LevelMax, sc, p, tol)
                    Num_Points=np.shape(points)[1]
                     
                    integral=0.0
                    for i in range(Num_Points):
                        integral += intfunc_uni(points[:,i],exponent) * weights[i]
                         
                    self.assertTrue(np.allclose(integral, intfunc_uni_analytic(exponent, num_dim), rtol=1e-10, atol=1e-10))
  
    def testComputeGrid_gauss(self):
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
  
        for q in q_rules:
            print 'rule = ', q
            for exponent in range(1,2):
                for num_dim in range(1,7):
                    
                    if q==10: LevelMax=4
  
                    rules=q*np.ones(num_dim,int)
                    growth=np.zeros(num_dim,int)
                    q_max=LevelMax*num_dim
                    sc=np.zeros(num_dim,int)
                    p=[]
                    tol=np.sqrt( sys.float_info.epsilon )
                      
                    #exponent=0
                      
                    Idx=CalcSparseSet(0, LevelMax, num_dim)
                    Coeff= Sandia.calculate_coefficients(Idx, q_max)
                    points, weights =Compute_Grid(Idx, Coeff, q_max, rules, growth, LevelMax, sc, p, tol)
                    Num_Points=np.shape(points)[1]
                     
                    integral=0.0
                    for i in range(Num_Points):
                        integral += intfunc_gauss(points[:,i],exponent) * weights[i]
                     
                    self.assertTrue(np.allclose(integral, num_dim*intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))
                      
                    integral=0.0
                    for i in range(Num_Points):
                        integral += intfunc_gauss2(points[:,i],exponent) * weights[i]
                     
                    self.assertTrue(np.allclose(integral, intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))


    def testComputeGrid_gauss_adapt(self):
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

        # Test for a single value
        vt='single'
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,7):
                    
                    Nominal = singleData(intfunc_gauss(np.zeros(num_dim,float),exponent))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=100, ValueType=vt)
                                
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
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,7):
                    
                    val=intfunc_gauss(np.zeros(num_dim,float),exponent)
                    Nominal = listData(val,val*np.ones(2,float))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=100, ValueType=vt)
                                
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            val=intfunc_gauss(newgrid.Points[:,p],exponent)
                            values.append( listData(val,val*np.ones(2,float)) )
                                
                        adapt_flag=newgrid.Adapt(values)
                    
                    integral = newgrid.Current_Integral[0].val 
                    self.assertTrue(np.allclose(integral, num_dim*intfunc_gauss_analytic2(exponent, num_dim), rtol=1e-6, atol=1e-6))
                    


    def testComputeGrid_uni_adapt(self):

  
        q_rules=[1,2,3,4]
        LevelMax=6

        # Test for a single value
        vt='single'
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,7):
                    
                    Nominal = singleData(intfunc_uni(np.zeros(num_dim,float),exponent))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=1000, ValueType=vt)
                                
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
        for q in q_rules:
            for exponent in range(1,3):
                for num_dim in range(1,7):
                    
                    val=intfunc_uni(np.zeros(num_dim,float),exponent)
                    Nominal = listData(val,val*np.ones(2,float))
              
                    rules=q*np.ones(num_dim,int)
                                  
                    newgrid = sparse_grid(rules, num_dim, LevelMax, Nominal,  Adaptive=True, AdaptTol=1.0E-8, MaxVals=1000, ValueType=vt)
                                
                    adapt_flag=True
                    while adapt_flag:
                        values=[]
                        for p in newgrid.ComputeIndex:
                            val=intfunc_uni(newgrid.Points[:,p],exponent)
                            values.append( listData(val,val*np.ones(2,float)) )
                                
                        adapt_flag=newgrid.Adapt(values)
                    
                    integral = newgrid.Current_Integral[0].val 
                    self.assertTrue(np.allclose(integral, intfunc_uni_analytic(exponent, num_dim), rtol=1e-6, atol=1e-6))
                    

  
    def testLargestErrorFowardNeighbours(self):
        """
             This method tests "LargestErrorFowardNeighbours" of the "sparse_grid"
             class from the LibSG module.
              
             "LargestErrorFowardNeighbours" finds the index with the largest error
             and computes its admissible forward neighbours.
              
        """
        num_dim=1
        Nominal = 1.0
        rules=3*np.ones(num_dim,int)
  
        # check that the maximum level is not exceeded (one dimension only)
        for level in range(2,7):
            newgrid = sparse_grid(rules, num_dim, level, Nominal,  Adaptive=True, AdaptTol=1.0E-10, MaxVals=20)
            LargestErrorFowardNeighbours_genidx(level, newgrid)
            newidx=newgrid.LargestErrorFowardNeighbours()
            self.assertEqual(len(newidx), 0)
          
        # check the full tensor grid is allowed and check levelMax is not exceeded (2D)
        num_dim=2
        rules=3*np.ones(num_dim,int)
        for level in range(2,7):
            newgrid = sparse_grid(rules, num_dim, level, Nominal,  Adaptive=True, AdaptTol=1.0E-10, MaxVals=100)
            LargestErrorFowardNeighbours_genidx2D(level,newgrid)
            newidx=newgrid.LargestErrorFowardNeighbours()
              
            check= newgrid.Idx[:,newidx][:,0] == ([level,level])
            self.assertEqual(all(check),True)
            newgrid.Errors[newgrid.Active[:newgrid.N_Active][0]] = 10.0
            newidx=newgrid.LargestErrorFowardNeighbours()
            self.assertEqual(len(newidx), 0)



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