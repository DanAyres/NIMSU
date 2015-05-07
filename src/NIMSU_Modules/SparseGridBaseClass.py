
import numpy as np
from FORTRAN_Modules.maxheap import maxheap
from FORTRAN_Modules.sandiainterface import reduced_sg as Sandia
from FORTRAN_Modules.ggmethods import gerstnergriebel as ggmeth

from NIMSU_Modules.DataType_Results import singleData, listData

import sys
import itertools
from operator import mul


class sparse_grid:

    def __init__(self, rules, num_dim, LevelMax, ValueType, Nominal=0.0,  Adaptive=False, AdaptTol=1.0E-3, MaxVals=20, wmix=0.0, TensorProd=False):
        """
        
            Initialise the sparse grid class. 
            
            If the NON-adaptive option is selected; a standard sparse grid will 
            be constructed with maximum level "LevelMax".
            
            When the adaptive option is selected all of the necessary data structures
            are initialised and the next set of points to be computed are stored.
            
            The adaptive procedure is as follows:
                1) create new sparse_grid instance
                   newgrid = sparse_grid(x,y,z..)
                
                2)  flag=True
                    while flag:    
                        values=[]
                        for p in newgrid.ComputeIndex:
                            values.append( ENGINE(newgrid.Points[:,p]) )
        
                        flag = newgrid.Adapt(np.asarray(values))

            The growth rate of the quadrature rules are hard wired but may be
            changed to any of the following (rule dependent):
            
                growth : 
                    0, "DF", default growth associated with this quadrature rule
                    1, "SL", slow linear, L+1;
                    2  "SO", slow linear odd, O=1+2((L+1)/2)
                    3, "ML", moderate linear, 2L+1;
                    4, "SE", slow exponential;
                    5, "ME", moderate exponential;
                    6, "FE", full exponential.
            
            Args:
        
                rules: list of length num_dim:
                
                     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
                     2, "F2",  Fejer Type 2, Open Fully Nested.
                     3, "GP",  Gauss Patterson, Open Fully Nested.
                     4, "GL",  Gauss Legendre, Open Weakly Nested.
                     5, "GH",  Gauss Hermite, Open Weakly Nested.
                     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
                     7, "LG",  Gauss Laguerre, Open Non Nested.
                     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
                     9, "GJ",  Gauss Jacobi, Open Non Nested.
                    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
                    11, "UO",  User supplied Open, presumably Non Nested.
                    12, "UC",  User supplied Closed, presumably Non Nested.
    
                num_dim: number of dimensions
                
                LevelMax: the maximum quadrature level:
                
                Nominal:  the value of the integrand at (0,0,0,...)
                
                Adaptive: Toggle the adaptive grid
                
                AdaptTol: Error tolerance of the integration
                
                MaxVals:  Largest number of indexes allowed
                
                wmix:    (1-wmix) * var + wmix * mean
                
                ValueType:    determine how the results are stored.
                              ='single': a single result is stored in the Values array.
                              ='list':  each result consists of a value and a list of values in a list
                              e.g.    [val, [a1,b2,...] ].
                              ='hdf5': each result consists of a value and an hdf5 filename.
                              e.g.    [val,filemane]
                    
        """
        
        np.set_printoptions(edgeitems=3,infstr='inf', \
                            linewidth=150, nanstr='nan', precision=8,\
                            suppress=False, threshold=1000, formatter=None)        
                
        self.rules=rules
        self.num_dim=num_dim
        self.growth=np.zeros(self.num_dim,int) 
        self.LevelMax=LevelMax
        self.AdaptTol=AdaptTol
        self.Nominal = Nominal
        self.wmix=wmix
        self.q_max = self.LevelMax * self.num_dim
        
        self.ValueType = ValueType.lower()
        
        # calculate the length of the results list
        if self.ValueType == 'list':
            self.Results_ListLen = len(self.Nominal.list)
        
        # Check that LevelMax is allowed for the quadrature rule
        if any([r==10 for r in rules]):
            if self.LevelMax > 4:
                self.LevelMax=4
        
        # Some required parameters
        self.sc=np.zeros(self.num_dim,int)
        self.p=[]
        self.tol=np.sqrt( sys.float_info.epsilon )
        
        
        self.Adaptive=Adaptive
        
        
        if TensorProd:
            self.Points, self.Weights, self.Num_Points = self.TensorProduct()
            return
        
        if self.Adaptive==True:
        
            # Maximum memory allowed for the indexes
            self.MaxVals=MaxVals
            
            
            # Maximum memory allowed for points, weights and values
            self.MaxPoints = 5*self.MaxVals
        
        
            # Allocate memory for the Gerstner-Griebel scheme
            self.Active=np.zeros(self.MaxVals,int)
            self.Old=np.zeros(self.MaxVals,int)
            self.Indicator=np.ones(self.MaxVals,int)
            self.N_Forward=-1*np.ones([self.num_dim,self.MaxVals],int)
            self.N_Backward=-1*np.ones([self.num_dim,self.MaxVals],int)
            self.Idx=np.zeros([self.num_dim,self.MaxVals],int)
            self.Errors=-1.0*np.ones(self.MaxVals,float)
            
            
            # Initiialise the sparse grid
            self.N_Idx = self.num_dim+1
            self.N_Old = 1
            self.N_Active = self.num_dim
            self.Init_Adaptive_Grid()


            # Allocate memory for points and values
            if self.ValueType == 'single':
                self.Values=[singleData(0.0) for dummy in range(self.MaxPoints)]
            elif self.ValueType == 'list':
                self.Values=[listData([0.0,np.zeros(self.Results_ListLen,float)]) for dummy in range(self.MaxPoints)]
                
            #self.Values=np.zeros(self.MaxPoints,float)
            #self.Values=[[] for dummy in range(self.MaxPoints)]
            self.Points=np.zeros([self.num_dim,self.MaxPoints],float)
            
            
            # Initialise points, weights, and values
            self.Values[0] = self.Nominal
            self.Points[:,:1], self.Weights = self.Compute_Grid(self.Idx[:,0])
            self.Index = np.zeros(1,int)
            self.Num_Points = 1
            
            
            # Calculate the "Active Points"
            Active_Points,dummy  = self.Compute_Grid(self.Idx[:,self.Active[:self.N_Active]])
            Index, Unique = Sandia.concatenate(self.Points[:,:self.Num_Points],Active_Points)
            self.ComputeIndex = self.CalcComputeIndex(Index, Unique, Active_Points)
            
            # The heap for the adaptive procedure is uninitialized 
            self.heap_init = False
            self.heap_length = 0
            
            
        else:
        
            self.Idx = CalcSparseSet(0, self.LevelMax, self.num_dim)
            
            self.Points, self.Weights = self.Compute_Grid(self.Idx)
            self.Num_Points = np.shape(self.Points)[1]
            self.Index=range(self.Num_Points)

    def TensorProduct(self):
        
        Idx = CalcSparseSet(0, self.LevelMax, 1)
        seed = 123456789
        growth=np.zeros(1,int)
        sc=np.zeros(1,int)
        pts=[]
        wts=[]
        Num_Points=1
        for r in self.rules:

            Coeff= Sandia.calculate_coefficients(Idx, self.q_max)
            new_np = Sandia.max_next_points(Idx, Coeff, [r], growth)
            points = Sandia.weights_and_points(new_np, self.LevelMax, Idx, Coeff, growth, [r], sc, self.p)
        
            N_Unique, sparse_index = Sandia.unique_points(seed, self.tol, points)
            Points, Weights=Sandia.reduce_points_and_weights(N_Unique, points, Idx, sparse_index, Coeff, growth, [r], sc, self.p)
                        
            pts.append(Points[0,:])
            wts.append(Weights)
            Num_Points *=np.shape(Points)[1]
            
        tmp=itertools.product(*pts)
        Points=np.zeros([self.num_dim,Num_Points],float)
        for i,j in enumerate(tmp):
            Points[:,i] = np.asarray(j)
            
        tmp=itertools.product(*wts)
        Weights=np.zeros([Num_Points],float)
        for i,j in enumerate(tmp):
            Weights[i] = reduce(mul,j)
            
        return Points, Weights, Num_Points

    def Adapt(self, values):
        
        
        # Compute the current integral
        self.Current_Integral = self.CalcIntegral(self.Weights, self.Index)
        
        # Add the new values to the old
        for pos, p in enumerate(self.ComputeIndex):
            self.Values[p] = values[pos]
            
        # Calculate the new error measures
        self.CalcNewErrorMeasures()


        # Check the error for convergence
        if sum(self.Errors[self.Active[:self.N_Active]]) < self.AdaptTol :
            return False
        
        
        # Remove the largest error and calculate the forward neighbours        
        new_idx=self.LargestErrorFowardNeighbours()    
        if len(new_idx) == 0: 
            self.ComputeIndex = []
            return True   
        
        # Calculate the points and the compute index for the new indices
        Active_Points,dummy  = self.Compute_Grid(self.Idx[:,new_idx])
        Index, Unique = Sandia.concatenate(self.Points[:,:self.Num_Points],Active_Points)
        self.ComputeIndex = self.CalcComputeIndex(Index, Unique, Active_Points)


        #Calculate the weights for the current integral
        Active_Points,self.Weights  = self.Compute_Grid(self.Idx[:,self.Old[:self.N_Old]])
        self.Index, dummy = Sandia.concatenate(self.Points[:,:self.Num_Points],Active_Points)
        self.Index = self.Index-1
    
        return True
    
    def LargestErrorFowardNeighbours(self):
        """
            Removes the largest error from the heap, removes the corresponding 
            index from the active set and adds it to the old set.
            
            All of the forward neighbours of this index are computed and added 
            to the active set.
           
           Args:
           
           Returns:
               new_idx:    the positions in self.Idx of the new active indices.
            
        """
        
        # Get the largest Error - remove from active and add to the old set
        self.Active, maxpos = maxheap.heap_delete(self.Active, self.Errors, self.heap_length)
        self.heap_length -=1
        self.N_Active -=1
        self.Old[self.N_Old] = maxpos
        self.N_Old += 1
        self.Indicator[maxpos] = 0
        
        # Calculate the forward neighbours and add them to the active set
        self.N_Backward,self.N_Forward,N_Idx, \
        self.Idx,self.Indicator = ggmeth.calculateneighbours(self.num_dim, self.N_Backward,\
                                                             self.N_Forward, maxpos+1, self.N_Idx,\
                                                             self.Idx, self.Indicator, self.LevelMax)
        new_idx=np.arange(self.N_Idx,N_Idx)
        self.N_Idx=N_Idx
        for i in new_idx:
            self.Active[self.N_Active] = i
            self.N_Active += 1        
    
        return new_idx
    
    def CalcNewErrorMeasures(self):
        """
            For all of the indexes in the self.Active set and
            error measure is computes (if it hasn't already).
            
            This new index is added to the heap.
            
            Args:
            
            Returns:
            
            Raises:
            
        """
        for p in self.Active[:self.N_Active]:
            if self.Errors[p] < 0.0:
                #print self.CalcErrorMeasure(p), p
                self.Errors[p] = self.CalcErrorMeasure(p)
                # Add new values to the heap
                self.Active[:self.heap_length+1],dummy= maxheap.heap_insert(self.Errors[:self.N_Idx], 
                                                                            p, self.Active[:self.heap_length+1],
                                                                            self.heap_length)
                self.heap_length +=1
                
        if self.heap_length != self.N_Active:
            raise ValueError
    
    def CalcErrorMeasure(self,Idxpos):
        """
           Given the position of a new index, a new integral is computed which
           includes this latest index. This is then used to compute an error 
           measure along with the current integral.
           
           Args:
               Idxpos: the position in self.Idx of the new index
               
           Returns:
               EM:    the error measure
            
        """
        
        pos=np.zeros([self.N_Old+1],int)
        pos[:self.N_Old] = self.Old[:self.N_Old]
        pos[-1] = Idxpos
        
        Active_Points,Weights  = self.Compute_Grid(self.Idx[:,pos])
            
        Index, dummy = Sandia.concatenate(self.Points[:,:self.Num_Points],Active_Points)
        Index = Index -1
        new_integral = self.CalcIntegral(Weights, Index)
        
        return self.ErrorMeasure(new_integral, self.Current_Integral) 

    def ErrorMeasure(self, new_int, old_int):
        """
            Error measure. Calculates an error measure
            based upon the old and new integrals
            
            Args:
                new_int:    the latest integral   (np array)
                old_int:    the previous integral (np array)
            
            Returns:
                EM:    The error measure (float)
        """
        
        
        mean=np.fabs(new_int[0].val - old_int[0].val)
        var=np.fabs(new_int[1].val - old_int[1].val)
        if old_int[1].val == 0.0:
            return (1.0-self.wmix) * var + self.wmix * mean
        else:
            return (1.0-self.wmix) * var / np.fabs(old_int[1].val) + self.wmix * mean / np.fabs(old_int[0].val)

    def CalcIntegral(self, Weights, Index):
        """
            Calculate the integral and square integral using quadrature.
            
            Args:
                Weights:    The quadrature weights
                Index:      The positions in the static "self.Points" array corresponding
                            to the weights. i.e. weights[j] -> self.Points[ Index[j] ]
            
            Returns:    
                list:       [integral, square integral]
        """
        
        
        if self.ValueType=='single':
            valm=singleData(0.0)
            valv=singleData(0.0)
        elif self.ValueType == 'list':
            valm=listData(0.0,np.zeros(self.Results_ListLen,float))
            valv=listData(0.0,np.zeros(self.Results_ListLen,float))
        
        for pos,p in enumerate(Index):
            #print 'int', valm.val, self.Values[p].val
            valm += Weights[pos] * self.Values[p]
            valv += Weights[pos] * self.Values[p] * self.Values[p]
        
        return [valm,valv]
    
    def CalcComputeIndex(self, Index, Unique, Points):
        """
            Calculate the index of the points in the static "self.Points"
            array which must be calculated before the next adaptive cycle.
            
            The self.Points array is updated and the index of the new point
            is added to the compute index array.
            
            Args:
                Index    : Potitions of the new points in the self.Points array.
                           It assumes the new, unique, points have already been added.
                Unique   : Boolean array corresponding to the new points array
                           which indicated if it already exists in self.Points.
                Points   : the new points array.
            
            Returns:
                Compute_Index    : Indexes of self.Points to compute.
        
        """
        
        ComputeIndex=[]
        for i in range( len(Index) ):
            if Unique[i]: 
                ComputeIndex.append(Index[i]-1)
                self.Points[:,Index[i]-1] = Points[:,i]
                self.Num_Points += 1
                
        return np.asarray(ComputeIndex)
    
    def Init_Adaptive_Grid(self,):
        """
            Initialise the adaptive sparse grid and all of its data types.
            
        
        """
        self.Indicator[0] = 0
        self.Old[0] = 0
    
        for i in range(self.num_dim):
            self.Active[i] = i+1
            self.Indicator[i+1] = 1
            self.N_Forward[i,0] = i+2
            self.N_Backward[i,i+1] = 1
            self.Idx[i,i+1] = 1

    def Compute_Grid(self,Idx):
        """
            Calculate a sparse grid quadrature scheme from a set of index
            vectors.
            
            Args:
                Idx:    Array of index vectors. e.g. Idx= [ [0,0], [0,1], [1,0] ] 
                
            Returns:
                Points:    Quadrature points
                Weight:    Quadrature weights
        """

        seed = 123456789
        Coeff= Sandia.calculate_coefficients(Idx, self.q_max)
        new_np = Sandia.max_next_points(Idx, Coeff, self.rules, self.growth)
        points = Sandia.weights_and_points(new_np, self.LevelMax, Idx, Coeff, self.growth, self.rules, self.sc, self.p)
        
        N_Unique, sparse_index = Sandia.unique_points(seed, self.tol, points)
        return Sandia.reduce_points_and_weights(N_Unique, points, Idx, sparse_index, Coeff, self.growth, self.rules, self.sc, self.p)


def Compute_Grid(Idx, Coeff, q_max, rules, growth, LevelMax, sc, p, tol, ):
    """
        Calculate a sparse grid quadrature scheme from a set of index
        vectors.
        
        Args:
            Idx:    Array of index vectors. e.g. Idx= [ [0,0], [0,1], [1,0] ] 
                
        Returns:
            Points:    Quadrature points
            Weight:    Quadrature weights
    """

    seed = 123456789
    #Coeff= Sandia.calculate_coefficients(Idx, q_max)
    new_np = Sandia.max_next_points(Idx, Coeff, rules, growth)
    points = Sandia.weights_and_points(new_np, LevelMax, Idx, Coeff, growth, rules, sc, p)
    N_Unique, sparse_index = Sandia.unique_points(seed, tol, points)
    return Sandia.reduce_points_and_weights(N_Unique, points, Idx, sparse_index, Coeff, growth, rules, sc, p)

def CalcSparseSet(level_min, level_max, dim_num):
    sparse_set=[]
    s=True
    for level in range(level_min,level_max+1):
        sparse_set.extend(level_sums(level,dim_num,s) )
    return np.asarray(map(list, zip(*sparse_set)))
    
def level_sums(n,r,s):
    """
        Return a vector of length r such that
        sum_i vec_i = n
    """
#    from numpy import zeros
    from copy import copy
    OutSet=[]
#    a=zeros(r,int)
    a=[0 for i in range(r)]
    if(s):
        t=n
        h=-1
        a[0]=n
        s=False
        OutSet.append(copy(a))
    else:
        return OutSet        
    while (a[r-1] != n):
        if(1<t):
            h=-1
        h+=1
        t=a[h]
        a[h]=0
        a[0]=t-1
        a[h+1]+=1
        OutSet.append(copy(a))
    else:
        return OutSet


