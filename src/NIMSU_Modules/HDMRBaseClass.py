'''
Created on 31 Mar 2015

@author: daiel
'''

from NIMSU_Modules.SparseGridBaseClass import sparse_grid
import numpy as np
from operator import mul, itemgetter
from itertools import combinations
from NIMSU_Modules.DataType_Results import singleData, listData
from FORTRAN_Modules.sandiainterface import reduced_sg as Sandia

class HDMR_Base():
    
    def __init__(self, control_info, engine, sampler):
        """
            Initialise the HDMR base class. The HDMR expansion is stored
            in a single data structure called data_struct. This data structure
            is a list of lists. The index of data_struct corresponds to the dimension
            of the HDMR terms. Each list has the following form:
            
            [set, error, mean, variance, sparse_grid, flag]  
            
            where
            
            @param set:     is the dimensions associated with the HDMR term. e.g. [1] or [2,5]
            @param error:   is the error measure for that term.  (float)
            @param mean:    "single" - float, "list" - [float,[float_1,float_2,..,float_{self.Results_ListLen}]]
                     "hdf5" - [float,hdf5_filename]
            var:    same as mean
            sparse_grid:    Sparse grid integration object
            flag:    identifies if the sparse grid require adapting
            
            
            Args:
                control_info:    run control information
                engine:          object used to run the physics code
                sampler:         object used to sample the input parameters
        """
        
        self.quadAdapt=control_info.hdmr_quad_adapt
        self.quadTol=control_info.hdmr_quad_tol
        self.quadLevels=control_info.hdmr_sg_level
        self.ValueType=control_info.ValueType
        self.R0=engine.Nominal
        self.num_dim=sampler.Total_stoch_dim
        self.Results_ListLen = control_info.Results_listlen
        self.tolerance=control_info.hdmr_metric
        
        self.data_struct=[]
        
        self.engine=engine
        
        
        self.Total_Realisations=0
        
        # Initialise the data structure with all of the 1D terms
        ds_temp=[]
        for i in range(self.num_dim):
            ds_temp.append([ [i], 0.0, 0.0, 0.0, \
                            None,\
                            True if self.quadAdapt else False ])
        self.data_struct.append(ds_temp)
        
        # Init the mean and variance
        if self.ValueType == 'single':
            self.totalMean = self.R0
            self.totalVar = singleData(0.0)
        elif self.ValueType == 'list':
            self.totalMean = self.R0
            self.totalVar = listData(0.0,np.zeros(self.Results_ListLen,float))
        self.localVar=0.0
            
        self.dim_split = 2

        # no return

    def compute(self, sampler):
        
        maximum = min(len(self.quadLevels)+1, self.num_dim+1)
        for dimension in range(1,maximum):
            self.localVar=0.0
            
            num_sets=nCr(self.num_dim, dimension)
    
            # We want to calculate more than one set at a time. This is a bottle-neck.        
            split=num_sets if dimension==1 else self.dim_split
            for i in range(0,num_sets,split):
                start=i
                end=min(i+split,num_sets)
            
                # Integrate
                adapt_Int_Complete=False
                while not adapt_Int_Complete: 
                    
                    # Generate the next set of samples
                    Samples = self.createSamples(dimension, sampler, self.data_struct[dimension-1][start:end])
                    
                    # Compute the results        
                    Results=self.getResults(Samples)
                    
                    # Compute the new integration points (adaptive procedure)
                    self.updateQuadIntegrals(Results, dimension, self.data_struct[dimension-1][start:end])
                    
                    # Check that the adaptive quadratures have finished
                    adapt_Int_Complete= all([data[5]==False for data in self.data_struct[dimension-1][start:end] ])
                    
                    

                # Compute HDMR functions, total statistics and local variance
                self.HDMR_functions(sampler.PDF, self.data_struct[dimension-1][start:end])
            
                self.CalcErrorMeasures(self.data_struct[dimension-1][start:end])
            
                # Test for "superposition" convergence
                if self.innerConvergenceTest(self.data_struct[dimension-1][start:end]):
                    del self.data_struct[dimension-1][end:]
                    break

            # Test for "dimension" convergence
            if self.outerConvergenceTest():
                break

            #for p,data in enumerate(self.data_struct[dimension-1]):
            #    print p+1,data
            
            # sort the list acording to error measure - largest first
            #self.data_struct[dimension-1] = sorted(self.data_struct[dimension-1],key=itemgetter(1), reverse=True)
            
            self.calcNextSets(dimension, sampler)
            
            #print 
        

    def outerConvergenceTest(self):
        
        if np.abs(self.localVar)/self.totalVar.val < self.tolerance:
            return True
        else:
            return False
        
    def innerConvergenceTest(self, HDMR_sets):
        """
            Test for convergence of the HDMR terms for a given dimension.
            The HDMR component functions for dimension >=2 are calculated 
            in batches of self.dim_split. The errors, which are defined as
            the relative variance, are averaged over the batch size. This 
            value is compared to a pre-defined error tolerance.
            
        """
        error=0.0
        for data in HDMR_sets:
            error += data[1]
        error /= min(self.dim_split,len(HDMR_sets))
             
        if error < self.tolerance:
            return True
        else: 
            return False
        
    def CalcErrorMeasures(self, HDMR_sets):
        """
            Calculate the error measures associated with each HDMR function.
            The error measure is defined as the relative variance associated with
            each function. i.e. the error, gk, for function k is given as
            
                gk = var_k / varTotal
                
            Args:
                HDMR_sets
        
        """
        
        # Compute the error measure - relative variance contribution
        for data in HDMR_sets:
            data[1] = np.abs(data[3].val)/np.abs(self.totalVar.val)
        
    def calcNextSets(self, dimension, sampler):
        """
            Calculate the next set of HDMR sets. One dimensional sets are combined with
            n dimensional sets to create n+1 dimensional sets. Fro example
            
            let n=2. We have the following one dimensional sets:
            
            [1],[2],[3]
            
            and two dimsnional sets:
            
            [1,2], [1,3], [2,3]
            
            The only admissible set is: [1,2,3]
            
            
            Args:
                dimension:    the current dimension
                sampler:      sampling object
            
            Returns:
                
        
        """
        
        
        ds_temp=[]
        
        for c_set in self.data_struct[dimension-1]:
            for b_set in self.data_struct[0]:
                
                # continue if the current set contains the base set
                if b_set[0][0] in c_set[0]:
                    continue
                else:
                    tmp=[]
                    for i in c_set[0]:
                        tmp.append(i)
                    for i in b_set[0]:
                        tmp.append(i)
                    tmp=sorted(tmp)
                
                # Check if this set has already been counted. Add if not
                if tmp in [row[0] for row in ds_temp]:
                    continue
                else:
                    ds_temp.append([tmp,b_set[1]*c_set[1], 0.0, 0.0, \
                                    None,\
                                    True if self.quadAdapt else False ])
        
        
        try:
            self.data_struct[dimension] = ds_temp
        except IndexError:
            self.data_struct.append(ds_temp)
    
        # Sort the new set according to the predicted error measures.
        self.data_struct[dimension] = sorted(self.data_struct[dimension],key=itemgetter(1), reverse=True)
        
        # no return
    
    def HDMR_functions(self, PDF, HDMR_sets):
        
        for data in HDMR_sets:
                
            pdf=reduce(mul,[PDF[p] for p in data[0]])
            grid=data[4]
                
            if self.ValueType=='single': 
                mean=singleData(0.0)
                var=singleData(0.0)
            elif self.ValueType=='list': 
                mean=listData(0.0,np.zeros(self.Results_ListLen,float))
                var=listData(0.0,np.zeros(self.Results_ListLen,float))
                
            # Sum the quadrature sets
            for w, val in zip(grid.Weights,grid.Values):
                mean+=w*val; var+=w*val*val
            
            # Compute the mean and variance
            mean*=pdf
            var=pdf*var-mean**2
                
            # sum all contributions from the subsets
            m,v=self.subsetContributions(data[0])
                
            # update the mean and variance
            data[2] = mean - m - self.R0
            data[3] = var - v
                
            # update the totals
            self.totalMean += data[2]
            self.totalVar += data[3]
            self.localVar += data[3].val
                
    def subsetContributions(self, currentSet):
        """
            Compute the subsets of a given set. For each HDMR function
            corresponding to these subsets, sum their means and variances.
            
            For example: consider the current set [1,2,5]. The subsets 
            are    [1],[2],[5],[1,2],[1,5] and [2,5].
            
            If all of these subsets exist, their means and variances are summed.
            
            Args:
                currentSet    : the current HDMR set
            Returns:
                mean:   sum of all subset means
                var:    sum of all subset variances
        """
        
        if self.ValueType=='single': 
            mean=singleData(0.0)
            var=singleData(0.0)
        elif self.ValueType=='list': 
            mean=listData(0.0,np.zeros(self.Results_ListLen,float))
            var=listData(0.0,np.zeros(self.Results_ListLen,float))
        
        subsets=sum(map(lambda r: list(combinations(currentSet, r)), range(1, len(currentSet))), [])
    
        if not subsets:
            return mean, var
        else:
                
            # Temporarily store the sets of the HDMR terms
            indx=[]
            for i in range(len(currentSet)-1):
                indx.append([data[0] for data in self.data_struct[i] ])
            
            # Loop through the subsets and add their mean and var (if they exist)
            for s in subsets:
                l=len(s)
                try:
                    pos=indx[l-1].index(list(s))
                except ValueError:
                    continue
                mean+=self.data_struct[l-1][pos][2]
                var+=self.data_struct[l-1][pos][3]
    
            return mean, var
        
    def updateQuadIntegrals(self, Results, dimension, HDMR_sets):
        """
            Add the recently computed values to the quadrature
            object. These new values are used to compute the next set of points.
            If the integral has sufficiently converged
        
        """
        
        pos=0
        for data in HDMR_sets:
            grid=data[4]
            if self.quadAdapt:
                if data[5] == False: continue
                values=[]
                for dummy in grid.ComputeIndex:
                    key='s'+str(pos)
                    pos+=1
                    values.append(Results[key]) 
                data[5]=grid.Adapt(values)
            else:
                grid.Values=[]
                for dummy in range(grid.Num_Points):
                    key='s'+str(pos)
                    pos+=1
                    grid.Values.append(Results[key])
                    
    def getResults(self, Samples):
        
        Results={ key:Samples[key] for key in Samples if not isinstance(Samples[key],dict) }
        Samples_n = { key:Samples[key] for key in Samples if isinstance(Samples[key],dict) }
        
        Results.update(self.engine.interface(Samples_n))
        #Results.update(EngineTemp(Samples_n, self.ValueType,self.Results_ListLen))
        
        return Results         
    
    def createSamples(self, dimension, sampler, HDMR_sets):
        """
            Create a dictionary of samples. Each sample is itself a dictionary of
            all the variable parameters. The key for each sample is a sequential
            integer value prefixed by an "s". e.g. s0, s1, s2,...,sn
            
            Args:
                dimension:  number of dimensions in the sample space.
                sampler:    sampling object created from the <SampleUncertainties> class
                            in the <samples> module.
            
            Returns:
                Samples: dictionaty of sampled values.
        
        """
        
        Samples={}
        rvs=np.zeros(self.num_dim,float)
        pos=0
        for data in HDMR_sets:
            
            # Init the sparse grids on-the-fly
            if data[4] == None:
                rules = [sampler.rules[i] for i in data[0]]
                data[4] = sparse_grid(rules, \
                                      dimension, \
                                      self.quadLevels[dimension-1], \
                                      self.R0, \
                                      Adaptive=self.quadAdapt, \
                                      AdaptTol=self.quadTol, \
                                      MaxVals=40, \
                                      ValueType=self.ValueType)
            
            grid=data[4]
            
            if self.quadAdapt:
                if data[5] == False: continue
                Point_Locs=grid.ComputeIndex
            else:
                Point_Locs = range(grid.Num_Points)
                
            for p in Point_Locs:
                key='s'+str(pos)
                rvs[data[0]] = grid.Points[:,p]

                # check for all zeros (Nominal)
                if np.allclose(rvs, np.zeros(self.num_dim,float), rtol=1e-10, atol=1e-10):
                    Samples[key] = self.R0
                    
                # Points may exist in one of the lower dimensional sub-spaces
                elif dimension > 1:
                    Samples[key]= self.findPreviousPoints(data[0], grid.Points[:,p], sampler, rvs)

                # One dimensional terms always have to be sampled
                else:
                    self.Total_Realisations+=1
                    Samples[key]=sampler.sample(rvs)

                rvs[data[0]] = 0.0
                pos +=1

        return Samples

    def findPreviousPoints(self, dims, points, sampler, rvs):
        
        if any(points==0.0): # These may exist somewhere in a lower subspace
            
            
            reduced_dims = [d for d,p in zip(dims,points) if p != 0.0  ]
            
            dimension = len(reduced_dims)
            
            dim_sets = [data[0] for data in self.data_struct[dimension-1] ]
            
            try:
                
                # Find the right HDMR set
                pos=dim_sets.index(reduced_dims)
                grid = self.data_struct[dimension-1][pos][4]
            
                Index, Unique = Sandia.concatenate(grid.Points,points[reduced_dims])
                if Unique[0] == 0:
                    return grid.Values[Index[0]-1]
                else:
                    # point doesn't exist
                    self.Total_Realisations+=1
                    return sampler.sample(rvs)


            except IndexError:
                # subset doesn't exist - generate new sample
                self.Total_Realisations+=1
                return sampler.sample(rvs)
        else:
            # No points lie of the axis - generate a new sample
            self.Total_Realisations+=1
            return sampler.sample(rvs)

def nCr(n,r):
    from math import factorial
    f = factorial
    return f(n) / f(r) / f(n-r)    
    
def EngineTemp(Samples, ValueType, ListLen):
    
    Results={}
    if ValueType == 'single':
        for key in Samples:
            Results[key] = Samples[key]['nu1'] * Samples[key]['Sig_f1'] / (Samples[key]['Sig_f1']+Samples[key]['Sig_c1'])
    elif ValueType=='list':
        for key in Samples:
            val=Samples[key]['nu1'] * Samples[key]['Sig_f1'] / (Samples[key]['Sig_f1']+Samples[key]['Sig_c1'])
            Results[key] = [val,val*np.ones(ListLen,float)]
    
    return Results 