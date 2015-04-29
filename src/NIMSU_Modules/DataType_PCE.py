'''
Created on 27 Apr 2015

@author: daiel
'''

from NIMSU_Modules.OrthPol_Interface import OrthPoly
from SparseGrids.SparseGridBaseClass import sparse_grid
import numpy as np

class PCE_Element():

    def __init__(self, ME_Dims, Total_Stoch_Dim, Domain, Inf_Bounds, Pmax, Rules, SG_LevelMax, ValueType):
        """
            Defines a polynomial chaos expansion (PCE) over 
            an 'element' of probability space.
            
            Args:
                ME_Dims:    the dimensions which are to be split into multiple elements
                Total_Stoch_Dim: the total number of stochastic dimensions
                Domain:    the domain of this element in probability space. For example;
                            [ [-1,1], [-1,1], [-1,0] ]  for a three dimensional element.
                Inf_Bounds: specifies if any of the domain boundaries are infinite. Example input;
                            [ [False,False], [False,False], [True,False] ]
                            
                Pmax:    Maximum polynomial order
                Rules:    The quadrature rules for each dimension as define in LibSG. This
                          is used to determine the distribution type for each dimension.
                SG_LevelMax: The maximum level to be used for the sparse grid quadratures
                ValueType:    The type of the results, single, list, hdf5
        """
        
        n_coeff = Pmax + 1
        method='sti'
        n_unions = 4
        dist_args=[0.0,0.0]
        scale=1.0
        
        R0=1.0  # Nominal value - not needed for non-adaptive grids
        
        
        # Calculate the quadrature scheme
        ME_Polys=[]
        if ME_Dims==None:
            # If there are no Multi-Element dimensions then return a sparse grid object over the 
            # whole stochastic space.
            sg=sparse_grid(Rules, \
                        Total_Stoch_Dim, \
                        SG_LevelMax, \
                        R0, \
                        Adaptive=False, \
                        AdaptTol=1E-5, \
                        MaxVals=40, \
                        ValueType=ValueType)
            
            self.Points=sg.Points
            self.Weights=sg.Weights
            self.Num_Points=sg.Num_Points
        
             
        else:
            
            # Create a set of orthogonal polynomials for each 'Multi-Element'
            for D,B,M in zip(Domain,Inf_Bounds,ME_Dims):
                if M in [1,2,3,4]:
                    dist='uniform'
                elif M in [2,6,10]:
                    dist='normal'
                elif M in [7,8]:
                    dist = 'gamma'
                elif M==9:
                    dist='beta'
                else:
                    print 'Rule not recognised'
                ME_Polys.append(OrthPoly(n_coeff, method, n_unions, D, B[0], B[1], dist, dist_args, scale))
        
        
            # Generate a sparse grid for all dimensions which are not going to be split
            ME_Rules = [r for i,r in enumerate(Rules) if i not in ME_Dims]
            sg=sparse_grid(ME_Rules, \
                        Total_Stoch_Dim-len(ME_Dims), \
                        SG_LevelMax, \
                        R0, \
                        Adaptive=False, \
                        AdaptTol=1E-5, \
                        MaxVals=40, \
                        ValueType=ValueType)

            #
            #    The final quadrature set will be a tensor product of the sparse grid and
            #    each of the one dimensional quadratures in the newly generated 
            #    orthogonal polynomials.
            #

            # Calculate the total number of points
            N_ME_Points = 1
            for p in ME_Polys:
                N_ME_Points*=len(p.weights)
            
            self.Num_Points = N_ME_Points * sg.Num_Points
            self.Points = np.zeros([Total_Stoch_Dim,self.Num_Points],float)
            self.Weights = np.zeros([self.Num_Points],float)
            
            
            