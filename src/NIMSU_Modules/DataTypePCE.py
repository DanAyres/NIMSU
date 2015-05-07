'''
Created on 27 Apr 2015

@author: Daniel Ayres

References:

-   Walter Gautschi
    On Generating Orthogonal Polynomials
    SIAM J Sci. Stat. Comput.
    Vol 3, No 3, pp 289-317 
    1982 

-   Xiaoliang Wan and George Karniadakis
    Multi-Element Generalized Polynomial Chaos for Arbitrary Probability Measures
    SIAM J Sci. Stat. Comput.
    Vol 28, No 3, pp 901-928
    2006

'''

from NIMSU_Modules.OrthPol_Interface import OrthPoly
from NIMSU_Modules.SparseGridBaseClass import sparse_grid, CalcSparseSet
from scipy.sparse import dia_matrix
import numpy as np

import itertools
from operator import mul


class PCE():
    
    def __init__(self, control_info, sampler, Nominal, 
                 Pmin=0, 
                 Dims=None, 
                 QuadIdx=None,
                 CalcQuadrature=True):
        """
        
            This data type defines a polynomial chaos expansion over
            a specified number of dimensions. Some, or all, of these
            dimensions may be designated "adaptive". In which case,
            the polynomial basis for these dimensions is constructed 
            using the Gram-Schmidt procedure.
            
            
            The initial element structure is a single element over the 
            entire probability space.
            
            Calling the adapt() method computes an error measure for each 
            adaptive dimension and splits the elements accordingly.
            The default splitting is by 2. 
             
            
            
            Args:
            
                Dims:    list.
                         This is used with the HDMR method. It contains the global dimensions
                         associated with this PCE. For example Dims = [1,10,12].
                         
                NDim:    int.
                         The total number of dimensions in the PCE.
                          
                Pmax:    int.
                         The maximum polynomial order.
                         
                Rules:   list.
                         The quadrature rules for each dimension as defined in LibSG. This
                         is used to determine the distribution type for each dimension as well
                         as the polynomial coefficients and norms for non-adaptive dimensions.
            
            Optional Args:
            
                ME_Adapt_Dims:    The members of Dims which are to be adaptively split
                
                Domain:    the domain of this element in probability space. For example;
                            [ [-1,1], [-1,1], [-1,0] ]  for a three dimensional element.
                Inf_Bounds: specifies if any of the domain boundaries are infinite. Example input;
                            [ [False,False], [False,False], [True,False] ]
                            
                Pmin:    Minimum polynomial order
                
                SG_LevelMax: The maximum level to be used for the sparse grid quadratures
                
                ValueType:    The type of the results, single, list, hdf5
        
        """
        
        self.control_info=control_info
        self.sampler=sampler
        self.Nominal=Nominal
        
        self.Dims=Dims
        self.CalcQuadrature=CalcQuadrature
        
        self.Elements=[]
        
        if self.control_info.scm_me_dims is not None:
            InfBoundaries=[]
            Domains=[]
            Distributions=[]
            
            # Assign the distributions, boundaries and domains for each adaptive dimension
            for m in self.control_info.scm_me_dims:
                if self.sampler.rules[m] in [1,2,3,4]:
                    Distributions.append('uniform')
                    Domains.append([-1.0,1.0])
                    InfBoundaries.append([False,False])
                elif self.sampler.rules[m] in [5,6,10]:
                    Distributions.append('normal')
                    Domains.append([-1.0,1.0])
                    InfBoundaries.append([True,True])
                elif self.sampler.rules[m] in [7,8]:
                    Distributions.append('gamma')
                    Domains.append([0.0,1.0])
                    InfBoundaries.append([False,True])
                elif self.sampler.rules[m]==9:
                    Distributions.append('beta')
                    Domains.append([-1.0,1.0])
                    InfBoundaries.append([False,False])
                else:
                    raise ValueError("Quadrature rule not recognised")
        
            
            self.Elements.append(PCE_Element(self.control_info, self.sampler, self.Nominal,
                             Pmin,
                             Distributions,
                             Domains,
                             InfBoundaries,
                             QuadIdx, 
                             self.Dims, 
                             CalcQuadrature))
            
        else:

            Distributions=None
            Domains=None
            InfBoundaries=None
            self.Elements.append(PCE_Element(self.control_info, self.sampler, self.Nominal,
                             Pmin,
                             Distributions,
                             Domains,
                             InfBoundaries,
                             QuadIdx, 
                             self.Dims, 
                             CalcQuadrature))
                             

        # Calculate the probability density function for this PCE
        self.PDF=1.0
        for i,r in enumerate(self.sampler.rules):
            
            # These dimensions are not present in the current (local) expansion
            if self.Dims is not None:
                if i not in self.Dims: continue
            else:
                if r in [1,2,3,4]:
                    self.PDF*=0.5
                elif r in [5,6,10]:
                    self.PDF/=np.sqrt(np.pi)
                elif r in [7,8]:
                    print 'Not implemented'
                    raise
                elif r==9:
                    print 'Not implemented'
                    raise    
                else:
                    raise ValueError("Quadrature rule not recognised")        
                    
    
    
    def CalculateSamples(self,):
        """
            Loop through each of the PCE Elements and generate samples
            using their quadrature points.
            
            Args:
            
            Returns:
                Samples:    dict
        """
        
        
        pos=0
        Samples={}
        for ele in self.Elements:
            
            if ele.control_info.scm_quad_adapt:
                if ele.AdaptiveQuadFlag==False: continue
                Point_Locs=ele.SparseGrid.ComputeIndex
            else:
                Point_Locs = range(ele.SparseGrid.Num_Points)
                
            for p in Point_Locs:
                key='s'+str(pos)
                Samples[key] = self.sampler.sample(ele.SparseGrid.Points[:,p])
                pos+=1

        return Samples

    def UpdateIntegrals(self, Results):
        """
            Add the computed values to the sparse grid object
            for each element. If the quadrature is adaptive:
            call the sparse_grid.adapt method.
            
            Args:
                
                Results:    dict
                
            Returns:    
            
                flag:    logical
                         True if the adaptive integration is complete
                         
        """
        
        pos=0
        for ele in self.Elements:
            
            if ele.control_info.scm_quad_adapt:
                if ele.AdaptiveQuadFlag==False: continue
                values=[]
                for dummy in ele.SparseGrid.ComputeIndex:
                    key='s'+str(pos)
                    pos+=1
                    values.append(Results[key]) 
                ele.AdaptiveQuadFlag=ele.SparseGrid.Adapt(values)
            else:
                ele.SparseGrid.Values=[]
                for dummy in range(ele.SparseGrid.Num_Points):
                    key='s'+str(pos)
                    pos+=1
                    ele.SparseGrid.Values.append(Results[key])
        
        return all([ele.AdaptiveQuadFlag==False for ele in self.Elements])

    def CalculateErrorMeasures(self):
        
        for ele in self.Elements:
            for poly in ele.pce:
                print poly[1:]
        
        pass

    def SplitElements(self):
        
        return all([ele.AdaptiveElementFlag==False for ele in self.Elements])


class PCE_Element():
    
    def __init__(self,  control_info, sampler, Nominal,
                             Pmin,
                             Distributions,
                             Domains,
                             InfBoundaries,
                             QuadIdx, 
                             Dims, 
                             CalcQuadrature): 
                        
        
        
        self.control_info=control_info
        self.sampler=sampler
        self.Nominal=Nominal
        
        # ME parameters
        self.Distributions=Distributions
        self.Domains=Domains
        self.InfBoundaries=InfBoundaries
        
        # Optional HDMR parameters
        self.Pmin=Pmin
        self.QuadIdx=QuadIdx
        self.Dims=Dims
        self.CalcQuadrature=CalcQuadrature
                
        
        self.EleProb=1.0
        
        if self.control_info.scm_me_dims==[]: self.control_info.scm_me_dims=None
        
        if self.control_info.scm_quad_adapt: self.AdaptiveQuadFlag=True
        else:self.AdaptiveQuadFlag=False
        
        if self.control_info.scm_me_dims==None: self.AdaptiveElementFlag=False
        else:self.AdaptiveElementFlag=True

        # Some hard-wired parameters for the generation of orthogonal polynomials
        self.Orthpol_method='sti'
        self.Orthpol_n_unions = 4
        self.Orthpol_dist_args=[0.0,0.0]
        self.Orthpol_scale=1.0


        self.ME_Polys=[]
        
        if self.control_info.scm_me_dims is not None:
            
            # Compute the quadrature rules for each of the adaptive dimensions
            for D,B,dist in zip(self.Domains,self.InfBoundaries,self.Distributions):
                self.ME_Polys.append(OrthPoly(self.control_info.scm_poly_order+1, 
                                         self.Orthpol_method, 
                                         self.Orthpol_n_unions, 
                                         D, B[0], B[1], dist, 
                                         self.Orthpol_dist_args, 
                                         self.Orthpol_scale))
            
            # Calculate the element probability
            self.EleProb = reduce(mul,[sum(p.Weights)*p.PDF for p in self.ME_Polys])
            
            # Generate a sparse grid for all dimensions which are not going to be split
            # This is a non adaptive grid. Tensor Product may also be specified.
            if len(self.control_info.scm_me_dims) != self.sampler.Total_stoch_dim:
                ME_Rules = [r for i,r in enumerate(self.sampler.rules) if i not in self.control_info.scm_me_dims]
                sg=sparse_grid(ME_Rules, \
                            self.sampler.Total_stoch_dim-len(self.control_info.scm_me_dims), \
                            self.control_info.scm_sg_level[0], \
                            self.control_info.ValueType)
            else:
                sg=None

            # Form the full quadrature scheme
            self.Points, self.Weights, Idx, self.Num_Points = ME_Quadrature(sg, self.ME_Polys, self.sampler.Total_stoch_dim, self.control_info.scm_poly_order, self.control_info.scm_me_dims)
            
            #Build a PCE from the quadrature index
            self.pce=IndextoBasis(Idx, self.sampler.Total_stoch_dim, self.control_info.scm_poly_order)
            
        
        else:

            # No quadrature just PCE
            if self.CalcQuadrature==False: 
                        
                # Build PCE using the usual method
                if self.QuadIdx==None:
                    self.pce=PceBasis(self.Pmin, self.control_info.scm_poly_order, self.sampler.Total_stoch_dim)
                
                # Compute a PCE with a user-specified quadrature index from a sparse grid
                else:
                    self.pce = IndextoBasis(self.QuadIdx, self.sampler.Total_stoch_dim, self.control_info.scm_poly_order)
                    
            
            # If there are no Multi-Element dimensions then return a sparse grid object over the 
            # whole stochastic space.
            else:
                if self.sampler.rules is None:
                    raise ValueError("Quadrature rules must be known when calculating quadratures")
                
                self.SparseGrid=sparse_grid(self.sampler.rules, \
                                            self.sampler.Total_stoch_dim, \
                                            self.control_info.scm_sg_level[0], \
                                            self.control_info.ValueType,\
                                            Nominal=self.Nominal,\
                                            Adaptive=self.control_info.scm_quad_adapt, \
                                            AdaptTol=self.control_info.scm_quad_tol, \
                                            MaxVals=40 )
        
                self.pce = IndextoBasis(self.SparseGrid.Idx, self.sampler.Total_stoch_dim, self.control_info.scm_poly_order)
                
                
        # Assign the polynomial coefficients and norms to each dimension.
        self.Coeffs,self.Norms = AssignCoeffsNorms(self.control_info.scm_poly_order, self.Dims, self.control_info.scm_me_dims, self.ME_Polys, self.sampler.rules)
        
    def PolyNorms(self, orders, dims):
        """
            Compute the product of the one dimensional polynomial normalisation
            coefficients. Only those multiplications are performed for the 
            dimensions specified by the list <dims>.
            
            Args:
            
                orders:    list
                           the one dimensional polynomial orders.
                           
                dims:      list
                           the dimensions to be multiplied.
           
           
            Returns:
                
                Product of one dimensional norms.
        """
        
        return reduce(mul, [ N[o] for N,o in zip(  [self.Norms[d] for d in dims]  ,orders) ]  )

    def PolyProd(self, orders, x, dims):
        """
            Evaluate the multidimensional polynomial P(x_1,x_2,...x_n).
            The multidimensional polynomial is a product of one dimensional
            polynomials  P(x_1,x_2,...x_n) = p(x_1)p(x_2)...p(x_n). 
            The coefficients (powers of x) of each polynomial are
            stored in self.Coeff.
            
            
            Args:
                orders:    list
                           one dimensional polynomial orders
                           
                x:         list
                           the value of x in each dimension
        
        
            Returns:
                
                prod:     float
                          Product of 1D polynomials of order <orders> evaluated at <x>
        """
        
        
        prod=1.0
        for o,val,d in zip(orders,x,dims):
            rv=np.power(val,range(self.control_info.scm_poly_order+1))
            P=self.Coeffs[d].dot(rv)
            prod*=P[o]
        return prod


def ME_Quadrature(sg, ME_Polys, Total_Stoch_Dim, Pmax, ME_Adapt_Dims):
    """
        The final quadrature set will be a Cartesian product of the sparse grid and
        each of the one dimensional quadratures in the newly generated 
        orthogonal polynomials.
        
           
        Args:
           
            sg:    sparse_grid
                    A sparse grid object defined over all dimensions which are not adaptive.
                    
            ME_Polys:   list
                        A list of orth_pol objects. One for each adaptive dimension.
                        
            Total_Stoch_Dim:    int
                                Total number of dimensions
            
            Pmax:    int
                     maximum polynomial order
                     
            ME_Adapt_Dims:    list[int]
                              A list of the adaptive multi-element dimensions.
                              
           
        Returns:
        
            Points:    ndarray[Total_Stoch_Dim, Num_Points]
                       The quadrature points.
            
            Weights:   ndarray[Num_Points]
                       The quadrature weights.
                       
            Idx:       ndarray[Total_Stoch_Dim, :]
                       The indexs associated with the construction of the quadrature grid.
            
            Num_Points:    int
                           Number of points in the quadrature rule.
        
        """ 
    if sg is None:
        SG_Num_Points = 1
        SG_Idx = np.zeros([1,1],float) 
        SG_Weights = 1.0
        SG_Points = 1.0
    else:
        SG_Num_Points = sg.Num_Points
        SG_Weights = sg.Weights
        SG_Points = sg.Points
        SG_Idx = sg.Idx
                
    # Calculate the total number of points
    N_ME_Points = 1
    for p in ME_Polys:
        N_ME_Points*=len(p.Weights)
    Num_Points = N_ME_Points * SG_Num_Points

    # Calculate the quadrature scheme
    Points = np.zeros([Total_Stoch_Dim, Num_Points],float)
    Weights = np.zeros([Num_Points],float)
    lenidx=len(SG_Idx[0,:])*np.power(Pmax,len(ME_Adapt_Dims))
    Idx = np.zeros( [Total_Stoch_Dim,lenidx] ,int)
    Not_ME_dims=[i for i in range(Total_Stoch_Dim) if i not in ME_Adapt_Dims]
            
    # Cartesian product
    # Points
    pts=[p.Points for p in ME_Polys]
    tmp=itertools.product(*pts)
    for counter,j in enumerate(tmp):
        if len(ME_Adapt_Dims) != Total_Stoch_Dim:
            Points[Not_ME_dims,counter*SG_Num_Points:(counter+1)*SG_Num_Points] = SG_Points            
        Points[ME_Adapt_Dims,counter*SG_Num_Points:(counter+1)*SG_Num_Points] = np.asarray(j).reshape((len(ME_Adapt_Dims),1))

    # Weights
    wts=[p.Weights for p in ME_Polys]
    tmp=itertools.product(*wts)
    for counter,j in enumerate(tmp):
        Weights[counter*SG_Num_Points:(counter+1)*SG_Num_Points] = SG_Weights*reduce(mul,j)
            
    # Index
    idx=[range(Pmax) for p in ME_Polys]
    tmp=itertools.product(*idx)
    lenidx=len(SG_Idx[0,:])
    for counter,j in enumerate(tmp):
        if len(ME_Adapt_Dims) != Total_Stoch_Dim: 
            Idx[Not_ME_dims,counter*lenidx:(counter+1)*lenidx]=sg.Idx
        Idx[ME_Adapt_Dims,counter*lenidx:(counter+1)*lenidx] = np.asarray(j).reshape((len(ME_Adapt_Dims),1))

    return Points, Weights, Idx, Num_Points

def AssignCoeffsNorms(Pmax, Dims, ME_Adapt_Dims, ME_Polys, Rules):
    """
        For each dimension in a PCE_Element, we assign the norms and
        polynomial coefficients up to a maximum order of Pmax.
            
        The coefficients are those of the powers of x in the polynomials
        P(x). Realisations of the first Pmax polynomials may be computed
        as follows
            
        a=numpy.power(x,range(Pmax+1))
        P = Coeff.dot(x)
            
        The norms are defined as follows
            
        N_i  = int P_i(x)^2 w(x) dx
            
        where w(x) is the appropriate weighting function.
        
        
        Args:
                
            Pmax:    Int.
                    The maximum polynomial order
                        
        Returns:
            
            Coeffs:   list.
                        The coefficients of the powers of x for each polynomial.
                        Each list element is a numpy matrix.
                          
            Norms:    list
                        The normalisation coefficient for each polynomial.
                          
        Raise:
            
            ValueError:    If the quadrature rule is not recognised
                            
        
    """
        
    uni=Legendre_Coefficients(Pmax)
    uni_n=[LN(order) for order in range(Pmax+1)]
    gau=Hermite_Coefficients(Pmax)
    gau_n=[HN(order) for order in range(Pmax+1)]
        
    Coeffs=[]
    Norms=[]
    ctr=0
    for i,r in enumerate(Rules):
            
        # These dimensions are not present in the current (local) expansion
        if Dims is not None:
            if i not in Dims: continue
            
        # Adaptive dimensions have their own quadrature scheme
        if ME_Adapt_Dims is not None:
            if i in ME_Adapt_Dims:
                Coeffs.append(ME_Polys[ctr].coefficients*ME_Polys[ctr].PDF)
                Norms.append(ME_Polys[ctr].norms*ME_Polys[ctr].PDF)
                continue
                
        # Add the conventional     
        if r in [1,2,3,4]:
            Coeffs.append(uni)
            Norms.append(uni_n)
        elif r in [5,6,10]:
            Coeffs.append(gau)
            Norms.append(gau_n)
        elif r in [7,8]:
            print 'Not implemented'
            raise
        elif r==9:
            print 'Not implemented'
            raise    
        else:
            raise ValueError("Quadrature rule not recognised")        
        
    return Coeffs, Norms

def PceBasis(pmin, pmax, total_dim):
    """
        Compute a polynomial basis in multi-dimensions of the form
        
        P(x) = p_1(x_1)p_x(x_2)...p_n(x_n)
        
        where n=total_dim.
        
        The basis functions are restricted as follows:
        
          pmin <= sum( order(p_i) ) <= pmax
          
          
        Args:
            
            pmin:   int
                    The minimum total order of the basis
                    
            pmax:   int
                    The maximum total order of the basis
                    
            total_dim:  int
                        The total number of dimensions in the basis.
        
        Returns:
        
            pce:    list
                    The pce list has the folowing form:
                    [ coefficient_value, orders, dims ]
                    
                    coefficient_value is the value of the PCE coefficient. Initialised to zero
                    orders are the non-zero orders of the 1D polynomials.
                    dims are the dimensions corresponding to the non-zero orders.
        
    """
    
    pce=[]
    
    basis=CalcSparseSet(pmin, pmax, total_dim)
    N_poly=np.shape(basis)[1]
    for i in range(N_poly):
        poly=[val for val in basis[:,i] if val !=0 ]
        dims=[val for val in range(total_dim) if basis[val,i] !=0 ]
        if dims==[]:pce.append([0.0, [0], [0]])
        else: pce.append([0.0, poly, dims])
    
    return pce

def IndextoBasis(levels, dims, Pmax):
    """
        Convert the levels from a generalised quadrature scheme
        to a polynomial chaos basis.
            
    """
        
    # Create the PCE orders from the quadrature levels
    levels=Quad2Poly(levels, dims) 

        
    # Create the PCE list
    pce_list=[]        
    for lev in levels:
            
        temp_dim=[]
        temp_lev=[]
        for l,d in zip(lev,range(dims)):
        
            if l==0: continue
            temp_dim.append(d)
            temp_lev.append(l)
            
        if any([p>Pmax for p in temp_lev]): continue
        
        if temp_lev==[]: pce_list.append([0.0,[0],[0]])
#         if temp_lev==[]: pce_list.append([0,0.0,[0],[0]])
        else:
#             pce_list.append([max(temp_lev),0.0, temp_lev, temp_dim])
            pce_list.append([0.0, temp_lev, temp_dim])

    return pce_list

def Quad2Poly(levels, stoch_dim):
    from copy import deepcopy as copy
        
    levels=map(list, zip(*levels))

    flevels=copy(levels)
    for lev in flevels:

        """ Calculate the forward neighbours """
        kk=copy(lev) # Be careful with bloody pointers.
        for j in range(stoch_dim):
            if kk[j]==0:continue
#             if kk[j]==1:continue
            if kk[j]==1 or kk[j]==2:
                kk[j]+=1
                if (kk not in levels):
                    if Admissible(kk, levels, stoch_dim): 
                        levels.append(copy(kk))
                kk[j]-=1 
            else:
                iterate= (np.power(2,kk[j]-1) - kk[j] ) + 1
#                 print kk[j], iterate
                for l in range(1,iterate):
                    kk[j]+=l
                    if kk[j]>11: continue
                    if (kk not in levels):
                        if Admissible(kk, levels, stoch_dim): 
                            levels.append(copy(kk))
                    kk[j]-=l
                    
    return levels

def Admissible( k, Active, num_dim):
    admis=[]
    for q in range(num_dim):
        if k[q]==0:
            admis.append(True)
            continue
        k[q]-=1
        if k in Active: admis.append(True)
        else:admis.append(False)
        k[q]+=1
            
    return all(admis)

def Hermite_Coefficients(pmax):
    """
    Return the coefficient matrix for the first pmax
    Hermite ploynomials. Return matrix in diagonal storage
    format.
    """
    
    if pmax>10:
        print "Hermite order %2i not supported"%pmax
        raise ValueError
    
    Herm_Coef=np.zeros([13,13],float)


#   x^i                   x^{i-2}                  x^{i-4}                  x^{i-6}       
    Herm_Coef[0,0]= 1.0
    Herm_Coef[1,1]= 1.0
    Herm_Coef[2,2]= 1.0  ; Herm_Coef[2,0]= -1.0
    Herm_Coef[3,3]= 1.0  ; Herm_Coef[3,1]= -3.0
    Herm_Coef[4,4]= 1.0  ; Herm_Coef[4,2]= -6.0   ; Herm_Coef[4,0]= 3.0
    Herm_Coef[5,5]= 1.0  ; Herm_Coef[5,3]= -10.0  ; Herm_Coef[5,1]= 15.0
    Herm_Coef[6,6]= 1.0  ; Herm_Coef[6,4]= -15.0  ; Herm_Coef[6,2]= 45.0   ; Herm_Coef[6,0]= -15.0
    Herm_Coef[7,7]= 1.0  ; Herm_Coef[7,5]= -21.0  ; Herm_Coef[7,3]= 105.0  ; Herm_Coef[7,1]= -105.0
    Herm_Coef[8,8]  = 1.0  ; Herm_Coef[8,6]  = -28.0  ; Herm_Coef[8,4]  = 210.0  ; Herm_Coef[8,2]  = -420.0     ; Herm_Coef[8,0]  = 105.0
    Herm_Coef[9,9]  = 1.0  ; Herm_Coef[9,7]  = -36.0  ; Herm_Coef[9,5]  = 378.0  ; Herm_Coef[9,3]  = -1260.0    ; Herm_Coef[9,1]  = 945.0
    Herm_Coef[10,10]= 1.0  ; Herm_Coef[10,8] = -45.0  ; Herm_Coef[10,6] = 630.0  ; Herm_Coef[10,4] = -3150.0    ; Herm_Coef[10,2]  = 4725.0 ; Herm_Coef[10,0]  = -945.0
    Herm_Coef[11,11]= 1.0  ; Herm_Coef[11,9]  = -55.0  ; Herm_Coef[11,7] = 990.0   ; Herm_Coef[11,5] = -6930.0     ; Herm_Coef[11,3]  = 17325.0 ; Herm_Coef[11,1]  = -10395.0
    Herm_Coef[12,12]= 1.0  ; Herm_Coef[12,10] = -66.0  ; Herm_Coef[12,8] = 1485.0  ; Herm_Coef[12,6] = -13860.0    ; Herm_Coef[12,4]  = 51975.0 ; Herm_Coef[12,2]  = -62370.0 ; Herm_Coef[12,0]  = 10395.0

    temp=np.zeros([pmax+1,pmax+1],float)
    for i in range(pmax+1):
        temp[i,:] =  Herm_Coef[i,0:pmax+1]  
        
    return dia_matrix(temp)

def Legendre_Coefficients(pmax):
    """
    Return the coefficient matrix for the first pmax
    Legendre ploynomials. Return matrix in diagonal storage
    format.
    """
    
    Leg_Coef=np.zeros([9,9],float)

#   x^i                          x^{i-2}                       x^{i-4}                       x^{i-6}              
    Leg_Coef[0,0]= 1.0
    Leg_Coef[1,1]= 1.0
    Leg_Coef[2,2]= 3.0/2.0          ; Leg_Coef[2,0]= -1.0/2.0
    Leg_Coef[3,3]= 5.0/2.0          ; Leg_Coef[3,1]= -3.0/2.0
    Leg_Coef[4,4]= 35.0/8.0         ; Leg_Coef[4,2]= -30.0/8.0          ; Leg_Coef[4,0]= 3.0/8.0
    Leg_Coef[5,5]= 63.0/8.0         ; Leg_Coef[5,3]= -70.0/8.0          ; Leg_Coef[5,1]= 15.0/8.0
    Leg_Coef[6,6]= 231.0/16.0       ; Leg_Coef[6,4]= -315.0/16.0        ; Leg_Coef[6,2]= 105.0/16.0       ; Leg_Coef[6,0]= -5.0/16.0
    Leg_Coef[7,7]= 429.0/16.0       ; Leg_Coef[7,5]= -693.0/16.0        ; Leg_Coef[7,3]= 315.0/16.0       ; Leg_Coef[7,1]= -35.0/16.0
    Leg_Coef[8,8]= 6435.0/128.0     ; Leg_Coef[8,6]= -12012.0/128.0     ; Leg_Coef[8,4]= 6930.0/128.0     ; Leg_Coef[8,2]= -1260.0/128.0   ; Leg_Coef[8,0]= 35.0/128.0

#    Leg_Coef[9,9]= 12155.0/128.0     ; Leg_Coef[9,7]= -25740.0/128.0     ; Leg_Coef[9,5]= 18018.0/128.0     ; Leg_Coef[9,3]= -4620.0/128.0   ; Leg_Coef[9,1]= 315.0/128.0
#    Leg_Coef[10,10]= 6435.0/256.0     ; Leg_Coef[10,8]= -12012.0/256.0     ; Leg_Coef[10,6]= 6930.0/256.0     ; Leg_Coef[10,4]= -1260.0/256.0   ; Leg_Coef[10,2]= 35.0/256.0    ; Leg_Coef[10,0]= 35.0/256.0
#    Leg_Coef[11,11]= 6435.0/256.0     ; Leg_Coef[11,9]= -12012.0/256.0     ; Leg_Coef[11,7]= 6930.0/256.0     ; Leg_Coef[11,5]= -1260.0/256.0   ; Leg_Coef[11,3]= 35.0/256.0    ; Leg_Coef[11,1]= 35.0/256.0

    pmax1=min(pmax,8)
    temp=np.zeros([pmax1+1,pmax1+1],float)
    for i in range(pmax1+1):
        temp[i,:] =  Leg_Coef[i,0:pmax1+1]  
        
    return dia_matrix(temp)

def HN(order):
    facts=[1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800]
    return facts[order]

def LN(order):
    return float(1.0/(2.0*order+1.0))