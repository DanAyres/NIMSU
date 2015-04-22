'''
Created on 9 Feb 2015

@author: daniel
'''


import numpy as np
from UserInterface.Errors_and_Warnings import ProcessingError


def Decompose_Covariance_Matrix(Covariance_Matrix, depend, decomp, NParams, kl_theta, output):


    if depend=='no':
        
        if Covariance_Matrix.shape != (NParams,NParams):
            msg="Shape of Covariance_Matrix should be (%i, %i)\n Not %s"%(NParams,NParams (Covariance_Matrix.shape))
            ProcessingError(output, "Sampling.Decompose.Decompose_Covariance_Matrix", msg)
            raise ValueError
    
        stoch_dim = NParams
        Sampling_Matrix = np.sqrt(Covariance_Matrix)
    
    elif depend=='yes':
        
        if Covariance_Matrix.shape != (NParams,1):
            msg="Shape of Covariance_Matrix should be (%i,1)\n Not %s"%(NParams, (Covariance_Matrix.shape))
            ProcessingError(output, "Sampling.Decompose.Decompose_Covariance_Matrix", msg)
            raise ValueError
        
        stoch_dim = 1    
        Sampling_Matrix = np.sqrt(Covariance_Matrix) 
            
    elif depend=='corr':
            
        if decomp==None:
            msg="No decomposition specified!"
            ProcessingError(output, "Sampling.Decompose.Decompose_Covariance_Matrix", msg)
            raise ValueError
        elif decomp=='kl':
            Sampling_Matrix, stoch_dim = Karhunen_Loeve(Covariance_Matrix, NParams, kl_theta)
        elif decomp=="chol":
            stoch_dim = NParams
            Sampling_Matrix=Cholesky(Covariance_Matrix)
        else:
            msg="Decomposition %s not recognised"%decomp
            ProcessingError(output, "Sampling.Decompose.Decompose_Covariance_Matrix", msg)
            raise ValueError
            
    return Sampling_Matrix, stoch_dim



def Cumulative_Eigenvalues(eig_vals,tol):
    """
        Return the number (i) of eigenvalues that satisfy
        sum(eig_vals_i) >= tol
        
        Args:
            eig_vals: Vector of eigen values
            tol:      Relative contribution of i eigenvalues
            
        Returns
            i    :  number of eigenvalues
            
        Raises:
        
    """    
    eig_vals = eig_vals/sum(eig_vals)
    eig_sum=0.0
    for i in range(len(eig_vals)):
        eig_sum+=eig_vals[i]
        if(eig_sum>=tol): break
    return i+1



def Cholesky(C):
    """
        Compute the Cholesky factorisation of 
        a matrix.
        
        @b Arguments
        @arg C: A symmetric, positive definate matrix.
        
        @b Returns 
        @arg L: The lower Cholesky factorization of C
        
        @b Raises
        @arg numpy.linalg.LinAlgError

        @b Description
        
        Compute the lower Cholesky factorization of a symmetric, positive definite 
        matrix C which is defined as follows
        
        @f[ C = LL^T 
        @f]
    
    """
    
    try:
        return np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        print "Cholesky error"
        
        
        

def Karhunen_Loeve(C, NParams, kl_theta):
    """
        Form the Karhunen Loeve decomposition
        of a matrix.
        
        Args:
            C          : Covariance matrix
            NParams    : Size of C
            kl_theta   : 
        
        Returns:
            Sampling_Matrix    : KL decomposition of C
            stoch_dim          : Number of stochastic dimensions
        
        Raises:
    
    """            
    
    # KL decomp
    eigs,U= np.linalg.eig(C)
    
    # Determine number of dimensions that satisfy kl_theta
    stoch_dim=Cumulative_Eigenvalues(eigs, kl_theta)
    
    # Form the sampling matrix
    eigs = np.sqrt(eigs)
    Sampling_Matrix=np.zeros([NParams,stoch_dim])
    for g1 in range(stoch_dim):
        Sampling_Matrix[:,g1] = eigs[g1]*U[:,g1]
            
    return Sampling_Matrix, stoch_dim
