'''
Created on 9 Feb 2015

@author: daniel
'''

import numpy as np
from UserInterface.Errors_and_Warnings import ProcessingError, WarningMessage



def Form_Covariance_Matrix(Corr_Matrix, Sigma, NParams, depend, kernel, kernel_args, output):
    """
        From the user-specified parameters, form a covariance matrix.
        
        If no correlation matrix is specified - generate one using the
        specified kernel.
        
        Test for positive-definiteness. If not, force
        negative eigenvalues = 0.
        
        Args:
        
            output        : Location to direct any Errors/Warnings
        Returns:
        
        
        Raises:
    
    """
    
    
    
    # All parameters in group are independent
    if depend=='no':
        Covariance_Matrix=np.zeros([NParams,NParams],float)
        np.fill_diagonal(Covariance_Matrix, Sigma**2)
                
    # All parameters are dependent
    elif depend=='yes':
        Covariance_Matrix=np.zeros([NParams,1],float)
        Covariance_Matrix[:,0] = Sigma * Sigma

    elif depend=='corr':
        # Generate correlation matrix using group.kernel
        if Corr_Matrix==None:
            if kernel=='exponential':
                Corr_Matrix = Exponential_kernel(NParams, kernel_args)
            else:
                msg="Kernel '%s' not recognised"%kernel
                ProcessingError(output, "Form_Covariance_Matrix", msg)
                raise ValueError
        
        
        # Form Covariance matrix from the correlation matrix
        Covariance_Matrix=Calculate_Covariance_Matrix(Corr_Matrix, Sigma, NParams)
        
        
        # Test for positive-definite. If not, then force all eigvals > 0.            
        Covariance_Matrix, flag=PosDef(Covariance_Matrix, NParams, output)
    
        if flag:
            print 'Need to write modified covariance to file'
    
    
    return Covariance_Matrix


def Calculate_Covariance_Matrix(Correlation_Matrix, Sigma, NParams):
    """
        Calculate the covariance from the correlation.
        
        Args:
            Correlation_Matrix
            Sigma                : Vector of standard deviations
            NParams              : Size of covariance matrix
            
        Returns:
            Covariance_Matrix
            
        Raises:
    
    """
    Covariance_Matrix=np.zeros([NParams,NParams],float)
    for i in range(NParams):
        for j in range(i+1): 
            Covariance_Matrix[i,j] = Correlation_Matrix[i,j]*Sigma[i]*Sigma[j]
            if j!=i: Covariance_Matrix[j,i] = Covariance_Matrix[i,j]
            
    return Covariance_Matrix


def Exponential_kernel(NParams, cor_len):
    """
        Exponential correlation kernel.
        
        rho_ij = exp( -|i-j|/cor_len )
        
        Args:
            NParams    : Number of parameters in correlation matrix
            cor_len    : Correlation length
        
        Returns
            Corr_Matrix
            
        Raises:
        
    """
    Corr_Matrix=np.zeros([NParams,NParams])
    for i in range(NParams):
        for j in range(i+1):
            Corr_Matrix[i,j] = np.exp(-np.fabs(i-j)/cor_len)
            Corr_Matrix[j,i] = Corr_Matrix[i,j]
    return Corr_Matrix


def PosDef(Covariance_Matrix, NParams, output):
    """
        Test if a matrix is positive definite.
        If not, then set all negative eigenvalues
        to zero and recombine the matrix.
        
        Args:
            Covariance_Matrix
            NParams
            output        : Location to direct any Errors/Warnings
        
        Returns:
            Covariance_Matrix
            flag                : True if not PosDef and Covariance_Matrix
                                  matrix has been modified
            
        Raises:
    """
    
    flag=False
    eigs,R= np.linalg.eig(Covariance_Matrix)
    
    temp=np.zeros([NParams,NParams],float)
    for i in range(NParams):
        if eigs[i] < 0.0:
            flag=True
            msg="Eigenvalue %i is less than zero"%(i+1)
            WarningMessage(output, "PosDef", msg) 
            temp[i,i] = 0.0
        else: temp[i,i] = eigs[i]
    #print eigs
    return R.dot(temp).dot(R.T), flag
    
    
