'''
Created on 9 Feb 2015

@author: daniel
'''

import numpy as np

def Transform_Statistics(Covariance_Matrix, Corr_Matrix, Sigma, Means, NParams, dist, depend):
    """
        Transform distribution X to Gaussian
    
        Args:
            Covariance_Matrix
            Corr_Matrix
            Sigma       : Standard deviations
            Means
            NParams     : Number of parameters
            dist        : Statistical distribution
            depend      : correlation dependency 
            
        Returns:
            Covariance_Matrix
            Corr_Matrix
            Sigma
            Means
            
        Raises:
    """

    if dist=='normal':
        return Covariance_Matrix, Corr_Matrix, Sigma, Means
    elif dist=='lognormal':
        return Transform_LogNormal(Covariance_Matrix, Corr_Matrix, Sigma, Means, NParams, depend)
    elif dist=='uniform':
        if depend=='corr':
            print "Correlated sampling not support for uniform yet"
            raise ValueError
        else:
            return Covariance_Matrix, Corr_Matrix, Sigma, Means
    else:
        print 'Distribution not recognised'
        

def Transform_LogNormal(Covariance_Matrix, Sigma, Means, NParams, depend):
    """
        Transform LogNormal statistics to Gaussian Statistics
        
        Args:
            Covariance_Matrix
            Corr_Matrix
            Sigma
            Means
            NParams
            
        Returns:
            N_Cov
            N_Corr
            N_Sigma
            N_Means
            
        Raises:
    """
    
    N_Means = LN2N_Mean(Means, Sigma)
    
    if depend == 'corr':
        N_Cov, N_Corr, N_Sigma = LN2N_Covar(Means, Covariance_Matrix, NParams)
        
    elif depend=='no':
        N_Corr = None
        N_Cov = np.zeros([NParams,NParams],float)
        for i in range(NParams):
            N_Cov[i,i] = LogCovar(Sigma[i]**2, Means[i], Means[i])
        N_Sigma = np.sqrt(N_Cov.diagonal())
        
    elif depend=='yes':
        N_Corr = None
        N_Cov = np.zeros([NParams,1],float)
        for i in range(NParams):
            N_Cov[i,0] = LogCovar(Sigma[i]**2, Means[i], Means[i])
        N_Sigma = np.sqrt(N_Cov[:,0])
        
    return N_Cov, N_Corr, N_Sigma, N_Means


def LN2N_Mean(LN_Mean, LN_SD):
    """
        Calculate the Normally distributed equivalent mean
        from the Log-Normally distributed mean and standard deviation
    
        Args:
            LN_Mean    : Log-Normal mean - (numpy vector)
            LN_SD      : Log-Normal standard deviation - (numpy vector)
            
        Returns:
            N_Mean     : Normally distributed mean - (numpy vector)
    """
    
    return  np.log(  LN_Mean**2 / np.sqrt( LN_SD**2 + LN_Mean**2  )   )

def LN2N_Covar(LN_Mean, LN_Covar, NParams):
    """
        Calculate the Normally distributed equivalent covariance
        from the Log-Normally distributed mean and standard deviation
       
       Args:
           LN_Mean    : Log-Normal mean - (numpy vector)
           LN_Covar   : Log_Normal covaraince matrix - (2D numpy vector)
           
        Returns:
            N_Covar   : Normally distributed covariance
            N_Corr    : Normally distributed correlation
            N_SD      : Normally distributed standard deviations
    """
    
    N_Covar=np.zeros([NParams,NParams],float)
    
    for i in range(NParams):
        for j in range(i+1):
            N_Covar[i,j] = LogCovar( LN_Covar[i,j], LN_Mean[i], LN_Mean[j] )
            N_Covar[j,i] = N_Covar[i,j]
    
    N_SD = np.sqrt(N_Covar.diagonal())

    N_Corr=FormCorrelation(N_Covar, N_SD, NParams)
    
    return N_Covar, N_Corr, N_SD


def LogCovar(LN_Covar, LN_Mean_i, LN_Mean_j):
    return np.log( LN_Covar / ( LN_Mean_i*LN_Mean_j ) + 1.0 )

def FormCorrelation(Covariance_Matrix, Sigma, NParams):
    """
    
    """
    Corr_Matrix=np.zeros([NParams,NParams],float)
    for i in range(NParams):
        for j in range(i+1):
            Corr_Matrix[i,j] = Covariance_Matrix[i,j]/(Sigma[i]*Sigma[j])
            Corr_Matrix[j,i] = Corr_Matrix[i,j]
            
    return Corr_Matrix