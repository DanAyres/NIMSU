'''
Created on 7 Feb 2015

@author: daniel
'''

from Sampling.SampFormCovariance import Form_Covariance_Matrix
from Sampling.SampTransforms import Transform_Statistics
from Sampling.SampDecompose import Decompose_Covariance_Matrix
from Sampling.SampDistributions import Assign_Sampling_Routine, Assign_Rules_PDFS

class SampleUncertainties():

    def __init__(self, input_params, control_info, output):
        """
            @brief Base class used for sampling sets of random variables
        
            @b Description
            Prepares the input data for sampling. This consists of five operations:
            - Formation of a covariance matrix
            - Transformation to Gaussian random variables (if required)
            - Decomposition of the covariance matrix using Karhunen Loeve or Cholesky factorization
            - Assignment of the sampling method
            - Assignment of the quadrature rules and the probabity density functions
            
            @b Arguments
            @arg input_params: 
            @arg control_info
            @arg output
        
        """
    
        Total_Params = 0
        Total_stoch_dim = 0
        self.rules=[]
        self.PDF=[]
        
        for group in input_params:
            
            # Form the covariance matrix from the correlation and SDs or generate using kernel
            Covariance_Matrix = Form_Covariance_Matrix(group.Corr_Matrix, 
                                                       group.Sigma, 
                                                       group.NParams, 
                                                       group.Depend, 
                                                       group.Kernel, 
                                                       group.Kernel_args,
                                                       output)
              
            # Transform to Gaussian space (if necessary)
            Covariance_Matrix, group.C, group.Sigma, group.Means = Transform_Statistics(Covariance_Matrix, 
                                                                                        group.Corr_Matrix, 
                                                                                        group.Sigma, 
                                                                                        group.Means, 
                                                                                        group.NParams, 
                                                                                        group.Dist,
                                                                                        group.Depend)
            
            # Create the sampling matrix
            group.Sampling_Matrix, group.Stoch_dim = Decompose_Covariance_Matrix(Covariance_Matrix, 
                                                                                 group.Depend, 
                                                                                 group.Decomp, 
                                                                                 group.NParams, 
                                                                                 group.Kl_theta,
                                                                                 output)
                
            # Assign the sampling routine
            group.Sample= Assign_Sampling_Routine(group.Dist, group.Depend)
            
            # Assign the quadrature rules and probability density functions
            rules,PDF = Assign_Rules_PDFS(group.Dist, group.Stoch_dim)
            self.rules.extend(rules)
            self.PDF.extend(PDF)
            
         
            Total_Params += group.NParams
            Total_stoch_dim += group.Stoch_dim
            
        self.Total_Params = Total_Params
        self.Total_stoch_dim = Total_stoch_dim
        
        
        self.input_params = input_params
     
    # No return   
    
    
    def sample(self, rvs):
        """
            Generate a set of samples using a vector of unit random variables.
            
            @b Arguments
            @arg rvs: A vector or unit random variables
                
            @b Returns
            @arg  samples : Dictionary of sampled parameters (key = parameter name) .
            
        """
    
        samples={}
        Stoch_Start=0
        Stoch_End=0
                
        for group in self.input_params:
            Stoch_End+=group.Stoch_dim
            values = group.Sample(group.Means, group.Sampling_Matrix, rvs[Stoch_Start:Stoch_End])
            for pos,name in enumerate(group.Names):
                samples[name]=values[pos]
            
            Stoch_Start+=group.Stoch_dim
                        
        return samples


#def sample_withrules(rvs,input_params):
#    """
#        Generate a set of samples using the array of variables rvs. 
#        With these samples use the sampling rules to form a final
#        set. 
#        
#        For example, sampling using rvs create the following:
#        
#        var1, var2, var3
#        
#        The sampling rules then enforce the following:
#        
#        param_a = var1 * var2
#        param_b = var3 / var2
        
#    """        
#    
#    samples=sample_norules(rvs,input_params)
    
    

