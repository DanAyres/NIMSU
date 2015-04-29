'''
Created on 3 Feb 2015

@author: daniel
'''

import os
import numpy as np
from NIMSU_Modules.UIErrorsandWarnings import ProcessingError, ReadError
from NIMSU_Modules.SampTransforms import FormCorrelation

class param:
    """
        Defines a set of data parameters. Can be mutually correlated / uncorrelated.
        All parameters in the class are described using a common PDF.
        
        ARGS:
            NParams:  Number of parameters
            Names:    Name tag for each parameter
            Means:    Mean value of each parameter
            Sigma:    Standard deviations
            C:        Covariance matirx
            Depend:   Specify the relation between parameters (yes/no/corr)
            dist:     Specify the distribution (normal/uniform/etc..)
            decomp:   The method used to decompose the covariance matrix
            kernel:   Correlation kernel
            Cov_Len:  Correlation length of kernel.
            kl_theta: Variance ratio represented by the KL expansion (default=1.0)
            
        Returns:
        
        Raises:
    """
    def __init__(self, NParams, Names, Means, Sigma, C, Depend, dist, decomp, kernel, Cov_Len, kl_theta=1.0):
        
        self.NParams=NParams     
        self.Names=Names         
        self.Means=Means         
        self.Sigma=Sigma         
        self.Corr_Matrix= C                
        self.Depend=Depend       
        self.Dist=dist           
        self.Decomp=decomp       
        
        # Optional parameters        
        self.Kernel=kernel
        self.Kernel_args=Cov_Len     
        self.Kl_theta=kl_theta  

        # May include -supported distributions
        #             -supported covariance kernels
        
        
        
class Control:
    """
        Contains the control information for the program.
    """

    def __init__(self):
        self.run_ref='test_run'
        self.distribution='normal'
        self.param_infile='test.paramfile'
        self.template_file='test.template'
        self.input_file='test.infile'
        self.output_file='test.outfile'
        self.run_script='run_test.sh'
        self.setup_script=None
        self.setup_files=None
        self.num_procs=1
        self.scale_sampling=False
#       MC
        self.use_mc=True
        self.mc_method='latinhyper'
        self.mc_runs=100
        self.mc_pdf=False
        self.mc_pdf_bins=100
        self.mc_dump=100
        self.mc_seed=121
#       SCM
        self.use_scm=False
        self.scm_quad_adapt=False
        self.scm_quad_tol=1E-3
        self.scm_sg_level=[2]
        self.scm_pce=True
        self.scm_specify_order=False
        self.scm_poly_order=3
        self.scm_pdf=False
        self.scm_pdf_bins=100
        self.scm_pdf_samp=10000
#       HDMR
        self.use_hdmr=False
        self.hdmr_sg_level=[2,2]
        self.hdmr_pce=True
        self.hdmr_poly_order=3
        self.hdmr_pdf=False
        self.hdmr_adapt=True
        self.hdmr_quad_adapt=False
        self.hdmr_quad_tol=1E-3
        self.hdmr_metric=1E-1
        self.hdmr_pdf_bins=100
        self.hdmr_pdf_samp=10000
#       Sensitivity
        self.use_sen=False
        self.sen_order=1
        self.sen_pce=False
        
                
def ReadControlFile(ctrl_file, out):
    """
        Read the control file.
        
        ARGS:
            ctrl_file: open control file
            out      : open file or IOString for warnings and errors
            
        Returns:
            Control() object
            
        Raises:
            ValueError
    """
    
    control_info=Control()
    
    for line in ctrl_file:
        data=line.partition('#')
        
        if data[0].strip() != ('' and '\n'):
            
            if not data[0].strip().find(':'):
                msg="Delimiter ':' missing from parameter assignment.\n\n   %s"%data[0].strip()
                ProcessingError(out, "UserInterface.interface.ReadControlFile", msg)
                raise ValueError
            
            keyword,value=data[0].strip().split(':')
            
            
            name = keyword.strip().replace(' ','_').lower()
            value=value.strip()
            
            try:
                getattr(control_info,name)
            except AttributeError:
                msg='Incorrect Keyword "%s" in Control File'%keyword
                ProcessingError(out, "UserInterface.interface.ReadControlFile", msg)
                raise ValueError
            else:
                if value != '':
                    if name in ['num_procs',
                                'mc_runs' ,
                                'mc_seed' ,
                                'mc_pdf_bins', 
                                'scm_pdf_bins', 
                                'hdmr_pdf_bins',
                                'scm_pdf_samp', 
                                'hdmr_pdf_samp',
                                'hdmr_poly_order',
                                'scm_poly_order',
                                'mc_dump']:
                        value=int(value)
                    if name in ['sen_order']:
                        value=int(value)
                        if value not in [1,2,4]:
                            print 'Sensitivity accuracy must be 1,2 or 4. Not:', value
                            print 'using first order accurate'
                            value = 1
                    if name in ['use_mc',
                                'use_scm',
                                'use_hdmr',
                                'use_sen',
                                'mc_pdf',
                                'scm_pdf',
                                'hdmr_pdf',
                                'hdmr_adapt',
                                'hdmr_quad_adapt',
                                'scm_quad_adapt',
                                'scm_pce',
                                'scm_specify_order',
                                'hdmr_pce',
                                'scale_sampling',
                                'sen_pce']:
                        if value.lower() in['yes','true']:
                            value = True
                        elif value.lower() in['no','false']:
                            value = False
                        else:
                            raise ValueError('Expecting yes/no/true/false. Not', value)
                    if name in ['scm_sg_level',
                                'hdmr_sg_level']:
                        values=value.split(',')
                        value=[int(v) for v in values]
                    if name in ['hdmr_metric',
                                'scm_quad_tol',
                                'hdmr_quad_tol']:
                        value = float(value)
                    if name in ['setup_files']:
                        values = value.split(',')
                        value = [v.strip() for v in values]
                        
                    setattr(control_info,name,value)
                
                else:
                    msg= "No value given for '%s' keyword"%keyword
                    ProcessingError(out, "UserInterface.interface.ReadControlFile", msg)
                    raise ValueError
                    
                    
    return control_info


def ReadParamFile(param_file, out):
    
    parameter_list=[]
       
    for line in param_file:
        data=line.partition('#')[0].split()
        if data:
            
            #######  FILE  #######
            #Read uncertainty parameters from a separate file. 
            if data[0].lower() == "file":
                if len(data) < 5:
                    msg="Insufficient values. Expecting at least 5"
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                
                # Name of external subfile
                fname=data[1]
                
                # Number of parameters
                try:
                    NParams=int(data[2])
                except ValueError:
                    msg="The third argument '%s' should be of type integer."%data[2]
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                
                # Probability distribution
                dist=data[3].lower()
                supported_dists=["normal", "uniform", "lognormal"]
                if dist not in supported_dists:
                    msg="Distribution of type '%s' not supported.\n   Supported distributions are: %s"%(dist,supported_dists)
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                
                # Type of correlation dependency
                depend=data[4].lower()
                if depend not in ["yes", "no", "corr"]:
                    msg="Dependency of type '%s' not supported.\n   Dependency should be one of: %s"%(depend,["yes", "no", "corr"])
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                
                # Read covariance data or generate using specified kernel
                if depend == "corr":
                    
                    if len(data) < 7:
                        msg="Expecting further arguments : <Decomp> and <CovMat>"
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                    
                    # Method for decomposing the covariance matrix
                    decomp=data[5].lower()
                    if decomp not in ["kl","chol"]:
                        msg="Allowed decompositions are 'kl' or 'chol'. Not '%s'"%decomp
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError

                    # Is the covariance matrix provided in the external file                    
                    covmat=data[6].lower()
                    if covmat not in ["yes","true","no","false"]:
                        msg="Allowed values of covmat are yes/true/no/false. Not '%s'"%covmat
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                    
                    if covmat in ["yes","true"]:
                        ReadCov=True
                    else:
                        # Read additional parameters describing covariance kernel
                        ReadCov=False
                        
                        if len(data)<9:
                            msg="Expecting further arguments : kernel and kernel_params"
                            ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                            raise ValueError
                        
                        # Name of covariance kernel
                        kernel=data[7].lower()
                        cov_kernels=['exponential']
                        if kernel not in cov_kernels:
                            msg="Covariance kernel '%s' not defined."%kernel
                            ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                            raise ValueError
                        
                        # Correlation length
                        try:
                            kernel_args=float(data[8])
                        except ValueError:
                            msg="The value of kernel_args: '%s' should be of type float"%kernel_args
                            ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                            raise ValueError
                        
                    #Start reading sub-file
                    if not os.path.isfile(fname):
                        msg="File '%s' does not exist."%fname
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                 
                # depend=yes/no
                else:
                    decomp='kl'
                    kernel='exponential'
                    kernel_args=1.0
                    
                Names,Means,Sigma,C=ReadSubFileData(open(fname,'r'), ReadCov, NParams)
                    
                    
                parameter_list.append(param(NParams, Names, Means, Sigma, C, depend, dist, decomp, kernel, kernel_args) )
                    
                    
            #######  DATA  #######
            elif data[0].lower() == "data":
                
                if len(data) < 7:
                    msg="Insufficient values. Expecting at least 7"
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                
                # Data tag
                name=data[1]
                
                # Number of parameters
                try:
                    NParams=int(data[2])
                except ValueError:
                    msg="The third argument '%s' should be of type integer."%data[2]
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                names = [name+str(i+1) for i in range(NParams)]
                
                # Mean values
                if  not ( not data[3].find('[')<0 and not data[3].find(']')<0):
                    msg="Mean values not enclosed with brackets '[ ]'. \n-- Mean values read: %s\n-- Remove all white space."%data[3]
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                means=np.asarray([float(x) for x in data[3].split('[')[-1].split(']')[0].split(',') ])
                if means.shape[0] !=NParams:
                        msg="Number of mean values must be equal to %s. Not %s"%(NParams,means.shape[0])
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                
                # Standard deviations
                scale=data[4].partition('%')
                if(scale[1]):
                    Sigma=np.asarray([float(scale[0])/100.0 for i in range(NParams)])
                    if Sigma.shape[0] !=NParams:
                        msg="Number of SDs values must be equal to %s. Not %s"%(NParams,Sigma.shape[0])
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                    Sigma = Sigma.__mul__(means)
                else:
                    if  not ( not data[4].find('[')<0 and not data[4].find(']')<0):
                        msg="SD values not enclosed with brackets '[ ]'\n-- or white space between RSD and %% symbol. \n-- SD values read: %s\n-- Remove all white space."%data[4]
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                    Sigma=np.asarray([float(x) for x in scale[0].split('[')[-1].split(']')[0].split(',') ])
                   
                # Type of distribution 
                dist=data[5].lower()
                supported_dists=["normal", "uniform", "lognormal"]
                if dist not in supported_dists:
                    msg="Distribution of type '%s' not supported."%dist
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError                       
                
                # Type of correlation dependency
                depend=data[6].lower()
                if depend not in ["yes", "no", "corr"]:
                    msg="Dependency of type '%s' not supported"%depend
                    ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                    raise ValueError
                
                C=None
                if depend=='corr':
                    
                    if len(data) < 10:
                        msg="Expecting further arguments : <Decomp> , <kernel> and <kernel_args>"
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                    
                    # Method for decomposing the covariance matrix
                    decomp=data[7].lower()
                    if decomp not in ["kl","chol"]:
                        msg="Allowed decompositions are 'kl' or 'chol'. Not '%s'"%decomp
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError

                    # Name of covariance kernel
                    kernel=data[8].lower()
                    cov_kernels=['exponential']
                    if kernel not in cov_kernels:
                        msg="Covariance kernel '%s' not defined\nAvailable kernels are: %s."%(kernel,cov_kernels)
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError
                        
                    # Correlation length
                    try:
                        kernel_args=float(data[9])
                    except ValueError:
                        msg="The value of kernel_args: '%s' should be of type float"%kernel_args
                        ReadError(out, "UserInterface.interface.ReadParamFile", line, msg)
                        raise ValueError

                # depend=yes/no
                else:
                    decomp='kl'
                    kernel='exponential'
                    kernel_args=1.0
                    
                parameter_list.append(param(NParams, names, means, Sigma, C, depend, dist, decomp, kernel, kernel_args) )

    return parameter_list




                
def ReadSubFileData(datafile, ReadCov, NParams, out):
    """
        Read data from a pre-formatted file.
        The names, means and standard deviations of each parameter 
        are given in the following format
        
        name1  mean1 sd1
        or
        name2  mean2 %rsd
        
        where rsd is the relative standard deviation given in percent.
        
        If ReadCov is true the correlation matrix is read immediately 
        after the names and means etc. The full (NParams,NParams) 
        correlation matrix is expected. 
        
        It is also expected that the diagonals of the correlation matrix
        are equal to unity with an absolute tolerance of 1E-10. If this
        is not the case it is assumed that the user has entered the
        COVARIANCE matrix which will be normalised. This induces a 
        consistency check with the standard deviations read previously.
        If these are different, the SDs determined from the covariance
        matrix are used. 
        
        Args:
            datafile: open file object to read
            ReadCov: boolean, toggle reading of CORRELATION matrix.
            NParams: Number of parameters in the data set/CORRELATION matrix
            out    : where to direct output (warnings/errors etc.)
        Returns:
            names: name corresponding to each parameter.
            means: the mean value of each parameter.
            std  : the standard deviation of each parameter.
            C    : The CORRELATION matrix if ReadCov=True. None if ReadCov=False.
        Raises:
            ValueError    : For dubious input
    """
    
    counter=1
    names=[]
    means=np.zeros(NParams,float)
    std=np.zeros(NParams,float)
    C=None
    
    while(True):
        
        if counter>NParams:
            break
        
        line=datafile.readline()            
        if not line: break    
        split=line.partition('#')[0].split()  # Remove comments
        
        if split:
                   
            if len(split)!=3:
                msg="--Expecting 3 values: name mean RSD"
                ReadError(out, "UserInterface.interface.ReadSubFileData", line, msg)
                raise ValueError
                
            names.append(split[0])
            means[counter-1]=float(split[1])
            
            scale=split[2].partition('%')
            if(scale[1]):
                std[counter-1] = means[counter-1] * float(scale[0])/100
            else:
                std[counter-1] = float(scale[0])
            
            counter+=1
    
    
    if ReadCov:
        counter=1
        C=np.zeros([NParams,NParams],float)
        while(True):
        
            if counter>NParams:
                break
            
            line=datafile.readline()    
            if not line: break            
            split=line.partition('#')[0].split()  # Remove comments
        
            if split:
                
                if len(split)!=NParams:
                    msg="Expecting %s values in row %s of correlation matrix. Not %s"%(NParams,counter,len(split))
                    ReadError(out, "UserInterface.interface.ReadSubFileData", line, msg)
                    raise ValueError
                
                C[counter-1,:] = np.asarray([float(val) for val in split])
                
                counter+=1
                
        if(counter-1 < NParams):
            print "WARNING: Only %s out of %s rows of the correlation matrix read from file"%((counter-1),NParams)
    
        # Check diagonals of correlation matrix
        if not np.allclose(C.diagonal(), np.ones(NParams,float), rtol=1e-05, atol=1e-08):
            std=np.sqrt(C.diagonal())
            C=FormCorrelation(C, std, NParams)
            
        if ( np.amin(C) < -1.0 ) or ( np.amax(C) > 1.0 ) :
            msg="Correlation matrix is not bounded!!\n   Max: %12.5E  Min: %12.5E"%(np.amax(C),np.amin(C))
            ProcessingError(out, "UserInterface.interface.ReadSubFileData", msg)
            raise ValueError
            
            
    
    return names,means,std,C

