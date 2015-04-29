'''
Created on 4 Feb 2015

@author: daniel
'''
import unittest
from UserInterface.UIBaseClass import ReadControlFile, ReadParamFile, ReadSubFileData
from StringIO import StringIO
import numpy as np
import sys

class TestUserInterface(unittest.TestCase):
    
    
    def testReadControlFile_General(self):
        """
            Check that the general control parameters
            have been properly assigned to the
            control_info class.
        """
        out=sys.stdout
        instring = StringIO("Run_Ref : test_run\n" +
                           "Param Infile : test.paramfile\n " +
                           "template File :test.template\n" +
                           "input file :test.infile #Trailing comment\n" +
                           "output file :test.outfile\n" +
                           "run script :run_test.sh\n" +
                           "setup script : setup.sh\n" + 
                           "setup files : mycode.setup\n" +
                           "num procs : 3\n" +
                           "scale sampling : False\n")

        run_info = ReadControlFile(instring,out)
        
        
        self.assertEqual(run_info.run_ref, "test_run")
        self.assertEqual(run_info.param_infile, "test.paramfile")
        self.assertEqual(run_info.template_file, "test.template")
        self.assertEqual(run_info.input_file, "test.infile")
        self.assertEqual(run_info.output_file, "test.outfile")
        self.assertEqual(run_info.run_script, "run_test.sh")
        
        self.assertEqual(run_info.setup_script, "setup.sh")
        self.assertEqual(run_info.setup_files, ["mycode.setup"])
        self.assertEqual(run_info.run_script, "run_test.sh")
        self.assertEqual(run_info.run_script, "run_test.sh")
        
        
    def testReadControlFile_MonteCarlo(self):
        """
            Check that the Monte Carlo control parameters
            have been properly assigned to the
            control_info class.
        """
        out=sys.stdout
        instring = StringIO( "use_mc : True\n"+
                             "mc method : latinhyper\n"+
                             "mc runs : 100\n"+
                             "mc pdf : False\n"+
                             "mc pdf bins : 100\n"+
                             "mc dump :100\n"+
                             "mc seed :121")

        run_info = ReadControlFile(instring,out)
        
        
        self.assertEqual(run_info.use_mc, True)
        self.assertEqual(run_info.mc_method, "latinhyper")
        self.assertEqual(run_info.mc_runs, 100)
        self.assertEqual(run_info.mc_pdf, False)
        self.assertEqual(run_info.mc_pdf_bins,100)
        self.assertEqual(run_info.mc_dump, 100)
        self.assertEqual(run_info.mc_seed, 121)


    def testReadControlFile_SCM(self):
        """
            Check that the SCM control parameters
            have been properly assigned to the
            control_info class.
        """
        out=sys.stdout
        instring = StringIO("use scm: False\n"+
                            "scm quad adapt: False\n"+
                            "scm quad tol: 1E-3\n"+
                            "scm sg level: 2\n"+
                            "scm pce: True\n"+
                            "scm specify order: False\n"+
                            "scm poly order: 3\n"+
                            "scm pdf : False\n"+
                            "scm pdf bins: 100\n"+
                            "scm pdf samp: 10000")

        run_info = ReadControlFile(instring,out)
        
        
        self.assertEqual(run_info.use_scm, False)
        self.assertEqual(run_info.scm_quad_adapt, False)
        self.assertEqual(run_info.scm_quad_tol, 1E-3)
        self.assertEqual(run_info.scm_sg_level, [2])
        self.assertEqual(run_info.scm_pce,True)
        self.assertEqual(run_info.scm_specify_order, False)
        self.assertEqual(run_info.scm_pdf, False)
        self.assertEqual(run_info.scm_pdf_bins, 100)
        self.assertEqual(run_info.scm_pdf_samp, 10000)        
        
        
    def testReadControlFile_HDMR(self):
        """
            Check that the HDMR control parameters
            have been properly assigned to the
            control_info class.
        """
        out=sys.stdout
        instring = StringIO("use_hdmr:       False\n"+
                            "hdmr sg level  : 2, 2\n"+
                            "hdmr pce       : True\n"+
                            "hdmr poly order: 3\n"+
                            "hdmr pdf       : False\n"+
                            "hdmr adapt     : True\n"+
                            "hdmr quad adapt: False\n"+
                            "hdmr quad tol :  1E-3\n"+
                            "hdmr metric   :  1E-1\n"+
                            "hdmr pdf bins :  100\n"+
                            "hdmr pdf samp :  10000")

        run_info = ReadControlFile(instring,out)
        
        
        self.assertEqual(run_info.use_hdmr, False)
        self.assertEqual(run_info.hdmr_sg_level, [2,2])
        self.assertEqual(run_info.hdmr_pce, True)
        self.assertEqual(run_info.hdmr_poly_order, 3)
        self.assertEqual(run_info.hdmr_pdf,False)
        self.assertEqual(run_info.hdmr_adapt, True)
        self.assertEqual(run_info.hdmr_quad_adapt, False)
        self.assertEqual(run_info.hdmr_quad_tol, 1E-3)
        self.assertEqual(run_info.hdmr_metric, 1E-1)
        self.assertEqual(run_info.hdmr_pdf_bins, 100)
        self.assertEqual(run_info.hdmr_pdf_samp, 10000)
        
        
    def testReadControlFile_Misspelt(self):
        
        out=sys.stdout
        instring = StringIO("""Parameters : test.paramfile""")
        
        self.assertRaises(ValueError, ReadControlFile, instring,out)
        
        
    def testReadControlFile_MissingDelimiter(self):
        
        out=sys.stdout
        instring = StringIO("Param infile : test.infile\n" + 
                            "MC Runs  100 ")
        
        self.assertRaises(ValueError, ReadControlFile, instring, out)
        
        
    def testReadControlFile_MissingValue(self):
        
        out=sys.stdout
        instring = StringIO("Param infile : test.infile\n" + 
                            "MC Runs : ")
        
        self.assertRaises(ValueError, ReadControlFile, instring, out)    
        

    def testReadParamFile_FILE(self):
        """
            Check that ValueErrors are raised when 
            incorrect values are input.
        """
        
        out=sys.stdout
        # Requires 5 values
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat 9 Normal")
        self.assertRaises(ValueError, ReadParamFile, instring, out)
        
        # 'nine' should be of type int
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat nine Normal corr")
        self.assertRaises(ValueError, ReadParamFile, instring, out)

        # 'dave' is not a supported distribution        
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat 9 dave corr")
        self.assertRaises(ValueError, ReadParamFile, instring, out)

        #         
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat 9 uniform uncorr")
        self.assertRaises(ValueError, ReadParamFile, instring, out)
        
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat 9 uniform corr")
        
        self.assertRaises(ValueError, ReadParamFile, instring, out)
        
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat 9 uniform corr KL dave")
        
        self.assertRaises(ValueError, ReadParamFile, instring, out)
        
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "FILE fname.dat 9 uniform corr KL True")
        
        self.assertRaises(ValueError, ReadParamFile, instring, out)
        
        
    def testReadParamFile_DATA(self):
        out=sys.stdout
        
        # Too few values
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "DATA sig_f 3 [1.0,2.0,3.0] 1% normal")
        
        self.assertRaises(ValueError, ReadParamFile, instring, out)
        
        
        #Check param() object contains the correct values
        instring = StringIO("#Comment\n"+
                            "#Comment\n"+
                            "DATA sig_f 3 [1.0,2.0,3.0] 1% normal corr kl exponential 1.0")
        
        plist=ReadParamFile(instring, out)
        
        self.assertEqual(plist[0].Corr_Matrix, None, "Wrong covariance value")
        self.assertEqual(plist[0].Kernel_args, 1.0, "Wrong correlation length")
        self.assertEqual(plist[0].Depend, 'corr', "Wrong dependency")
        self.assertEqual(plist[0].NParams, 3, "Wrong number of parameters")
        self.assertEqual(plist[0].Names, ['sig_f1', 'sig_f2', 'sig_f3'], "Wrong names")
        self.assertEqual(plist[0].Decomp, 'kl', "Wrong decomposition")
        self.assertEqual(plist[0].Dist, 'normal', "Wrong distribution")
        self.assertEqual(plist[0].Kernel, 'exponential', "Wrong kernel")
        self.assertEqual(plist[0].Kl_theta, 1.0, "Wrong kernel argument")
        
        
        self.assertTrue(np.allclose(plist[0].Means, np.asarray([1.0,2.0,3.0]), rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(plist[0].Sigma, np.asarray([0.01,0.02,0.03]), rtol=1e-05, atol=1e-08))
        
        
    def testReadSubFileData(self):
        
        out=sys.stdout
        # Test return values and absolute uncertainties
        instring = StringIO("#Comment line\n" +
                            "sig_f1  1.0  0.1\n" +
                            "sig_f2  2.0  0.2\n" +
                            "sig_f3  3.0  0.3")
        
        ReadCov = False
        NParams = 3
        
        names,means,std,C = ReadSubFileData(instring, ReadCov, NParams, out)
        
        self.assertEqual(names, ["sig_f1","sig_f2","sig_f3"], "Name list not equal")
        self.assertEqual(C, None, "Name list not equal")
        self.assertTrue(np.allclose(means, np.asarray([1.0,2.0,3.0]), rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(std, np.asarray([0.1,0.2,0.3]), rtol=1e-05, atol=1e-08))
        

        # Test return values and relative uncertainties in %
        instring = StringIO("#Comment line\n" +
                            "sig_f1  1.0  1%\n" +
                            "sig_f2  2.0  10%\n" +
                            "sig_f3  3.0  3%")
        
        ReadCov = False
        NParams = 3
        
        names,means,std,C = ReadSubFileData(instring, ReadCov, NParams, out)
        
        self.assertTrue(np.allclose(std, np.asarray([0.01,0.2,0.09]), rtol=1e-05, atol=1e-08))
        
        
        # Test missing value
        instring = StringIO("#Comment line\n" +
                            "sig_f1  1.0  0.1\n" +
                            "sig_f2  2.0  0.2\n" +
                            "sig_f3  3.0 ")
        
        self.assertRaises(ValueError, ReadSubFileData,instring, ReadCov, NParams, out)
        
    def testReadSubFileData_correlation(self):
        
        out=sys.stdout    
        
        # Test unbounded correlation matrix
        instring = StringIO("#Comment line\n" +
                        "sig_f1  1.0  0.1  #inline comment\n" +
                        "sig_f2  2.0  0.2\n" +
                        "sig_f3  3.0  0.3\n" + 
                        "#Comment line\n" +
                        "1.0 2.0 3.0 #inline comment\n" +
                        "4.0 1.0 6.0\n" +
                        "#comment line\n" +
                        "7.0 8.0 1.0\n" )
        
        ReadCov = True
        NParams = 3
        self.assertRaises(ValueError, ReadSubFileData,instring, ReadCov, NParams, out)
        
        # Test form correlation from covariance
        instring = StringIO("#Comment line\n" +
                        "sig_f1  1.0  1.0  #inline comment\n" +
                        "sig_f2  1.0  1.0\n" +
                        "sig_f3  1.0  1.0\n" + 
                        "#Comment line\n" +
                        "1.0 0.0 0.0 #inline comment\n" +
                        "0.0 4.0 0.0\n" +
                        "#comment line\n" +
                        "0.0 0.0 9.0\n" )    
        
        names,means,std,C = ReadSubFileData(instring, ReadCov, NParams, out)
        test_C=np.zeros([NParams,NParams],float)
        test_sd=np.zeros([NParams],float)
        test_sd[0]=1.0;test_sd[1]=2.0;test_sd[2]=3.0
        
        test_C[0,0]=1.0;test_C[1,0]=0.0;test_C[2,0]=0.0
        test_C[0,1]=0.0;test_C[1,1]=1.0;test_C[2,1]=0.0
        test_C[0,2]=0.0;test_C[1,2]=0.0;test_C[2,2]=1.0
        self.assertTrue(len(names),NParams)
        self.assertTrue(len(means),NParams)
        self.assertTrue(np.allclose(C, test_C, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(std, test_sd, rtol=1e-05, atol=1e-08))
        
        
        
        