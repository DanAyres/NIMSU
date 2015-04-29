'''
Created on 27 Apr 2015

@author: daiel
'''



from FORTRAN_Modules import orthpol
import numpy as np

class OrthPoly():

    def __init__(self, n_coeff, solver, n_unions, domain, leftinf, rightinf, dist, dist_args, scale):
        
        
        self.dist=dist.lower()
        if  self.dist not in ['uniform', 'normal', 'beta']:
            print "Dist not recognised"
            raise
        
        if self.dist =='uniform':
            dist=1
            self.PDF=0.5
        elif self.dist=='normal':
            self.PDF=1.0/np.sqrt(np.pi)
            dist=2
        elif self.dist=='beta':
            dist=3
        self.solver=solver.lower()
        if self.solver not in ['sti','lanc']:
            print "Solver not recognised"
            raise

        if self.solver=='sti':
            solver = 1
        else:
            solver = 2
                    
        self.n_poly=n_coeff-1
        self.domain=domain
        self.leftinf=leftinf
        self.rightinf=rightinf
        
        alpha, beta=orthpol.orthpol_interface.recursion_coefficients(n_coeff, 
                                                                     solver, 
                                                                     n_unions, 
                                                                     self.domain, 
                                                                     self.leftinf, 
                                                                     self.rightinf, 
                                                                     dist, 
                                                                     dist_args, 
                                                                     scale)


        eps = 1.0E-10
        self.abcissa, self.weights = orthpol.orthpol_interface.quadrature_rule(alpha,beta,eps)
        
        self.coefficients=orthpol.orthpol_interface.polynomial_coefficients(alpha,beta)

        self.norms=orthpol.orthpol_interface.normalisation_factors(self.coefficients,self.abcissa,self.weights)
        
        self.coefficients = self.coefficients.T

