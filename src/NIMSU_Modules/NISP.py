'''
Created on 24 Apr 2015

@author: daiel
'''

import numpy as np
from NIMSU_Modules.SparseGridBaseClass import sparse_grid
from NIMSU_Modules.DataTypePCE import PCE
from operator import mul


def Project2Basis(ele, PDF):
    for poly in ele.pce:
        for pos,p in enumerate(ele.SparseGrid.Index):
            pts=ele.SparseGrid.Points[poly[2],p]
            poly[0] += PDF * ele.SparseGrid.Weights[pos] * ele.PolyProd(poly[1], pts, poly[2]) * ele.SparseGrid.Values[p]
        poly[0] /= ele.PolyNorms(poly[1], poly[2])
        
        #print  poly[0].val, poly[0].list
    

def NISP(control_info, sampler, engine):
    
    newPCE=PCE(control_info, sampler, engine.Nominal)                 
        
    # Integrate over each elements    
    adapt_Int_Complete=False
    while not adapt_Int_Complete:
        
        Samples = newPCE.CalculateSamples()
            
        Results=engine.interface(Samples)
            
        adapt_Int_Complete = newPCE.UpdateIntegrals(Results)
            
    # Project onto the polynomial basis   
    for ele in newPCE.Elements:
        Project2Basis(ele, newPCE.PDF)
        
    #newPCE.CalculateErrorMeasures()
            