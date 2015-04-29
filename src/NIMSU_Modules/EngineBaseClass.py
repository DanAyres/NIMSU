'''
Created on 6 Feb 2015

@author: daniel
'''

import os
from string import Template
from process_queue.multi_process import ProcessQueue

class EngineInterfaceBase:
    """
        Base class for generic engine interface functions
    
    """
    
    def __init__(self, control_info):
        self.control_info=control_info
        
    
    def open_read_file(self,filename):
        return open(filename,'r')
    
    def open_write_file(self,filename):
        return open(filename,'w')
    
    def write_input_file(self):
        
    def results_template(self):
        
    def read_output_files(self):
        
    def interface(self, samples):
        """
            Main engine interface function
            
            For each engine run required:
                -Create a separate subdirectory for run files
                -Write the engine input files using write_input_file
                -Spawn the engine run using ProcessQueue
                -Wait for all runs to complete
                -Read the required results using read_output_files()
            
            ARGS:
                samples: a list of sampled parameters
                
            Return:
                results: a list of results which maps directly to the input samples
        
        """
        
        NSamples = len(samples)
        
        run_dirs=[]
        for n in range(NSamples):
            
            # create unique run-directory identifier.
            run_DIR = 'something'
            run_dirs.append(run_DIR)
            
            
            # create run directory and necessary files to execute code
            # Input file, run script, etc.
            
            
            # Generate the process queue
            queue=ProcessQueue(NProcs)
            for d in run_dirs:
                queue.append(run_script,d)
            queue.execute()
            
            
        # Extract results
        results=[] 
        

        
        return results