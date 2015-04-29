'''
Created on 6 Feb 2015

@author: Some one else
'''

import os
import subprocess
import time

class ProcessQueue:
    """
        Manages a queue of jobs which are spawned as sub-processes
        in controlled numbers as defined by pool_size
        
        ARGS:
            pool_size: number of jobs in the pool
            
        Returns:
        
        Raises:
    
    """
    
    def __init__(self, pool_size):
        self.queue=[]
        self.pool_size=pool_size
        
    def append(self, job, wkdir):
        self.queue.append((job,wkdir))
        
    def execute(self):
        """
            Run the jobs in the queue.
            Block until they are all finished.
        
        """
        
        
        open_jobs=set()
        
        
        # While a job remains in the queue, extract one from the queue.
        while len(self.queue) > 0:
        
            job, wd = self.queue.pop()
        
            while len(open_jobs) == self.pool_size:
                poll_jobs(open_jobs)
                
            cwd=os.getcwd()
            os.chdir(wd)
            try:
                open_jobs.add(subprocess.Popen(job,shell=True))
            finally:
                os.chdir(cwd)
            
            while len(open_jobs) > 0:
                poll_jobs(open_jobs)
                
                
def poll_jobs(open_jobs):
    
    time.sleep(1)
    done_jobs = set([job for job in open_jobs if job.poll() != None])
    open_jobs -= done_jobs