'''
Created on 29 Apr 2015

@author: daiel
'''

import unittest
import os

suite=unittest.TestSuite()

# Create a list of all of the tests for the HDMR_Modules
tests=[]
for r,d,files in os.walk('NIMSU_Modules/Test'):
    for f in files:
        if f[:4] == 'Test' and f[-3:]=='.py':
            tests.append('NIMSU_Modules.Test.' + f[:-3])
            
# Add the tests to the test suite
for t in tests:
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))

# Run all the tests
unittest.TextTestRunner(verbosity=2,buffer=True).run(suite)