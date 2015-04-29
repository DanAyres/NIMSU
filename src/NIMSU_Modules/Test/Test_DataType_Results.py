'''
Created on 14 Apr 2015

@author: daiel
'''
import unittest
from NIMSU_Modules.DataType_Results import singleData, listData, hdf5Data
import numpy as np

class Test(unittest.TestCase):

    def testsingleData(self):
        
        #
        # ----------- Addition methods -----------
        #
        a=singleData(4)
        b=singleData(4)
        
        # Basic add
        c=a+b
        self.assertEqual(c.val, 8)
        
        # check the exception handling
        a+=3
        self.assertEqual(a.val, 7)
        
        # check reverse add
        c=sum([a,b])
        self.assertEqual(c.val, 11)
        
        #
        # ----------- Subtraction methods -----------
        #
        a=singleData(4)
        b=singleData(4)
        
        c=a-b
        self.assertEqual(c.val, 0)
        
        a-=3
        self.assertEqual(a.val, 1)
        
        #
        # ----------- Multiplication methods -----------
        #
        a=singleData(4)
        b=singleData(4)
        
        c=a*b
        self.assertEqual(c.val, 16)
        
        
        a*=3
        self.assertEqual(a.val, 12)
        
        
        c=3*b
        self.assertEqual(c.val, 12)
        
        #
        # ----------- Division methods -----------
        #
        a=singleData(4)
        b=singleData(3)
        
        c=a/b
        self.assertEqual(c.val, 1)
        
        c=b/3.0
        self.assertEqual(c.val, 1.0)
        
        #
        # ----------- Exponent methods -----------
        #
        a=singleData(4)
        b=singleData(3)
        
        c = a**2
        self.assertEqual(c.val, 16)
        c=b**3
        self.assertEqual(c.val, 27)
 
    def testlistData(self):
         
        #
        # ----------- Addition methods -----------
        #
        a=listData(4, listt=np.ones(4,float))
        b=listData(4, listt=3.5*np.ones(4,float))
         
        # Basic add
        c=a+b
        self.assertEqual(c.val, 8)
        self.assertTrue( np.allclose( c.list, 4.5*np.ones(4,float), rtol=1e-10, atol=1e-10))
         
        # check the exception handling
        a+=3
        self.assertEqual(a.val, 7)
        self.assertTrue( np.allclose( a.list, 4.0*np.ones(4,float), rtol=1e-10, atol=1e-10))
         
        # check reverse add
        c=sum([a,b])
        self.assertEqual(c.val, 11)
        self.assertTrue( np.allclose( c.list, 7.5*np.ones(4,float), rtol=1e-10, atol=1e-10))
         
        #
        # ----------- Subtraction methods -----------
        #
        a=listData(4, np.ones(4,float))
        b=listData(4, 3.5*np.ones(4,float))
          
        c=a-b
        self.assertEqual(c.val, 0)
        self.assertTrue( np.allclose( c.list, -2.5*np.ones(4,float), rtol=1e-10, atol=1e-10))
          
        a-=3
        self.assertEqual(a.val, 1)
        self.assertTrue( np.allclose( a.list, -2.0*np.ones(4,float), rtol=1e-10, atol=1e-10))
          
        #
        # ----------- Multiplication methods -----------
        #
        a=listData(4, np.ones(4,float))
        b=listData(4, 3.5*np.ones(4,float))
          
        c=a*b
        self.assertEqual(c.val, 16)
        self.assertTrue( np.allclose( c.list, 3.5*np.ones(4,float), rtol=1e-10, atol=1e-10))
          
          
        a*=3
        self.assertEqual(a.val, 12)
        self.assertTrue( np.allclose( a.list, 3.0*np.ones(4,float), rtol=1e-10, atol=1e-10))
          
          
        c=3*b
        self.assertEqual(c.val, 12)
        self.assertTrue( np.allclose( c.list, 10.5*np.ones(4,float), rtol=1e-10, atol=1e-10))
          
        #
        # ----------- Division methods -----------
        #
        a=listData(4, np.ones(4,float))
        b=listData(3, 2.0*np.ones(4,float))
          
        c=a/b
        self.assertEqual(c.val, 1)
        self.assertTrue( np.allclose( c.list, 0.5*np.ones(4,float), rtol=1e-10, atol=1e-10))
          
        c=b/3.0
        self.assertEqual(c.val, 1.0)
        self.assertTrue( np.allclose( c.list, (2.0/3.0)*np.ones(4,float), rtol=1e-10, atol=1e-10))
  
  
        #
        # ----------- Power methods -----------
        #
        a=listData(4, np.ones(4,float))
        b=listData(3, 2.0*np.ones(4,float))
          
        c=a**2
        self.assertEqual(c.val, 16)
        self.assertTrue( np.allclose( c.list, np.ones(4,float), rtol=1e-10, atol=1e-10))
          
        c=b**3
        self.assertEqual(c.val, 27)
        self.assertTrue( np.allclose( c.list, 8*np.ones(4,float), rtol=1e-10, atol=1e-10))

#     def testhdf5Data(self):
#         
#         a=hdf5Data(0.0,hdf5='test')
#         
#         print a.val
#         print a.list
#         print a.hdf5
#         
#         pass 
