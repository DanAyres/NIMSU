'''
Created on 14 Apr 2015

@author: daiel
'''


class Data():
    
    def __init__(self,val,listt=[],hdf5='none'):
        self.val=val
        self.list=listt
        self.hdf5=hdf5
    
    def __add__(self,other):
        pass
    
    def __sub__(self,other):
        pass

    def __mul__(self,other):
        pass
    
    def __div__(self,other):
        pass
    
    def __radd__(self,other):
        pass    
    
    
class singleData(Data):
    
    def __add__(self,other):
        try:
            return singleData( self.val + other.val)
        except AttributeError:
            return singleData( self.val + other)
    
    def __radd__(self,other):
        return singleData( self.val + other)
    
    
    def __sub__(self,other):
        try:
            return singleData( self.val - other.val)
        except AttributeError:
            return singleData( self.val - other)
        
    def __mul__(self,other):
        try:
            return singleData( self.val * other.val)
        except AttributeError:
            return singleData( self.val * other)
        
    def __rmul__(self,other):
        return singleData( self.val * other)
    
    def __div__(self,other):
        try:
            return singleData( self.val / other.val)
        except AttributeError:
            return singleData( self.val / other)
        
    def __pow__(self,other):
        return singleData(self.val**other)
    
class listData(Data):
    
    def __add__(self,other):
        try:
            return listData( self.val + other.val, listt=(self.list + other.list) )
        except AttributeError:
            return listData( self.val + other, listt=(self.list + other) )
    
    def __radd__(self,other):
        return listData( self.val + other, listt=(self.list + other) )
    
    
    def __sub__(self,other):
        try:
            return listData( self.val - other.val, listt=(self.list - other.list) )
        except AttributeError:
            return listData( self.val - other, listt=(self.list - other) )
        
    def __mul__(self,other):
        try:
            return listData( self.val * other.val, listt=(self.list * other.list) )
        except AttributeError:
            return listData( self.val * other, listt=(self.list * other) )
        
    def __rmul__(self,other):
        return listData( self.val * other, listt=(self.list * other) )
    
    def __div__(self,other):
        
        try:
            return listData( self.val / other.val, listt=(self.list / other.list) )
        except AttributeError:
            return listData( self.val / other, listt=(self.list / other) )     
        
    def __pow__(self,other):
        return singleData(self.val**other,listt=(self.list**other))   
    
class hdf5Data(Data):
    
    def __add__(self):
        pass 