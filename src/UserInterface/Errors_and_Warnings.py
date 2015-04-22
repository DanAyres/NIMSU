'''
Created on 9 Feb 2015

@author: daniel
'''


def ReadError(out, method_name, line, msg):
    out.write("#############################################################################\n")
    out.write("                           READ ERROR\n\n")
    out.write("   In: %s \n\n"%method_name)
    out.write("   At line: %s \n\n"%line)
    out.write("   %s \n\n"%msg)
    out.write("#############################################################################\n")

def ProcessingError(out, method_name, msg):
    out.write("#############################################################################\n")
    out.write("                           PROCESSING ERROR\n\n")
    out.write("   In: %s \n\n"%method_name)
    out.write("   %s \n\n"%msg)
    out.write("#############################################################################\n")
    
              
def WarningMessage(out, method_name, msg):
    out.write("#############################################################################\n")
    out.write("                           Warning\n\n")
    out.write("   In: %s \n\n"%method_name)
    out.write("   %s \n\n"%msg)
    out.write("#############################################################################\n")