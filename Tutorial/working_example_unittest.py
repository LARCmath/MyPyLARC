#!/usr/bin/env python3

 #*################################################################
 #                                                                #
 # Copyright (C) 2014, Institute for Defense Analyses             #
 # 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           #
 # This material may be reproduced by or for the US Government    #
 # pursuant to the copyright license under the clauses at DFARS   #
 # 252.227-7013 and 252.227-7014.                                 #
 #                                                                #
 # LARC : Linear Algebra via Recursive Compression                #
 # Authors:                                                       #
 #   - Steve Cuccaro (IDA-CCS)                                    #
 #   - John Daly (LPS)                                            #
 #   - John Gilbert (UCSB, IDA adjunct)                           #
 #   - Mark Pleszkoch (IDA-CCS)                                   #
 #   - Jenny Zito (IDA-CCS)                                       #
 #                                                                #
 # Additional contributors are listed in "LARCcontributors".      #
 #                                                                #
 # Questions: larc@super.org                                      #
 #                                                                #
 # All rights reserved.                                           #
 #                                                                #
 # Redistribution and use in source and binary forms, with or     #
 # without modification, are permitted provided that the          #
 # following conditions are met:                                  #
 #   - Redistribution of source code must retain the above        #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer.                                      #
 #   - Redistribution in binary form must reproduce the above     #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer in the documentation and/or other     #
 #     materials provided with the distribution.                  #
 #   - Neither the name of the copyright holder nor the names of  #
 #     its contributors may be used to endorse or promote         #
 #     products derived from this software without specific prior #
 #     written permission.                                        #
 #                                                                #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         #
 # CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    #
 # INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       #
 # MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       #
 # DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        #
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   #
 # NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   #
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       #
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
 # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, #
 # EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             #
 #                                                                #
 #*################################################################

from __future__ import print_function

import os
import sys 

if 2 == len(sys.argv):
    sys.path.append(sys.argv[1])
    sys.argv = [sys.argv[0]]

sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
# from pylarc import *
import MyPyLARC as mypy
import numpy as np
from ctypes import *
import unittest


##
# \file working_example_unittest.py
#
# \brief An example unit test. (also see test_unittest_matmath.py)
#
# 


# THIS IS AN EXAMPLE OF USING THE PYTHON UNIT TEST LIBRARY.
# DOCUMENTATION LIVES ON THE OUTSIDE WEB AT
# https://docs.python.org/2/library/unittest.html .

class TestLARCmath(unittest.TestCase):
    def testAdd(self):
        #*# SCALAR ADDITION ###
        # 0 + 0 = 0
        self.assertEqual(mypy.matrix_add(mypy.cvar.packedID_scalar0,mypy.cvar.packedID_scalar0),
mypy.cvar.packedID_scalar0) 
        
        # 1 + 0 = 1
        self.assertEqual(mypy.matrix_add(mypy.cvar.packedID_scalar1, mypy.cvar.packedID_scalar0),
mypy.cvar.packedID_scalar1)
        
        # 0 + 1 = 1
        self.assertEqual(mypy.matrix_add(mypy.cvar.packedID_scalar0, mypy.cvar.packedID_scalar1),
mypy.cvar.packedID_scalar1)
        
        # 1 + 0 = 0 + 1
        self.assertEqual(mypy.matrix_add(mypy.cvar.packedID_scalar1, mypy.cvar.packedID_scalar0), mypy.matrix_add(mypy.cvar.packedID_scalar0, mypy.cvar.packedID_scalar1))
        
        # 1 + 1 = 2
        #constructing scalar 2 and adding it to the matrix store
# NOTE: the following has been replaced, partly because it's overly complicated,
# mostly because it assumes Complex type 
        # a = np.matrix([[2+0j]])
        # alist = a.reshape(-1).tolist()[0]
        # arr = buildComplexArray(alist)
        # level = 0
        # dim = 2**level
        # aser = mypy.row_major_list_to_store(arr, level, level, dim)
        aser = mypy.get_valID_from_valString("2")
        
        self.assertEqual(mypy.matrix_add(mypy.cvar.packedID_scalar1, mypy.cvar.packedID_scalar1), aser)
        
        #*# SQUARE MATRIX ADDITION ###
        # check that error is thrown when dimensions do not match;
        # (would need to implement Python functionality to throw exception)
        # e.g. self.assertRaises(exception, callable, args, keywords)
        
        # 4x4 identity matrix + 4x4 zero matrix = 4x4 identity matrix
        self.assertEqual(mypy.matrix_add(I2_packedID, Z2_packedID), I2_packedID)
        
        
        #*# NON-SQUARE MATRIX ADDITION ###
        # check that error is thrown when dimensions do not match


    def testScalarMult(self):
        # test scalar multiplication
        # 0 * 0
        self.assertEqual(mypy.scalar_mult(mypy.cvar.packedID_scalar0,mypy.cvar.packedID_scalar0),
mypy.cvar.packedID_scalar0) 
		
        # 1 * 0 = 0
        self.assertEqual(mypy.scalar_mult(mypy.cvar.packedID_scalar1,
mypy.cvar.packedID_scalar0), mypy.cvar.packedID_scalar0)
		
        # 1 * 1 = 1
        self.assertEqual(mypy.scalar_mult(mypy.cvar.packedID_scalar1,
mypy.cvar.packedID_scalar1), mypy.cvar.packedID_scalar1)
        
        # 0 * 4x4 identity matrix = 4x4 0 matrix
        self.assertEqual(mypy.scalar_mult(mypy.cvar.packedID_scalar0,I2_packedID), Z2_packedID)
		
        # 1 * 4x4 identity matrix = 4x4 identity matrix
        self.assertEqual(mypy.scalar_mult(mypy.cvar.packedID_scalar1,I2_packedID), I2_packedID)
		
    def testMult(self):
        # test matrix multiplication
        # 4x4 identity matrix * 4x4 identity matrix = 4x4 identity matrix
        self.assertEqual(mypy.matrix_mult(I2_packedID,I2_packedID), I2_packedID) 	
		
		
#    def testKronecker(self)
#    def testAdjoint(self)
    
        

if __name__ == '__main__':

    print("Running LARC math (example) unit tests...")
    verbose = 1
    mypy.initialize_larc(26,24,10,-1,-1,verbose)
    mypy.create_report_thread(1800)
    print("matrixID for 0 ", mypy.matrixID_from_packedID(mypy.cvar.packedID_scalar0))
    print("matrixID for 1 ", mypy.matrixID_from_packedID(mypy.cvar.packedID_scalar1))
    
        # Define string for using in formating filenames
    scalarTypeStr = mypy.cvar.scalarTypeStr

    # build Zero matrices
    level = 2
    dim_whole = 2**level
    Z2_arr = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    Z2_packedID = mypy.row_major_list_to_store(list(map(str, Z2_arr)), level, level, dim_whole)
    mypy.print_naive(Z2_packedID)
    print("\nmatrixID for the level 2 zero matrix\n", mypy.matrixID_from_packedID(Z2_packedID))

    # build Identity matrices
    I2_arr = [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]
    I2_packedID = mypy.row_major_list_to_store(list(map(str, I2_arr)), level, level, dim_whole)
    mypy.print_naive(I2_packedID)
    print("\nmatrixID for the level 2 identity matrix\n", mypy.matrixID_from_packedID(I2_packedID))
    
    unittest.main()

