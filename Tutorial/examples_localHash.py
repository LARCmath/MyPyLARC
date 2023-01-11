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
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import math

##
# \file examples_localHash.py
#
# \brief This sample program shows the effects of
# LARC's locality hash.

if __name__ == '__main__':

    print("This code tests matrixID interface\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 3
    regionbitparam = 15
    zeroregionbitparam = 20
    verbose = 1
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)
    
    mypy.create_report_thread(1800)

    scalarTypeStr = mypy.cvar.scalarTypeStr

    if scalarTypeStr not in ('Real', 'Complex', 'MPReal', 'MPComplex', 'MPRational', 'MPRatComplex', 'Clifford'):
        raise Exception('Scalar type %s is not supported for this example.' % scalarTypeStr)

    # MAKE a matrix to test Globbing
    # two by two tests
    level = 1
    dim_whole = 2**level

    arr_a = list(map(str,[.4, 0, 0, .3]))
    a_pID = mypy.row_major_list_to_store(arr_a,level,level,dim_whole)
    arr_b = list(map(str,[.8, 0, 0, .6]))
    b_pID = mypy.row_major_list_to_store(arr_b,level,level,dim_whole)
    arr_c = list(map(str,[-.4, 0, 0, -.3]))
    c_pID = mypy.row_major_list_to_store(arr_c,level,level,dim_whole)
    print("The matrixIDs of a, b, and c are %d %d %d\n"
        %(mypy.matrixID_from_packedID(a_pID),
          mypy.matrixID_from_packedID(b_pID),
          mypy.matrixID_from_packedID(c_pID)))
    
    d_pID = mypy.matrix_add(a_pID,a_pID)
    
    print("matrix a:")
    mypy.print_naive(a_pID)
     
    print("\nMatrix id of a + a: %d, should be that of b: %d \n"
        %(mypy.matrixID_from_packedID(d_pID),
          mypy.matrixID_from_packedID(b_pID)))
    if (d_pID == b_pID):
        print("a + a PASSED: \n")
        mypy.print_naive(d_pID)
    else: print("FAILED: [.4,0,0,.3] + [.4,0,0,.3] = [.8,0,0,.6]\n")

    arr_prod1 = list(map(str,[.16, 0, 0, .09]))
    prod1_pID = mypy.row_major_list_to_store(arr_prod1,level,level,dim_whole)
    arr_e = list(map(str,[1, -1, -1, 1]))
    e_pID = mypy.row_major_list_to_store(arr_e,level,level,dim_whole)
    print("\nmatrix e:")
    mypy.print_naive(e_pID)
    
    arr_prod2 = list(map(str,[.4, -.4, -.3, .3]))
    prod2_pID = mypy.row_major_list_to_store(arr_prod2,level,level,dim_whole)
    print("\nThe matrixIDs of prod1, e, and prod2 are %d %d %d\n"
    %(mypy.matrixID_from_packedID(prod1_pID),
      mypy.matrixID_from_packedID(e_pID),
      mypy.matrixID_from_packedID(prod2_pID)))

    m_pID = mypy.matrix_mult(a_pID,a_pID)
    n_pID = mypy.matrix_mult(a_pID,e_pID)
    print("The matrixIDs of m and n are %d %d\n"
        %(mypy.matrixID_from_packedID(m_pID),
          mypy.matrixID_from_packedID(n_pID)))
    print("Matrix id of a * a: %d, should be that of prod1: %d \n" %(
        mypy.matrixID_from_packedID(m_pID),
        mypy.matrixID_from_packedID(prod1_pID)))
    if (m_pID == prod1_pID): 
        print("a * a PASSED:")
        mypy.print_naive(m_pID)
    else: print("FAILED: [.4,0,0,.3] * [.4,0,0,.3] = [.16, 0, 0, .09]\n")

    print("\nMatrix id of a * e: %d, should be that of prod2: %d \n" %(
        mypy.matrixID_from_packedID(n_pID),
        mypy.matrixID_from_packedID(prod2_pID)))
    if (n_pID == prod2_pID): 
        print("a * e PASSED:")
        mypy.print_naive(n_pID)
    else: print("FAILED: [.4,0,0,.3] * [1, -1, -1, 1] = [.4, -.4, -.3, .3]\n")

    # scalar tests
    level = 0
    dim_whole = 2**level

    arr_f = list(map(str,[.4]))
    f_pID = mypy.row_major_list_to_store(arr_f,level,level,dim_whole)
    print("We input the value into f (%d) of .4, the matrix is:"
         %mypy.matrixID_from_packedID(f_pID))
    mypy.print_naive(f_pID)

    arr_g = list(map(str,[.4+1e-5]))
    g_pID = mypy.row_major_list_to_store(arr_g,level,level,dim_whole)
    print("We input the value into g (%d) of .4+1e-5, the matrix is:"
         %mypy.matrixID_from_packedID(g_pID))
    mypy.print_naive(g_pID)

    arr_h = list(map(str,[.4-1e-10]))
    h_pID = mypy.row_major_list_to_store(arr_h,level,level,dim_whole)
    print("We input the value into h (%d) of .4-1e-10, the matrix is:"
         %mypy.matrixID_from_packedID(h_pID))
    mypy.print_naive(h_pID)
