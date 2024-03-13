#!/usr/bin/env python3

 #*################################################################
 #                                                                #
 # Copyright (C) 2014-2024, Institute for Defense Analyses        #
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
# import pylarc
import MyPyLARC as mypy
from ctypes import *

# Python Interface version of test_math.py
# Uses packedIDs instead of addresses

##
# \file examples_math.py
#
# \brief This routine illustrates various matrix
# operations on both square and nonsquare matrices.

if __name__ == '__main__':

    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)

    scalarTypeStr = mypy.cvar.scalarTypeStr

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_packedID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix file\n")
    mypy.print_naive(samp_packedID)
    
    print("does scalarM1_val print?")
    scalarM1_val = '-1'
    scalarM1_packedID = mypy.get_valID_from_valString(scalarM1_val)
    mypy.print_naive(scalarM1_packedID)
    
    print("testing scalar_mult:")
    samp2_packedID = mypy.scalar_mult(scalarM1_packedID,samp_packedID)
    mypy.print_naive(samp2_packedID)
    
    print("testing addition:")
    samp3_packedID = mypy.matrix_add(samp_packedID,samp2_packedID)
    mypy.print_naive(samp3_packedID)
    
    print("testing adjoint:")
    samp3_packedID = mypy.adjoint(samp_packedID)
    mypy.print_naive(samp3_packedID)
    
    print("testing non-square matrix mult:")
    samp4_packedID = mypy.matrix_mult(samp3_packedID,samp_packedID)
    mypy.print_naive(samp4_packedID)
    print("")
    samp4_packedID = mypy.matrix_mult(samp_packedID,samp3_packedID)
    mypy.print_naive(samp4_packedID)
    print("testing kron product:")
    samp4_packedID = mypy.kronecker_product(samp_packedID,samp_packedID)
    mypy.print_naive(samp4_packedID)
    print("testing join:")
    samp4_packedID = mypy.join(samp_packedID,samp_packedID)
    mypy.print_naive(samp4_packedID)
    print("testing stack:")
    samp4_packedID = mypy.stack(samp_packedID,samp_packedID)
    mypy.print_naive(samp4_packedID)
