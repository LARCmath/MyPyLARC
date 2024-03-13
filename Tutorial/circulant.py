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
import MyPyLARC as mypy
import numpy as np
import random
from ctypes import *

##
# \file circulant.py
#
# \brief This sample program creates a random circulant matrix.
#
# Circulant matrices compress nicely using LARC.

if __name__ == '__main__':
	# This version references matrices by packedID instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value


    #*  EXPLAIN VERBOSITY
    # verbose = 0 # SILENT
    verbose = 1 # BASIC
    # verbose = 2 # CHATTY
    # verbose = 3 # DEBUG
    # verbose = 4 # ALL
    
    if (verbose > 1):
        mypy.create_report_thread(1800)

    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)

    scalarTypeStr = mypy.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    # parameters for entering the python array into the store
    level = 7
    dim_whole = 2**level

    if scalarTypeStr in ('Integer', 'MPInteger'):
        randVals = [ random.randrange(0,10001) for i in range(dim_whole)]
    elif scalarTypeStr == 'Boolean':
        randVals = [ random.randrange(0,2) for i in range(dim_whole)]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        randVals = [ np.complex(random.random(),random.random()) for i in range(dim_whole)]
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational', 'Clifford', 'Upper', 'Lower'):
        randVals = [ random.random() for i in range(dim_whole)]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # create circulant matrix from the dim_whole random numbers
    a = []
    b = randVals
    for i in range(dim_whole):
        a.append(list(b))
        b.insert(0,b.pop(-1))
    #print("a: ", a)
    amat = np.matrix(a)
    serial = mypy.add_numpy_matrix_to_matrix_store(amat)
    filename_json = "Data/Out/circulant.lev%d.%s.json" %(level,scalarTypeStr)
    # mypy.print_naive(serial)
    # print("\n")
    mypy.fprint_larcMatrixFile(serial,os.path.join(os.path.dirname(__file__),filename_json))
