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
#import pylarc
import MyPyLARC as mypy
import numpy as np
import random
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by packedID instead of pointers
	
    print("This code illustrates the use of the bounding scalar types UPPER and LOWER.\n", flush=True)

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)

    # Print out a scalar 0 to synchronize the output buffers.
    mypy.print_naive(mypy.row_major_list_to_store(["0"], 0, 0, 1))
    print("", flush=True)

    print("This code illustrates the use of the bounding scalar types UPPER and LOWER.\n", flush=True)

    # In the Makefile you can compile with different scalarType values
    # Define string for using in formating filenames
    scalarTypeStr = mypy.cvar.scalarTypeStr
    print("The current Scalar Type is: ", scalarTypeStr, flush=True)
    if scalarTypeStr not in ("Upper", "Lower"):
        print("ERROR - This is only useful for bounding types.", flush=True)
        print("Please recompile LARC with TYPE=UPPER or TYPE=LOWER.", flush=True)
        sys.exit(1)

    print("\n\nFIRST EXAMPLE: Powers of a 2x2 column stochastic matrix.", flush=True)

    # Matrix A (first power)
    npa = np.array([[0.75, 0.5], [0.25, 0.5]])
    print("\n-------- Matrix A (from numpy):", flush=True)
    print(npa, flush=True)
    avals = ["0.75", "0.5", "0.25", "0.5"]
    a = mypy.row_major_list_to_store(avals, 1, 1, 2)
    print("\n-------- Matrix A (from LARC):", flush=True)
    mypy.print_naive(a)

    # Matrix B (second power)
    npb = npa @ npa
    print("\n-------- Matrix B = A*A (from numpy):", flush=True)
    print(npb, flush=True)
    b = mypy.matrix_mult(a,a)
    print("\n-------- Matrix B = A*A (from LARC):", flush=True)
    mypy.print_naive(b)

    # Matrix C (third power)
    npc = npb @ npa
    print("\n-------- Matrix C = A*A*A (from numpy):", flush=True)
    print(npc, flush=True)
    c1 = mypy.matrix_mult(b,a)
    print("\n-------- Matrix C1 = B*A (from LARC):", flush=True)
    mypy.print_naive(c1)
    c2 = mypy.matrix_mult(a,b)
    print("\n-------- Matrix C2 = A*B (from LARC):", flush=True)
    mypy.print_naive(c2)


    print("\n\nSECOND EXAMPLE: Probabilities through a Markov chain.", flush=True)

    np.set_printoptions(linewidth=200)
    nptm = np.array([[0, 0.25, 0.25,    0, 0.25,    0, 0.25,    0],
                     [0, 0.25, 0.25,    0, 0.25,    0, 0.25,    0],
                     [0,    0, 0.25, 0.25,    0, 0.25,    0, 0.25],
                     [0, 0.25,    0, 0.25, 0.25,    0, 0.25,    0],
                     [0,    0, 0.25,    0, 0.25, 0.25,    0, 0.25],
                     [0, 0.25,    0, 0.25,    0, 0.25, 0.25,    0],
                     [0,    0, 0.25,    0, 0.25,    0, 0.25, 0.25],
                     [0, 0.25,    0, 0.25,    0, 0.25,    0, 0.25] ])
    print("\n-------- Markov Matrix (from numpy):", flush=True)
    print(nptm, flush=True)

    vals = [ "0", "0.25", "0.25",    "0", "0.25",    "0", "0.25",    "0",
             "0", "0.25", "0.25",    "0", "0.25",    "0", "0.25",    "0",
             "0",    "0", "0.25", "0.25",    "0", "0.25",    "0", "0.25",
             "0", "0.25",    "0", "0.25", "0.25",    "0", "0.25",    "0",
             "0",    "0", "0.25",    "0", "0.25", "0.25",    "0", "0.25",
             "0", "0.25",    "0", "0.25",    "0", "0.25", "0.25",    "0",
             "0",    "0", "0.25",    "0", "0.25",    "0", "0.25", "0.25",
             "0", "0.25",    "0", "0.25",    "0", "0.25",    "0", "0.25" ]
    tm = mypy.row_major_list_to_store(vals, 3, 3, 8)
    print("\n-------- Markov Matrix (from LARC):", flush=True)
    mypy.print_naive(tm)

    npv = np.array([[1, 0, 0, 0, 0, 0, 0, 0]])
    print("\n-------- Initial State Vector (from numpy):", flush=True)
    print(npv[0], flush=True)

    v = mypy.row_major_list_to_store(["1"] + 7*["0"], 0, 3, 8)
    print("\n-------- Initial State Vector (from LARC):", flush=True)
    mypy.print_naive(v)

    for i in range(1, 15):
        npv = npv @ nptm
        print("\n-------- Round #{0} State Vector (from numpy):".format(i), flush=True)
        print(npv[0], flush=True)
        v = mypy.matrix_mult(v,tm)
        print("\n-------- Round #{0} State Vector (from LARC):".format(i), flush=True)
        mypy.print_naive(v)

    print("\nConverging towards 1/7 = {0}, which is LARC value:".format(1/7), flush=True)
    oneseventh = mypy.row_major_list_to_store(["1/7"], 0, 0, 1)
    mypy.print_naive(oneseventh)
    print("", flush=True)


