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

import numpy
import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
# import pylarc
import MyPyLARC as mypy
from ctypes import *


##
# \file preloading_mults_pi.py
#
# \brief This routine illustrates the way that
# region size for the locality hash
# is used by LARC to minimize the number of scalars
# that are creating by representing any scalar
# in a small region by a single representative.


if __name__ == '__main__':


    if 3 == len(sys.argv):
        # sys.path.append(sys.argv[1])
        # sys.argv = [sys.argv[0]]
        regionbitparam = int(sys.argv[1]) 
        num_mults = int(sys.argv[2])    


    else:
        print("This program needs two integer inputs: r n")
        print("   r:  the size of locality hash regions is 2^(-r)")
        print("   n:  the test is run for i*pi for i=1 to n")
        print("Sample Usage:")
        print("   python3 preloading_mults_pi.py 5 6\n")
        sys.exit()    


    #*#####################################
    #*    Print baseline usage report    ##
    #*#####################################
    # mypy.memory_and_time_report(0, "stdout")
        

    #*##################################################################
    #*    LARC  Initialization of Matrix Store and Operation Stores   ##
    #*##################################################################
    #* The routine initialize_larc() does the following:              ##
    #* * creates the matrix and op stores                             ##
    #* * preloads matrix store with: standard scalars and gates,      ##
    #*   and with all zero, identity, and (integer) Hadamard matrices ##
    #*   left to max matrix size                                      ##
    #*##################################################################
    max_level = 8     
    matrix_exponent = 22
    op_exponent = 19   
    zeroregionbitparam = 50
    # regionbitparam = 30
    verbose = 0
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,
                         zeroregionbitparam,verbose)
    # mypy.create_report_thread(180)

    print("\nFinished creating LARC matrix and op stores and loading basic matrices.")
    # print("stopHogging check to see if program is too large to occur once every 10 minutes.")
    print("\n***************************************")

    #*###################
    # see the code in roots of unity and in fft for other cases which can fail
    # with some values of regionbitparam and zeroregionbitparam
    #*###########

    #*  FIND out type being used
    scalarTypeStr = mypy.cvar.scalarTypeStr

    if ((scalarTypeStr == "Complex") or (scalarTypeStr == "MPRatComplex")):
        var_pi = complex(numpy.pi ,0)
    elif ((scalarTypeStr == "Integer") or (scalarTypeStr == "MPInteger")):
        print("cannot have pi in an integer type!")
        sys.exit()
    else:
        var_pi = numpy.pi
    print("\nUsing numpy.pi we find that pi is")
    print(var_pi)

    matpi_pID = mypy.get_valID_from_valString(mypy.value_to_string(var_pi,scalarTypeStr))
    print("\nWe have loaded pi into the matrix store and it has matrixID")
    print(mypy.matrixID_from_packedID(matpi_pID))

    # var_2pi = complex(2*numpy.pi ,0)
    # var_2pi = 2*numpy.pi
    # print("2 pi is")
    # print(var_2pi)
    #
    # mat2piID = mypy.get_valID_from_valString(mypy.value_to_string(var_2pi,scalarTypeStr))
    # print("We have loaded pi into the matrix store and it has matrixID")
    # print(mat2piID)
    #
    # pIDtwicepi = mypy.matrix_add(matpi_pID,matpi_pID)
    # print("We retrieved pi in python with numpy")
    # print("We loaded this value into the larc matrix store and it had matrixID matpi_pID")
    # print("We asked LARC to add the scalar values of two copies of the given matpi_pID")
    # print("LARC returned to us the matrixID")
    # print(mypy.matrixID_from_packedID(pIDtwicepi))

    n = num_mults

    print("\nCreating array int_pID[i] of matrixIDs for i for i=0 to n=%d" %n)
    int_pID = [0]*(n+1)
    for k in range(n+1):
        int_pID[k] = mypy.get_valID_from_valString(mypy.value_to_string(k,scalarTypeStr))
    print(list(map(mypy.matrixID_from_packedID,int_pID)))

    print("\nCreating array python_pi_i in python of pi*i.")
    python_pi_i = [0]*(n+1)
    for k in range(n+1):
        python_pi_i[k] = var_pi * k
    print(python_pi_i)   

    print("\nCreating array python_pi_i_pID of matrixIDs for pi*i.")
    python_pi_i_pID = [0]*(n+1)
    for k in range(n+1):
        python_pi_i_pID[k] = mypy.get_valID_from_valString(mypy.value_to_string(
            python_pi_i[k],scalarTypeStr))
    print(list(map(mypy.matrixID_from_packedID,python_pi_i_pID)))

    print("\nCalculating larc_pi_i_pID[k] = mat_mult(int_pID[k],python_pi_i_pID).")
    larc_pi_i_pID = [0]*(n+1)
    for k in range(n+1):
        larc_pi_i_pID[k] = mypy.matrix_mult(int_pID[k],python_pi_i_pID[1])
    print(list(map(mypy.matrixID_from_packedID,larc_pi_i_pID)))

    # user_input = input("Press Enter-key to continue.")
    print("\nRetrieving associated scalars from matrix store.")
    larc_pi_i_val = [0]*(n+1)
    for k in range(n+1):
        larc_pi_i_val[k] = mypy.get_readableString_scalar_from_pID_and_coords(larc_pi_i_pID[k],0,0)
    print(larc_pi_i_val)

    if (python_pi_i_pID == larc_pi_i_pID):
        print("\nBoth calculations produce the same matrixIDs.")
        print("If you would like to see a case where matrixIDs differ try the parameters 1 5.")
    else:
        print("\nThe two calculations do not produce identical matrixIDs.")
        print("If you would like to see an example which seems to work try the parameters 8 5.")

    print("\nEven if both calculations produce the same matrixIDs, there can")    
    print("still be a problem.  Try parameters 5 and 22 and you will notice that")
    print("larc tries to store 7*pi (= 21.9911485751) it returns the matrixID") 
    print("previously stored for the integer 22 since these fall in the same")
    print("locality hash region.  This can cause future inaccuracies in calculations.\n")
