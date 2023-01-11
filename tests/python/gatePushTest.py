#!/usr/bin/env python3

 #*##############################################################*#
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
 # Additional contributors are listed in "LARCContributors".      #
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
 #*##############################################################*#


from __future__ import print_function, division

import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *
import datetime
import json

## \file practiceGates.py
#
#  \brief Check some matrix build ops, larcMatrix write, and compare work.
#
# This program creates some matrices using routines from MyPyLARC/src/gates.c.
# It is a fairly simple demonstration of how to create a reversible circuit
# using the universal gate library for reversible computing
# (NOT, CNOT, CCNOT) and calculating results from such a circuit.
#
if __name__ == '__main__':

    verbose = 0
    debug = 0

    #*################################*#
    #*    Basic Parameter Setting     *#
    #*################################*#

    #*##################################################################*#
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*##################################################################*#
    scalarTypeStr = mypy.cvar.scalarTypeStr

    matrix_exponent = 17
    op_exponent = 15
    max_level= 8
    regionbitparam = -1
    zeroregionbitparam = -1
    report_interval_seconds = 3600
    min_memGiB_required = 10
    p_verbose = 0

    # initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)

    half_level = max_level//2     


    #*#########################################*#
    #*    Output paths and filename suffix     *#
    #*#########################################*#
    # make directory to hold output, if it doesn't already exist
    # (if parent directories don't exist, will make them too)
    # if directory already exists, print error message and abort

    output_path = "temp"

    #   TODO: Add some way to fail correctly for a make before pushing test
    if not os.path.isdir(output_path):
        os.umask(7) # corresponds to "chmod 770"
        os.makedirs(output_path)
    else:
        # delete old files so new reports are not appended to old ones
        files=glob.glob(output_path+'/circuit*')
        for f in files:
            os.remove(f)          
            
    # create base suffix for file names
    file_suffix = "_gates%d" %max_level

    # Get a large identity matrix
    Imax_mID = mypy.get_identity_pID(max_level)

    # pick offsets for CCnot and Control Not Gates
    CntrlNot_offset = 2  # offset for the Control Not gate
    CCnotA_offset = 3  # offset for CCnot gate A
    CCnotB_offset = 4  # offset for CCnot gate B
              
    circuit_matrix1_mID = Imax_mID
    for i in range(half_level):
        # CCnot gates, with controls on control wires circularly shifted
        # by CCnotB_offset and CCnotA_offset
        circuit_matrix1_mID = mypy.matrix_mult(circuit_matrix1_mID,
           mypy.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
        # Cnot gates with control on control wires circularly shifted by CntrlNot_offset
        circuit_matrix1_mID = mypy.matrix_mult(circuit_matrix1_mID,
           mypy.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
    base_name = "/circuit_matrix1"+file_suffix+".json"
    out_name1 = output_path+base_name
    mypy.fprint_larcMatrixFile(circuit_matrix1_mID,out_name1)

    # Clean out the matrix store of everything but the basic preloaded matrices
    mypy.clean_matrix_storage()


    # Start second circuit (which we hope will end up being identical matrix)
    circuit_matrix2_mID = Imax_mID

    for i in range(half_level):
        # Cnot gates with control on control wires circularly shifted by CntrlNot_offset
        circuit_matrix2_mID = mypy.matrix_mult(circuit_matrix2_mID,
           mypy.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
        # CCnot gates, with controls on control wires circularly shifted
        # by CCnotB_offset and CCnotA_offset
        circuit_matrix2_mID = mypy.matrix_mult(circuit_matrix2_mID,
           mypy.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
    base_name = "/circuit_matrix2"+file_suffix+".json"
    out_name2 = output_path+base_name
    mypy.fprint_larcMatrixFile(circuit_matrix2_mID,out_name2)

    # We verify that two circuits produce identical matrices by
    # reading in the larcMatrix file from the stored circuit
    # to get its assigned matrixID, and compare to the second ID
    circuit_matrix3_mID = mypy.read_larcMatrixFile(out_name1)

    scalarTypeStr = mypy.cvar.scalarTypeStr
    print("%s scalarType: Test for whether commuting circuits have same IDs:" 
          %scalarTypeStr)

    if (circuit_matrix3_mID == circuit_matrix2_mID):
        print("  PASSED.")
    else:
        print("  FAILED.")

    # Delete files and directory we created
    files=glob.glob(output_path+'/circuit*')
    for f in files:
        os.remove(f)
    os.rmdir(output_path)        
        
