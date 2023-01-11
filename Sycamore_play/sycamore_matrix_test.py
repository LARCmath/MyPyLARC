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
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *
import datetime
import json

## \file practiceGates.py
#
#  \brief Demonstrates creation of reversible circuits.
#
# This program creates some matrices using routines from MyPyLARC/src/gates.c.
# It is a fairly simple demonstration of how to create a reversible circuit
# using the universal gate library for reversible computing
# (NOT, CNOT, CCNOT) and calculating results from such a circuit.
#
if __name__ == '__main__':

    verbose = 1
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

    if scalarTypeStr not in ('Complex','MPComplex','MPRatComplex','Clifford'):
        print("This program requires a complex scalarType - please recompile")
        print("with TYPE = one of COMPLEX, MPCOMPLEX, MPRATCOMPLEX, CLIFFORD.")
        sys.exit(0)

    #*##################################################*#
    #*   Find out if machine has a large amount of      *#
    #*   memory available so we can make bigger tables  *#
    #*##################################################*#
    memory_available = mypy.memory_available_GiB()
    if (verbose > 0):
        print("\nThe memory available is %ld GiB" %memory_available)
        print("We will use this to select which computing_env to read from parameter file.")
        print("You could write code to select computing_env automatically.")

    if (memory_available > 200):
        if (verbose > 0):
            print("\nThis memory is more than 200 GiB\n")
        computing_env = 'large'
    else:    
        if (memory_available > 50):
            if (verbose > 0):
                print("\nThis memory is between 50 and 200 GiB\n")
            computing_env = 'medium'
        else:
            if (verbose > 0):
                print("\nThis memory is less than 50 GiB\n")
            computing_env = 'small'
            
    if (verbose > 0):
        print("This program believes the computing_environment is %s" %computing_env)

    #*###################################*#
    #*    Print baseline usage report    *#
    #*###################################*#
    if (verbose > 0):
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mypy.memory_and_time_report(0, "stdout")

    # read the parameter file into a python dictionary
    #with open('../InitParams/REPLACE_THIS.init_params','r') as init_file:
    with open('../InitParams/tutorial.init_params','r') as init_file:
        init_param = json.load(init_file)
        for p in init_param[computing_env]:
            if (verbose > 1):
                print('MatrixExponent: %d' %(p['matrix_exponent']))
                print('OpExponent: %d' %(p['op_exponent']))
                print('MaxLevel: %d' %(p['max_level']))
                print('RegionBitParam: %d' %(p['regionbitparam']))
                print('ZeroRegionBitParam: %d' %(p['zeroregionbitparam']))
                print('ReportIntervalSecs: %d' %(p['report_interval_seconds']))
                print('MinMemRequiredGiB: %d' %(p['min_memGiB_required']))
                print('Verbose: %d' %(p['verbose']))
                print('')
            matrix_exponent = p['matrix_exponent']
            op_exponent = p['op_exponent']
            max_level= p['max_level']
            regionbitparam = p['regionbitparam']
            zeroregionbitparam = p['zeroregionbitparam']
            report_interval_seconds = p['report_interval_seconds']
            min_memGiB_required = p['min_memGiB_required']
            p_verbose = p['verbose']

    # warn if the commandline value for verbose differs from the parameter file value for verbose        
    if (verbose > 0):
        if (verbose != p_verbose):
            print("NOTE: This program uses commandline (verbose = %d) " %verbose)
            print("      rather than the parameter file (verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    # initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)
    mypy.create_report_thread(report_interval_seconds)
    if (verbose):
        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("Seppuku check to see if program is too large to occur once every %d seconds.\n" %report_interval_seconds)

    half_level = max_level//2     


    
    #*#########################################*#
    #*    Matrices                             *#
    #*#########################################*#
    
    # parameters for entering the one qubit gates
    level = 1
    dim_whole = 2**level

    # create entries for gate_sqrtX
    sqrtX_list = [.5+.5j, .5-.5j,
                 .5-.5j, .5+.5j]
    sqrtX_str_list = mypy.map_to_str(sqrtX_list,scalarTypeStr)
    sqrtX_ID = mypy.row_major_list_to_store(
        sqrtX_str_list, level, level, dim_whole)
    print("sqrtX (from Python):")
    mypy.print_naive(sqrtX_ID)
    print()
    c_sqrtX_ID = mypy.build_sycamore_gate_sequence(["sqrtX"], 1)
    print("sqrtX (from C):")
    mypy.print_naive(c_sqrtX_ID)
    print("\n")

    # create entries for gate_sqrtY
    sqrtY_list = [.5+.5j, -.5-.5j,
                 .5+.5j, .5+.5j]
    sqrtY_str_list = mypy.map_to_str(sqrtY_list,scalarTypeStr)
    sqrtY_ID = mypy.row_major_list_to_store(
        sqrtY_str_list, level, level, dim_whole)
    print("sqrtY (from Python):")
    mypy.print_naive(sqrtY_ID)
    print()
    c_sqrtY_ID = mypy.build_sycamore_gate_sequence(["sqrtY"], 1)
    print("sqrtY (from C):")
    mypy.print_naive(c_sqrtY_ID)
    print("\n")

    # create entries for gate_sqrtW using panel method on panels A,B,C,D
    # bottom left =  sqrt(1/2) 
    panelC_ID = mypy.get_pID_for_enum_const(mypy.SCALAR_ENUM_INV_SQRT2)
    # top right  = - sqrt(1/2) * I
    negI_ID = mypy.cvar.packedID_scalar0iM1
    panelB_ID = mypy.scalar_mult(negI_ID,panelC_ID)
    # top left  and bottom right = .5 + .5 j 
    #panelA_ID = mypy.get_valID_from_valString(".5+I*.5")
    panelA_ID = mypy.get_valID_from_valString(mypy.value_to_string(0.5+0.5j, scalarTypeStr))
    panelD_ID = panelA_ID
    sqrtW_ID = mypy.get_pID_from_four_sub_pIDs(
        panelA_ID,panelB_ID,panelC_ID,panelD_ID,1,1)
    print("sqrtW (from Python):")
    mypy.print_naive(sqrtW_ID)
    print()
    c_sqrtW_ID = mypy.build_sycamore_gate_sequence(["sqrtW"], 1)
    print("sqrtW (from C):")
    mypy.print_naive(c_sqrtW_ID)
    print("\n")

    # Do products of single-qubit gates.
    two_single_ID = mypy.build_sycamore_gate_sequence(["sqrtW","sqrtW"], 2)
    print("sqrtW tensor sqrtW:")
    mypy.print_naive(two_single_ID)
    print("\n")

    # ORIGINAL COMMENT THAT I NOW THINK IS WRONG:
    # supreme gate is 4 by 4 with the entries
    # 1  0 0 0
    # 0  0 i 0
    # 0 -i 0 0
    # 0  0 0 sqrt(.75)-.5*i
    # was the last number supposed to be  e^ (-i \pi /6) or e^ (-i \pi /12)?

    # NEW COMMENT THAT I NOW THINK IS RIGHT:
    # supreme gate is 4 by 4 with the entries
    # 1 0 0 0
    # 0 0 i 0
    # 0 i 0 0
    # 0 0 0 sqrt(.75)+.5*i
    # where the last number is e^(i\pi / 6) = (-1)^(1/6)

    if (verbose > 1):
        call_verbose = 1
    else:
        call_verbose = 0
    root3_ID = mypy.get_pID_for_enum_const(mypy.SCALAR_ENUM_SQRT3)
    root3plusI_ID = mypy.matrix_add(root3_ID, mypy.cvar.packedID_scalar0i1)
    omega_ID = mypy.scalar_mult(mypy.cvar.packedID_scalar0_5, root3plusI_ID)
    print("omega:"); mypy.print_naive(omega_ID)
    print()
    omega2_ID = mypy.k_th_power_of_n_th_root_of_unity_pID(1,12,call_verbose)
    print("omega2:"); mypy.print_naive(omega2_ID)
    print("\n")

    zero_ID = mypy.get_valID_from_valString("0")
    one_ID = mypy.get_valID_from_valString("1")
    imag_ID = mypy.get_valID_from_valString("0+I*1")
    supreme_panelA_ID = mypy.get_pID_from_four_sub_pIDs(
        one_ID, zero_ID, zero_ID, zero_ID, 1, 1)
    supreme_panelB_ID = mypy.get_pID_from_four_sub_pIDs(
        zero_ID, zero_ID, imag_ID, zero_ID, 1, 1)
    supreme_panelC_ID = mypy.get_pID_from_four_sub_pIDs(
        zero_ID, imag_ID, zero_ID, zero_ID, 1, 1)
    supreme_panelD_ID = mypy.get_pID_from_four_sub_pIDs(
        zero_ID, zero_ID, zero_ID, omega_ID, 1, 1)
    supreme_ID = mypy.get_pID_from_four_sub_pIDs(
        supreme_panelA_ID, supreme_panelB_ID,
        supreme_panelC_ID, supreme_panelD_ID, 2, 2)
    print("supreme (from Python):")
    mypy.print_naive(supreme_ID)
    print("\n")

    # Using the C code to build the supremacy matrix instead
    supreme2_ID = mypy.build_sycamore_2gate(0, 1, 2)
    print("supreme2 (from C):")
    mypy.print_naive(supreme2_ID)
    print("\n")

