#!/usr/bin/env python3

 ##################################################################
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
 ##################################################################

from __future__ import print_function

import os
import glob
import sys
import random # needed for toeplitz
import numpy as np # needed for toeplitz
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
# import pylarc
import MyPyLARC as mypy
from ctypes import *

if __name__ == '__main__':

    #################################
    ##   SET THESE PARAMETERS      ##
    #################################
    max_level = 8          ##  problem_size is always power of two!


    ####################################################
    ##   Find out if machine is desktop workstation   ##
    ##   or a CPU-cycle servers (cs1-cs6)             ##
    ####################################################
    machine = os.uname()[1]
    cs = 0        # on desktop workstation, with smaller memory
    if (machine.find('cs') >= 0):
        cs = 1    # on CPU-cycle server cs1-cs6, with larger memory
        print("This machine is a CPU-cycle server")
    else:
        print("This machine is a desktop work station")


    #######################################
    ##    Print baseline usage report    ##
    #######################################
    mypy.rusage_report(0, "stdout")


    ####################################################################
    ##    LARCt Initialization of Matrix Store and Operation Stores   ##
    ####################################################################
    ## The routine initialize_larc() does the following:              ##
    ## * creates the matrix and op stores                             ##
    ## * preloads matrix store with: standard scalars and gates,      ##
    ##   and with all zero, identity, and (integer) Hadamard matrices ##
    ##   left to max matrix size                                      ##
    ####################################################################

    ####################################################################
    ##    Testing to see what approximation functions.
    ##    The zerobit thresh parameters that work for F_3 are:
    ##      -z 53 and smaller
    ##    The zerobit thresh parameters that fail for F_3 are:
    ##      -z 54 and larger  
    ## 
    ##    The rounding function does not effect whether it
    ##    works in ranges -s 10 to -s 1000
    ## 
    ## 
    ####################################################################


    ## SMALL STORES for working on desktop
    if max_level <= 8:    
        matrix_exponent = 22
        op_exponent = 19   


        trunc_to_zero_bits = 54
        rnd_sig_bits = 1000

        ######################################################
        ##  Sample failure values for LARC approximation    ##
        ##  are to set         rnd_sig_bits = 60            ##
        ##  and                trunc_to_zero_bits = 60      ##
        ######################################################
        ##  Default values for LARC approximation are       ##
        ##  both equal to DBL_MANT_DIG -2                   ##
        ##  rnd_sig_bits = -1    # default is 53 bits       ##
        ##  trunc_to_zero_bits = -1 # OLD default is 1074 bits  ##
        ##  trunc_to_zero_bits = -1 # OLD default is 1074 bits  ##
        ##  NOTE:  testing shows -z 47 will work            ##
        ######################################################
        # trunc_to_zero_bits = 52


        ######################################################

        ##  TODO: find out the space in which this test fails!!!

        ##  DBL_MANT_DIG is the number of digits in FLT_MANT  ##
        ##  why aren't we using DBL_MANT_BITS  ??????? the number of bits
        ##  used in the mantissa


        ######################################################
        verbose = 1
        mypy.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)
        mypy.create_report_thread(180)
        print_naive = 0
        print_nonzeros = 0
        print("Problem size is small enough to run on desktop")
        if print_naive:
           print("  will print files of naive matrices")
        else: 
           print("  not printing files of naive matrices")
        if print_nonzeros:
           print("  will print files of nonzero matrices\n")
        else: 
           print("  not printing files of nonzero matrices\n")
    ## LARGE STORES for cs1l,cs4l,cs9l
    else:      
        ## matrix_exponent = 26
        ## op_exponent = 24
        matrix_exponent = 30
        op_exponent = 31   
        rnd_sig_bits = -1 # default value
        trunc_to_zero_bits = -1 # default value
        ## trunc_to_zero_bits = 20 # truncate to zero if value is less than 2**(-threshold)
        ## trunc_to_zero_bits = 16 # truncate to zero if value is less than 2**(-threshold)
        verbose = 1
        mypy.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)
        mypy.create_report_thread(3600)   # once per hour 3600
        print_naive = 0      
        print_nonzeros = 0
        print("Problem size is NOT small enough to run on desktop")
        if print_naive:
           print("  WARNING: will try to print files of naive matrices!!!")
        else: 
           print("  not printing files of naive matrices")
        if print_nonzeros:
           print("  WARNING: will print files of nonzero matrices!!!\n")
        else: 
           print("  not printing files of nonzero matrices\n")

    print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
    print("Seppuku check to see if program is to large to occur once every 10 minutes.\n")


    
    # Define string for using in formating filenames
    scalarTypeStr = mypy.cvar.scalarTypeStr


    ##  START OF CODE STOLEN FROM toeplitz.py

    # build array in C from Python list of scalars
    # print("Using row_major_list_to_store on data entered from python\n")

    # parameters for entering the python array into the store
    level = 8
    dim_whole = 2**level

    if scalarTypeStr in ('Real', 'MPReal', 'MPRational'):
        randVals = [ random.random() for i in range(2*dim_whole-1)]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        randVals = [ np.complex(random.random(),random.random()) for i in range(2*dim_whole-1)]
    elif scalarTypeStr in ('Integer', 'MPInteger'):
        randVals = [ random.randrange(0,10001) for i in range(2*dim_whole-1)]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # create toeplitz matrix from the 2*dim_whole-1 random numbers
    a = []
    b = randVals
    # print("b: ", b)
    for i in range(dim_whole):
        a.append(list(b[dim_whole-i-1:2*dim_whole-i-1]))
    # print("a: ", a)
    amat = np.matrix(a)
    alist = amat.reshape(-1).tolist()[0]
    #print alist
    arr = mypy.map_to_str(alist, scalarTypeStr)
    # print 'arr:', mypy.str_scalarTypeArray(arr, len(alist))

    # creating or finding the matrix associated with the array
    top_ID = mypy.row_major_list_to_store_matrixID(arr, level, level, dim_whole)
    # filename_json = "../dat/out/toeplitz.lev%d.%s.json" %(level,scalarTypeStr)
    # mypy.print_naive_by_matID(serial)
    # print("\n")
    # mypy.write_larcMatrix_file_by_matID(serial,os.path.join(os.path.dirname(__file__),filename_json))

    ##  END OF CODE STOLEN FROM toeplitz.py

    # some operations to put in the store
    A_ID = mypy.matrix_add_matrixID(top_ID,top_ID)
    B_ID = mypy.matrix_mult_matrixID(A_ID,top_ID)
    



    # #############################
    # # test reporting
    # #############################
    # mypy.matrix_store_report("../dat/out/temp.mat_report")
    # mypy.op_store_report("../dat/out/temp.op_report")
    # mypy.rusage_report(0,"../dat/out/temp.rusage0_report")
    # mypy.rusage_report(1,"../dat/out/temp.rusage1_report")

    print("\n")
    mypy.matrix_store_report("stdout")
    print("\n")
    mypy.op_store_report("stdout")
    print("\n")
    mypy.rusage_report(0,"stdout")
    print("\n")
    mypy.rusage_report(1,"stdout")
