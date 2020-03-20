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


##########################################################################
## test_infoStore1 writes a small LARCMatrix file with metadata to the current
## directory
##########################################################################

from __future__ import print_function

import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *


if __name__ == '__main__':

    ###################################
    ##     SET THESE PARAMETERS      ##
    ###################################
    log_psize = 3          ##  problem_size is always power of two!

    #################################
    ##    Calculated Parameters    ##
    #################################
    problem_size = 2**(log_psize)  ## level of the matrices (2^level by 2^level matrices)
    max_level = problem_size       ## the maximum size of matrices in the matrix store

    # print("The globbing constants are: SIGHASH %d, ZEROBITTHRESH %d" %(cvar.py_SIGHASH, cvar.py_ZEROBITTHRESH))

    scalarTypeStr = mypy.cvar.scalarTypeStr
 
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
    ## SMALL STORES for working on desktop 
    if log_psize <= 3:    
        matrix_exponent = 22
        op_exponent = 19            
        rnd_sig_bits = -1 # default value
        trunc_to_zero_bits = -1 # default value is 20
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

    #############################
    # create a matrix in python #
    #############################
    if scalarTypeStr in ('Integer', 'MPInteger'):
        a_str = [1, 3, 5, 6,
                 8, 6, 3, 1,
                 -9, 11, 13, 15,
                 16, 13, 12, 10]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a_str = [1+2j, 3+4j, 5+6j, 7+8j,
                 8+7j, 6+5j, 3+4j, 1+2j,
                 9+10j, 11+12j, 13+14j, 15+16j,
                 16+15j, 14+13j, 12+11j, 10+9j]
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational'):
        a_str = [1, 3, 5, 6,
                 8, 6, 3, 1,
                 -9, 11, 13, 15,
                 16, 13, 12, 10]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # turn the string into an array then use the reader for  row major format
    # a_arr = list(map(str,a_str))
    a_arr = mypy.map_to_str(a_str,scalarTypeStr)

    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    serial = mypy.row_major_list_to_store_matrixID(a_arr, level, level, dim_whole)

    #########################
    # add metadata to store #
    #########################
    print("The info store several predefined info_names.")
    print("Each matrix can have a value for each info_name")
    print("stored in the info_store, and these values will")
    print("be printed in the LARCMatrix file as metadata.")
    print("To add more info_names you would need to edit larc/src/larc.h.")
    print("Here are the info_names that are currently defined:")
    mypy.list_info_names();
    info_name = "DATE"
    info_type = mypy.get_info_type_from_string_name(info_name)
    info_data = "20170707"
    print("Loading %s of %s to info_store" %(info_name,info_data))
    fail = mypy.info_set(info_type, serial, info_data)
    if fail:
       print("Unable to write data to info store")
    else:
       print("Wrote data to info_store")

    info_name = "COMPUTER"
    info_type = mypy.get_info_type_from_string_name(info_name)
    info_data = "a226"
    print("Loading %s of %s to info_store" %(info_name,info_data))
    fail = mypy.info_set(info_type, serial, info_data)
    if fail:
       print("Unable to write data to info store")
    else:
       print("Wrote data to info_store")

    info_name = "COMMENT"
    info_type = mypy.get_info_type_from_string_name(info_name)
    info_data = "matrix_exponent is " + str(matrix_exponent) + ", and rnd_sig_bits is default"
    print("Loading %s of %s to info_store" %(info_name,info_data))
    fail = mypy.info_set(info_type, serial, info_data)
    if fail:
       print("Unable to write data to info store")
    else:
       print("Wrote data to info_store")



    # write out LARCMatrix file
    local_name = "Data/Out/infoTest1.json"
    out_name = os.path.join(os.path.dirname(__file__),local_name)
    print("Writing LARCMatrix file with info attached to %s"  %out_name)
    mypy.write_larcMatrix_file_by_matID(serial,out_name)

