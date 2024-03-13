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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
# import pylarc
import MyPyLARC as mypy
from ctypes import *
import json


##
# \file examples_roots_unity.py
#
# \brief This routine generates the 24th roots of unity w,
# and then checks for the appropriate closure under
# multiplication of the various powers of w.
#
# That is, it verifies that
#    matrixProduct(packedID(w^i),packedID(w^j)) = packedID(w^k)
# where k congruent to i+j mod 24.

if __name__ == '__main__':


    #*###############################
    #*   SET THESE PARAMETERS      ##
    #*###############################
    max_level = 8          #*  problem_size is always power of two!
    verbose = 0

    
    #*####################################################################
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*####################################################################
    scalarTypeStr = mypy.cvar.scalarTypeStr

    if scalarTypeStr not in ('Complex', 'MPComplex', 'MPRatComplex'):
        raise Exception('Scalar type %s is not supported for this example.' % scalarTypeStr)

    #*####################################################
    #*   Find out if machine has a large amount of      ##
    #*   memory available so we can make bigger tables  ##
    #*####################################################
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

    #*#####################################
    #*    Print baseline usage report    ##
    #*#####################################
    if (verbose > 0):
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mypy.memory_and_time_report(0, "stdout")

    #* read the parameter file into a python dictionary
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

    #* warn if the commandline value for verbose differs from the parameter file value for verbose        
    if (verbose > 0):
        if (verbose != p_verbose):
            print("NOTE: This program uses commandline (verbose = %d) " %verbose)
            print("      rather than the parameter file (verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    #* initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)


    #*#######################
    # print roots of unity  #
    #*#######################
    # verbose = 0
    # n = 4
    # print("\nRunning code to produce the %d-th roots of unity." %n)
    # mypy.print_n_th_roots_of_unity(4,verbose)
    # # mypy.print_n_th_roots_of_unity(24,verbose)

    #*##################################################
    # load principal root of unity and return packedID #
    #*##################################################
    # verbose = 0
    # print("\nRunning code to produce the packedID for the n-th principal root of unity.")
    # p1_pID = mypy.principal_n_th_root_of_unity_matID(1, verbose)
    # print("principal nth root of unity for n=1 has packedID %d" %p1_pID)
    # p2_pID = mypy.principal_n_th_root_of_unity_matID(2, verbose)
    # print("principal nth root of unity for n=2 has packedID %d" %p2_pID)
    # p3_pID = mypy.principal_n_th_root_of_unity_matID(3, verbose)
    # print("principal nth root of unity for n=3 has packedID %d" %p3_pID)
    # p4_pID = mypy.principal_n_th_root_of_unity_matID(4, verbose)
    # print("principal nth root of unity for n=4 has packedID %d" %p4_pID)
    # # p24_pID = mypy.principal_n_th_root_of_unity_matID(24, verbose)
    # # print("principal nth root of unity for n=24 has packedID %d" %p24_pID)


    #*###############################################################
    # fill an array with the packedIDs of the nth roots of unity   #*
    # verbose = 0 is run quiet, verbose = 1 subroutines chatty, and  #*
    # verbose = 2 both local and subroutine calls chatty as well.   #*
    #*###############################################################
    #* verbose = 2
    verbose = 1
    if (verbose>1):
        call_verbose = 1
    else:    
        call_verbose = 0
    print("PRELOADING FOR IDENTITY PRESERVATION:")
    print("You may preload a set of scalars which you desire to act mathematically")
    print("within LARC as though they satisfy one or more identities.")
    print("For example, say you want multiplication to be closed for roots of unity.")
    print("This is equivalent to saying we want the function")
    print("  packedID = mypy.matrix_mult(rootj_ID,rootk_ID)")
    print("to return packedID = root(j+k mod n)_ID.")     
    print("When you load a scalar immediately after initialization, it becomes the")
    print("unique scalar value from its LSH region in LARC unless it happens to")
    print("lie in a region which already has a stored representative.")
    print("Two scalars that we care to differentiate should not end up in the same")
    print("LSH region of LARC.  You should test the identities that you care about,")
    print("and if necessary, alter your hash parameters to correct the")
    print("situation.\n")
    user_input = input("Press return to see example.")

    n = 24
    print("EXAMPLE:")
    print("Running code to produce a list of matrixIDs for all the %d-th roots of unity." %n)
    n_array = [0]*(n)
    if (verbose>1):
        print("The length of the array is %d" %len(n_array))
    for k in range(n):
        n_array[k]  = mypy.k_th_power_of_n_th_root_of_unity_pID(k,n,call_verbose)
    print("\nHere is the array of the LARC matrixIDs for the %dth roots of unity:" %n)   
    print(list(map(mypy.matrixID_from_packedID,n_array)))
    print(" ")
    user_input = input("Press return, then we will print the stored values.")
    
    if (verbose>0):
        print("\nThe stored values of these roots are:")
        print("")
        for k in range(n):
           mypy.print_naive(n_array[k])


    print("\nNow we verify that multiplication is closed and correct, by looking")
    print("at the matrixIDs of all products of pairs of these roots and confirming")
    print("that the results have the expected matrixIDs of the appropriate roots.")
    success = 1
    for k in range(n):
        for j in range(k+1):
            my_packedID = mypy.matrix_mult(n_array[j],n_array[k])
            flag = 0  # initialize as though there was a closure failure
            for i in range(n):
                if (my_packedID == n_array[i]):
                    flag = 1   # we found the matrixID of product in the preloaded roots
            if (flag == 0):
                success = 0
                print("\nMatrix multiplication of the %d-th roots of unity is not closed" %n)
                print("since the product of the %d-th power and %d-th power" %(j,k,))
                m = (j+k) % n
                print("which should have been the %d-th power which had value" %m )
                mypy.print_naive(n_array[m])

                print("The computed result has matrixID %g instead of %g" %(mypy.matrixID_from_packedID(my_packedID),
                                                                            n_array[m]))
                print("Which differs from the expected value by")
                # complex_dif = 0.1 + I*0.0
                # print(complex_dif)

    if success:
        print("\nPassed test! All products of pairs of roots of unity produce existing roots.\n")
    else:
        print("\nFailed test! At least one product of roots of unity produced a value")
        print("that was not a preloaded root of unity, so locality hash is not optimal.\n")
        

    # sys.exit(0)

    #* OLDER VERSION OF CODE
    # verbose = 1
    # if (verbose>1):
    #     call_verbose = 1
    # else:    
    #     call_verbose = 0
    # n = 5
    # print("\nRunning code to produce a list of matrixIDs for all the %d-th roots of unity." %n)
    # n5_array = [0]*(n)
    # if (verbose>1):
    #     print("The length of the array is %d" %len(n5_array))
    # for k in range(n):
    #     n5_array[k]  = mypy.k_th_power_of_n_th_root_of_unity_pID(k,n,call_verbose)
    # print("\nHere is the array of matrixIDs of the %dth roots of unity:" %n)
    # print(n5_array)
    # if (verbose>1):
    #     print("\nThe stored values of these roots are:")
    #     print("")
    #     for k in range(n):
    #        print("The matrixID is %d" %n5_array[k])
    #        print("The %d-th power of the %d-th root of unity is" %(k,n))
    #        mypy.print_naive(n5_array[k])
    #        print("")
    # print("\nNow we can look to see if multiplication is closed, by looking")
    # print("at the matrixIDs of products of pairs of these roots.")
    # for k in range(n):
    #     for j in range(k+1):
    #         my_matrixID = mypy.matrix_mult(n5_array[j],n5_array[k])
    #         print("The (%d,%d)th product has matrix ID" %(j,k))
    #         print(my_matrixID)
    # print("\n")

    # print("\nTODO: Fix these non desirable:")
    # print("\nFor n = 5 we see that the product of the 1st and 4th roots,")
    # my_matrixID = mypy.matrix_mult(n5_array[1],n5_array[4])
    # print("the (%d,%d)th product has matrix ID" %(1,4))
    # print(my_matrixID)
    # print("which has stored value")
    # mypy.print_naive(my_matrixID)
    # print("")
    # print("Whereas if the locality hash was perfect we expected to see")
    # print(n5_array[0])
    # print("which has stored value")
    # mypy.print_naive(n5_array[0])


    # pID = mypy.matrix_diff(int64_t A_pID, int64_t B_pID)
    # get numerical value pID
    # if val
    # i = 1
    # while value > 1/(2^i)for i in 
    #     see if the value is < 
    #         my_matrixID = mypy.matrix_mult(n5_array[j],n5_array[k])
    #         print("The (%d,%d)th product has matrix ID" %(j,k))
    #         print(my_matrixID)
    # print("\n")
