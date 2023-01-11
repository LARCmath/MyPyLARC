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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
from ctypes import *
import json    #* for loading parameter files


##
# \file fft_paramFile_init.py
#
# \brief This illustrates how to initialize LARC
# by running code that checks the available memory
# then selecting an appropriate set of parameters
# from a stored file to use in initializing LARC.
#
# The routine then has an example where matrices
# are generated that have components in
# a block sparse representation of the discrete
# Fourier transform.


if __name__ == '__main__':

    verbose = 0


    #*##################################################
    #*   Print description sparse block algorithm     ##
    #*   from Cooley-Tukey Radix-2 Factorization      ##
    #*   see Van Loan, Computational Frameworks for   ##
    #*   the Fast Fourier Transform, p.21             ##
    #*##################################################
    print("\nWe will create the matrices used in the Cooley-Tukey")
    print("radix-2 sparse block recursive FFT that is from Van Loan,")
    print("Computational Frameworks for the Fast Fourier Transform. p.21")
    if verbose:
        print("The level k, 2**k by 2**k Fourier matrix F_k can be")
        print("generated recursively by the equation")
        print("(subscripts represent levels, not size=2^level):")
        print("   F_k = C_k * (I_1 @ F_(k-1) ) * PI_k")
        print("where @ is used to indicate the Kronecker product,")
        print("   PI_k is the 2^k by 2^k inverse shuffle matrix ")
        print("     note: PI_k's transverse is its inverse and ")
        print("           would shuffle a column vector")
        print("     examples:")
        print("       PI_2 = (1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1)")  # swap 
        print("       PI_3 = (10000000; 00100000; 00001000; 00000010; ")
        print("               01000000; 00010000; 00000100; 00000001)")
        print("   C_k is the matrix constructed by having quadrant ")
        print("       submatrices [UL, UR; LL, LR] = ")
        print("       [I_(k-1), D_(k-1); I_(k-1), -D_(k-1)] ")   
        print("   D_k is the diagonal matrix which had ")
        print("       1, w, w^2, ....w^(2^k - 1) on the diagonal,") 
        print("       where w is the root of unity cooresponding")
        print("       to the next largest size fourier transform, ")
        print("       that is the 2^(k+1)th root unity.")
        print("       e.g. D_1 = (1, 0: 0 , i);")
        print("We end up getting a formula for the FFT in terms of")
        print("a product of block matrices:")
        print("F_k = (product_for j= 0 to k-1) I_j @ C_(k-j) * I_k @ F_0")
        print("       * (product_for j= k-1 to 0) I_j @ PI_(k-j)")
        print("Since F_0 = I_0 and PI_1 = I_1, two of these terms are identity")
        print("matrices (also interesting to note:  C_1 = H), so we have:")
        print("F_k = (product_for j= 0 to k-1) I_j @ C_(k-j) * ")
        print("      (product_for j= k-2 to 0) I_j @ PI_(k-j) ")
        print("where @ is the tensor product, and * is matrix multiply.")


    #*####################################################################
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*####################################################################
    scalarTypeStr = mypy.cvar.scalarTypeStr

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

    print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
    print("stopHogging check to see if program is too large to occur once every 10 minutes.\n")


    #*##############################
    # inverse permutation matrices #
    #*##############################
    print("\nPI_0 matrix is:")
    PI_0 = mypy.create_invShufMat(0)
    mypy.print_naive(PI_0)

    print("\nPI_1 matrix is:")
    PI_1 = mypy.create_invShufMat(1)
    mypy.print_naive(PI_1)

    print("\nPI_2 matrix is:")
    PI_2 = mypy.create_invShufMat(2)
    mypy.print_naive(PI_2)

    print("\nPI_3 matrix is:")
    PI_3 = mypy.create_invShufMat(3)
    mypy.print_naive(PI_3)

    # print("\nPI_4 matrix is:")
    # PI_4 = mypy.create_invShufMat(4)
    # mypy.print_naive(PI_4)


    #*#######################
    # print roots of unity  #
    #*#######################
    print("\n")
    mypy.print_pow2_roots_unity(1)
    mypy.print_pow2_roots_unity(2)
    mypy.print_pow2_roots_unity(3)


    #*#############################
    # create D matrices in python #
    #*#############################
    print("\nD_0 matrix is:")
    D_0 = mypy.create_FFT_DMat(0)
    mypy.print_naive(D_0)

    print("\nD_1 matrix is:")
    D_1 = mypy.create_FFT_DMat(1)
    mypy.print_naive(D_1)

    print("\nD_2 matrix is:")
    D_2 = mypy.create_FFT_DMat(2)
    mypy.print_naive(D_2)

    print("\nD_3 matrix is:")
    D_3 = mypy.create_FFT_DMat(3)
    mypy.print_naive(D_3)


    #*#############################
    # create C matrices in python #
    #*#############################
    print("\nC_1 matrix is:")
    C_1 = mypy.create_FFT_CMat(1)
    mypy.print_naive(C_1)

    print("\nC_2 matrix is:")
    C_2 = mypy.create_FFT_CMat(2)
    mypy.print_naive(C_2)

    print("\nC_3 matrix is:")
    C_3 = mypy.create_FFT_CMat(3)
    mypy.print_naive(C_3)


    #*###############################
    # create FFT matrices in python #
    #*###############################
    print("\nF_1 matrix is:")
    F_1 = mypy.create_FFTMat(1)
    mypy.print_naive(F_1)

    print("\nF_2 matrix is:")
    F_2 = mypy.create_FFTMat(2)
    mypy.print_naive(F_2)

    print("\nF_3 matrix is:")
    F_3 = mypy.create_FFTMat(3)
    mypy.print_naive(F_3)


    #*##############################
    # apply FFT matrix to a vector #
    #*##############################
    print("\nCreate a vector\n")
    A_arr = list(map(str,[1,0,0,0,0,1,0,0]))
    rowLevel = 3
    colLevel = 0
    dimWhole = 1 << colLevel
    A_pID = mypy.row_major_list_to_store(A_arr,rowLevel,colLevel,dimWhole)
    mypy.print_naive(A_pID)

    print("Now multiply FFT matrix by the vector to get the result\n")
    B_pID = mypy.matrix_mult(F_3,A_pID)
    mypy.print_naive(B_pID)


