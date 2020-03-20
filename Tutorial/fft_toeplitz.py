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
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
from ctypes import *
import json    ## for loading parameter files
import numpy as np
import random


if __name__ == '__main__':


    ######################################################################
    ## Set the level (matrices are 2**level by 2**level), the           ##
    ## LARC verbosity (0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL), ##
    ## and the LOCAL_verbosity.                                         ##
    ######################################################################

    if 3 == len(sys.argv):
        level = int(sys.argv[1]) 
        verbose = int(sys.argv[2])    
    else:
        print("\nThis program requires two commandline integer inputs: l v")
        print("   l:  level (matrices will be 2**level by 2**level)")
        print("   v:  verbosity")
        mypy.explain_verbosity()
        print("\nSample Usage:")
        print("   python fft_toeplitz 3 2\n")
        print("      would use 8x8 matrices and verbose=CHATTY.\n")
        print("   python fft_toeplitz 4 1\n")
        print("      would use 16x16 matrices and verbose=BASIC.\n")
        sys.exit()    

    dim = 2**level    

    if (verbose > 1):
       print("Program will use level=%d (%d by %d matrices)" %(level,dim,dim))
       mypy.explain_verbosity
    if (verbose == 2):
        print("The verbosity level is 2=CHATTY. Other options were: ")
        print("\t0=SILENT, 1=BASIC, 3=DEBUG, and 4=ALL.")
    if (verbose == 3):
        print("The verbosity level is 3=DEBUG. Other options were: ")
        print("\t0=SILENT, 1=BASIC, 2=CHATTY, and 4=ALL.")
    if (verbose == 4):
        print("The verbosity level is 4=ALL. Other options were :")
        print("\t0=SILENT, 1=BASIC, 2=CHATTY, and 3=DEBUG.")


    ######################################################################
    # Figure out the scalarType used during compiling of MyPyLARC        #
    ######################################################################
    # if (!mypy.cvar.scalarTypeStr):
    #        print("\nYou need to compile MyPyLARC.)
    #    explain_scalarType()
    # else:
     #   scalarType = return_scalarType()

    scalarTypeStr = mypy.cvar.scalarTypeStr

    if (verbose == 4):
    #    explain_scalarType()
        print("\nYou have compiled MyPyLARC with scalarType %s\n"%scalarTypeStr)
        

    ####################################################
    ##   Print description sparse block algorithm     ##
    ##   from Cooley-Tukey Radix-2 Factorization      ##
    ##   see Van Loan, Computational Frameworks for   ##
    ##   the Fast Fourier Transform, p.21             ##
    ####################################################
    if (verbose == 2):
        print("\nWe will create the matrices used in the Cooley-Tukey")
        print("radix-2 sparse block recursive FFT that is from Van Loan,")
        print("Computational Frameworks for the Fast Fourier Transform. p.21")
        if verbose:
            print("The level k, 2**k by 2**k Fourier matrix F_k can be")
            print("generated recursively by the equation")
            print("   F_k = C_k * (I_2 @ F_(k-1) ) * PI_k")
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



    ##   TODO: Hunting in uname for computing_env is a reference to local names cs2,cs3, ...
    ######################################################
    ##   Find out if machine has a large amount of      ##
    ##   memory available so we can make bigger tables  ##
    ######################################################
    memory_available = mypy.memory_available_GiB()
    print("\nThe memory available is %ld GiB\n" %memory_available)
    if (memory_available > 200):
        print("\n This memory is more than 200 GiB\n")
    else:
        print("\n This memory is less than 200 GiB\n")

    userInput= input("Press <Enter> to continue\n")



    
    ####################################################
    ##   Find out if machine is desktop workstation   ##
    ##   or a CPU-cycle servers (cs1-cs6)             ##
    ####################################################
    machine = os.uname()[1]
    cs = 0        # on desktop workstation, with smaller memory
    computing_env = 'desktop'
    if (machine.find('cs') >= 0):
        cs = 1    # on CPU-cycle server cs1-cs6, with larger memory
        computing_env = 'server'
        if (verbose > 1):
            print(" code on a CPU-cycle server")
    else:
        if (verbose > 1):
            print("This machine is a desktop work station")


    #######################################
    ##    Print baseline usage report    ##
    #######################################
    if (verbose > 1):
        mypy.rusage_report(0, "stdout")
    

    ####################################################################
    ##  For the FFT application we have found that the following parameters
    ##  work reasonably well:
    ##    
    ##  SMALL STORES for working on desktop and max_level <= 8:    
    ##      matrix_exponent = 22
    ##      op_exponent = 19   
    ##    
    ##    
    ##  LARGE STORES for working on cycle server for max_level > 8:    
    ##      matrix_exponent = 31
    ##      op_exponent = 30
    ##
    ##  To get F_3 to work we have a cutoff point around -z 53 or -z 54
    ##  so in order to succeed we set:
    ##      trunc_to_zero_bits = 54
    ##  The rounding function does not effect whether it
    ##  works in ranges -s 10 to -s 1000 so we set
    ##      rnd_sig_bits = 1000
    ##
    ####################################################################

    ## OLD NOTES
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

    ####################################################################
    ##    LARCt Initialization of Matrix Store and Operation Stores   ##
    ####################################################################
    ## The routine initialize_larc() does the following:              ##
    ## * creates the matrix and op stores                             ##
    ## * preloads matrix store with: standard scalars and gates,      ##
    ##   and with all zero, identity, and (integer) Hadamard matrices ##
    ##   left to max matrix size                                      ##
    ## To make life easier we can read the parameters for the         ##
    ## the initialization from a prestored file.                      ##
    ####################################################################

    ## read the parameter file into a python dictionary
    with open('../InitParams/fft.init_params','r') as init_file:
        init_param = json.load(init_file)
        for p in init_param[computing_env]:
            print('MatrixExponent: %d' %(p['matrix_exponent']))
            print('OpExponent: %d' %(p['op_exponent']))
            print('MaxLevel: %d' %(p['max_level']))
            print('RoundSigBits: %d' %(p['rnd_sig_bits']))
            print('TruncToZeroBits: %d' %(p['trunc_to_zero_bits']))
            print('ReportIntervalSecs: %d' %(p['report_interval_seconds']))
            print('Verbose: %d' %(p['verbose']))
            print('')
            matrix_exponent = p['matrix_exponent']
            op_exponent = p['op_exponent']
            max_level= p['max_level']
            rnd_sig_bits = p['rnd_sig_bits']
            trunc_to_zero_bits = p['trunc_to_zero_bits']
            report_interval_seconds = p['report_interval_seconds']
            p_verbose = p['verbose']

    ## warn if the commandline value for verbose differs from the parameter file value for verbose        
    if (verbose > 0):
        if (verbose != p_verbose):
            print("WARNING: Using commandline (verbose = %d) rather than param file (verbose = %d)!" %(verbose,p_verbose))
            print("         The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    ## initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)

    ## if in CHATTY OR DEBUG mode start a reporting thread              
    if (verbose > 1):              
        mypy.create_report_thread(report_interval_seconds)

    ## Finished with initializing LARC
    if (verbose > 1):             
        if computing_env == 'desktop':
            print("We think we are running on a desktop")
        else:  ## Large memory for cs1, cs4, cs9
            print("We think we are running on a cycle server.")

        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("Seppuku check to see if program is to large to occur once every 10 minutes.\n")


    ############################################################
    ## if matrices are too large do not allow naive printing  ##            
    ############################################################
    if (level < 4):
        print_naive = 1
    else:
        print_naive = 0
    if (verbose > 1):             
        if print_naive:
            print("  The level is small enough that we can print files of naive matrices to the screen.")
        else: 
            print("  The level= %d, is too big to reasonable print naive formated matrices to the screen." %level)
                  

    ################################
    # inverse permutation matrices #
    ################################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nPI_0 matrix is:")
        PI_0 = mypy.create_perm_inv_matrixID(0)
        if print_naive:
            mypy.print_naive_by_matID(PI_0)
            
        print("\nPI_1 matrix is:")
        PI_1 = mypy.create_perm_inv_matrixID(1)
        if print_naive:
            mypy.print_naive_by_matID(PI_1)
            
        print("\nPI_2 matrix is:")
        PI_2 = mypy.create_perm_inv_matrixID(2)
        if print_naive:
            mypy.print_naive_by_matID(PI_2)
                    
        print("\nPI_3 matrix is:")
        PI_3 = mypy.create_perm_inv_matrixID(3)
        if print_naive:
            mypy.print_naive_by_matID(PI_3)

    PI_name = "PI_"+str(level)        
    PI_ID = mypy.create_perm_inv_matrixID(level)
    print("\nFormed %s the inverse shuffle matrix for a DFT: it has matrixID %d." %(PI_name,PI_ID))
    filename = "Data/Out/%s_%s.json" %(PI_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive_by_matID(PI_ID)
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(PI_name,filename))
    PI_larcSize = mypy.write_larcMatrix_file_return_larcSize_by_matID(PI_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(PI_name,level,PI_larcSize));


    #########################
    # print roots of unity  #
    #########################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nTesting calculation for roots of unity:")
        mypy.print_pow2_roots_unity(1)
        mypy.print_pow2_roots_unity(2)
        mypy.print_pow2_roots_unity(3)


    ###############################
    # create D matrices in python #
    ###############################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nD_0 matrix is:")
        D_0 = mypy.create_fft_D_matrixID(0)
        if print_naive:
            mypy.print_naive_by_matID(D_0)

        print("\nD_1 matrix is:")
        D_1 = mypy.create_fft_D_matrixID(1)
        if print_naive:
            mypy.print_naive_by_matID(D_1)

        print("\nD_2 matrix is:")
        D_2 = mypy.create_fft_D_matrixID(2)
        if print_naive:
            mypy.print_naive_by_matID(D_2)

        print("\nD_3 matrix is:")
        D_3 = mypy.create_fft_D_matrixID(3)
        if print_naive:
            mypy.print_naive_by_matID(D_3)

    D_name = "D_"+str(level)        
    D_ID = mypy.create_fft_D_matrixID(level)
    print("\nFormed %s the D matrix for an DFT; it has matrixID %d." %(D_name,D_ID))
    filename = "Data/Out/%s_%s.json" %(D_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive_by_matID(D_ID)
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(D_name,filename))
    D_larcSize = mypy.write_larcMatrix_file_return_larcSize_by_matID(D_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(D_name,level,D_larcSize));


    ###############################
    # create C matrices in python #
    ###############################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nC_1 matrix is:")
        C_1 = mypy.create_fft_C_matrixID(1)
        if print_naive:
            mypy.print_naive_by_matID(C_1)

        print("\nC_2 matrix is:")
        C_2 = mypy.create_fft_C_matrixID(2)
        if print_naive:
            mypy.print_naive_by_matID(C_2)

        print("\nC_3 matrix is:")
        C_3 = mypy.create_fft_C_matrixID(3)
        if print_naive:
            mypy.print_naive_by_matID(C_3)

    C_name = "C_"+str(level)        
    C_ID = mypy.create_fft_C_matrixID(level)
    print("\nFormed %s the C matrix for a DFT; it has matrixID %d." %(C_name,C_ID))
    filename = "Data/Out/%s_%s.json" %(C_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive_by_matID(C_ID)
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(C_name,filename))
    C_larcSize = mypy.write_larcMatrix_file_return_larcSize_by_matID(C_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(C_name,level,C_larcSize));

        
    #################################
    # create FFT matrices in python #
    #################################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nF_1 matrix is:")
        F_1 = mypy.create_fft_matrix_matrixID(1)
        if print_naive:
            mypy.print_naive_by_matID(F_1)

        print("\nF_2 matrix is:")
        F_2 = mypy.create_fft_matrix_matrixID(2)
        if print_naive:
            mypy.print_naive_by_matID(F_2)

        print("\nF_3 matrix is:")
        F_3 = mypy.create_fft_matrix_matrixID(3)
        if print_naive:
            mypy.print_naive_by_matID(F_3)

    F_name = "F_"+str(level)        
    F_ID = mypy.create_fft_matrix_matrixID(level)
    print("\nFormed %s the DFT matrix; it has matrixID %d." %(F_name,F_ID))
    filename = "Data/Out/%s_%s.json" %(F_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive_by_matID(F_ID)  
    
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(F_name,filename))
    F_larcSize = mypy.write_larcMatrix_file_return_larcSize_by_matID(F_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(F_name,level,F_larcSize));


        
    ################################################
    # create random Toeplitz matrices with level 3 #
    ################################################
    ## level = 3
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
    for i in range(dim_whole):
        a.append(list(b[dim_whole-i-1:2*dim_whole-i-1]))
    amat = np.matrix(a)
    alist = amat.reshape(-1).tolist()[0]
    arr = mypy.map_to_str(alist, scalarTypeStr)
    if ((verbose == 3) and print_naive):
        print("\nCreating random Toeplitz matrix")
        print("\n  The random values are: ", b)
        print("\n  The list that will be reshaped for the matrix is: ", a)
        print("\n  Which gives the shaped list of values:")
        print(alist)
        print("\n  And the string array:")
        print(mypy.str_scalarTypeArray(arr, len(alist)))

    T_name = "T_"+str(level)
    # creating or finding the matrix associated with the array
    T_ID = mypy.row_major_list_to_store_matrixID(arr, level, level, dim_whole)
    filename = "Data/Out/%s_%s.json" %(T_name,scalarTypeStr)
    if (print_naive and (verbose > 1)):
        mypy.print_naive_by_matID(T_ID)
    print("\nPrinted the Toeplitz matrix %s to file %s." %(T_name,filename))
    T_larcSize = mypy.write_larcMatrix_file_return_larcSize_by_matID(T_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(T_name,level,T_larcSize));

    
    ###########################################################
    # apply inverse perm matrix PI_l to a Toeplitz matrix T_l #
    # where l is the level.                                   #
    ###########################################################
    print("\nNow multiplying the inverse shuffle matrix %s by the Toeplitz matrix %s" %(PI_name,T_name))
    PIT_name = "PIT_"+str(level)
    PIT_ID = mypy.matrix_mult_matrixID(PI_ID,T_ID)
    if (print_naive and (verbose > 1)):
        mypy.print_naive_by_matID(PIT_ID)
    filename = "Data/Out/%s_%s.json" %(PIT_name,scalarTypeStr)
    PIT_larcSize = mypy.write_larcMatrix_file_return_larcSize_by_matID(PIT_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(PIT_name,level,PIT_larcSize));

