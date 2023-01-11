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
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mypy
from ctypes import *
import json    #* for loading parameter files
import numpy as np
import random


if __name__ == '__main__':


    #*###############################################################
    #* Commandline arguments are:                                  #*
    #*   level = (matrices are 2**level by 2**level),              #*
    #*   verbose = the local verbosity for this program, and       #*
    #*   LARC_verbose = verbosity passed to LARC initialization.   #*
    #* where: 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL          #*
    #*###############################################################

    if 7 == len(sys.argv):
        level = int(sys.argv[1]) 
        use_Cbar = int(sys.argv[2])
        use_Pbar = int(sys.argv[3])
        verbose = int(sys.argv[4])
        LARC_verbose = int(sys.argv[5])
        vector_state = int(sys.argv[6])

    else:
        print("\nThis program requires five commandline integer inputs:")
        print("   level:  matrices will be 2**level by 2**level")
        print("   use_Cbar:  1 combine C_k terms in single matrix, 0 don't")
        print("   use_Pbar:  1 combine P_k terms in single matrix, 0 don't")
        print("   verbose:  the verbosity level for this program")
        print("   LARC_verbose:  the verbosity level for the LARC package")
        print("   vector_state:  0 new random, 1 use saved, 2 new and save")
        mypy.explain_verbosity()
        print("\nSample Usage:\n")
        print(" python strategies_fft.py 3 1 0 2 0 0")
        print("  level=3  final matrix will be 8x8,")
        print("  use_Cbar=1  group the C_k terms in a single matrix,")
        print("  use_Pbar=0  use the P_k terms individually,")
        print("  level=3  means we will use 8x8 matrices,")
        print("  verbose=2  prints CHATTY comments in this python routine,")
        print("  LARC_verbose=0  prints only any errors from LARC package\n")
        print("  vector_state=0  creates and uses random 0-1 vector\n")
        print(" python strategies_fft.py 4 0 1 3 1 1")
        print("  level=4  final matrix will be 16x16,")
        print("  use_Cbar=0  use the C_k terms individually,")
        print("  use_Pbar=1  group the P_k terms in a single matrix,")
        print("  verbose=3  prints DEBUG comments in this python routine,")
        print("  LARC_verbose=1  prints warnings/errors from LARC package\n")
        print("  vector_state=1  reads previously stored random 0-1 vector\n")
        # SECRET DEVELOPER code use_Cbar=use_Pbar = 2 does everything ...
        sys.exit()

    dim = 2**level    

    if (verbose > 1):
       print("\nYou have chosen parameters:")
       print("\t%d = maximum level of matrices (size up to %d by %d)"
              %(level,dim,dim))
       print("\t%d = use_Cbar (1 to group the C_k terms, 0 not)" %use_Cbar)
       print("\t%d = use_Pbar (1 to group the P_k terms, 0 not)" %use_Pbar)
       print("\t%d = verbosity level for this routine" %verbose)
       print("\t%d = verbosity level for LARC code" %LARC_verbose)
       mypy.explain_verbosity()
        

    #*####################################################################
    # Figure out the scalarType used during compiling of MyPyLARC        #
    #*####################################################################
    # if (!mypy.cvar.scalarTypeStr):
    #        print("\nYou need to compile MyPyLARC.)
    #    explain_scalarType()
    # else:
     #   scalarType = return_scalarType()

    scalarTypeStr = mypy.cvar.scalarTypeStr

    if (verbose == 4):
    #    explain_scalarType()
        print("\nYou have compiled MyPyLARC with scalarType %s\n"%scalarTypeStr)

        
    #*##################################################
    #*   Print description sparse block algorithm     #*
    #*   from Cooley-Tukey Radix-2 Factorization      #*
    #*   see Van Loan, Computational Frameworks for   #*
    #*   the Fast Fourier Transform, p.21             #*
    #*##################################################
    if (verbose == 2):
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


    #*####################################################
    #*   Find out if machine has a large amount of      #*
    #*   memory available so we can make bigger tables  #*
    #*####################################################
    #*## MAKE CHANGES IN OTHER CODE MAYBE or WRITE SUBROUTINE
    memory_available = mypy.memory_available_GiB()
    if (verbose > 0):
        print("\nThe memory available is %ld GiB" %memory_available)
        print("We will use this to select which computing_env to read from parameter file.")
        print("You could write code to select computing_env automatically.")

    if (memory_available > 200):
        if (verbose > 0):
            print("\nMemory more than 200 GiB is considered large\n")
        computing_env = 'large'
    else:    
        if (memory_available > 50):
            if (verbose > 0):
                print("\nMemory between 50 and 200 GiB is considered medium\n")
            computing_env = 'medium'
        else:
            if (verbose > 0):
                print("\nMemory less than 50 GiB is small\n")
            computing_env = 'small'
            
    if (verbose > 0):
        print("This program believes the computing_environment is %s" %computing_env)


    #*#####################################
    #*    Print baseline usage report    #*
    #*#####################################
    if (verbose > 0):
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and hash collision chain lengths")
        mypy.memory_and_time_report(0, "stdout")


    #*##################################################################
    #*  For the FFT application in 2019 we have found that the
    #* following parameters work reasonably well, but this work should
    #* be redone now that we have better numerical accuracy.
    #*    
    #*  SMALL STORES for working on desktop and max_level <= 8:    
    #*      matrix_exponent = 22
    #*      op_exponent = 19   
    #*    
    #*  LARGE STORES for working a larger memory machine for max_level > 8:    
    #*      matrix_exponent = 31
    #*      op_exponent = 30
    #*
    #*  To get F_3 to work we have a cutoff point around -z 53 or -z 54
    #*  so in order to succeed we set:
    #*      zeroregionbitparam = 54
    #*  The rounding function does not effect whether it
    #*  works in ranges -s 10 to -s 1000 so we set
    #*      regionbitparam = 1000
    #*
    #*##################################################################

    #*##################################################################
    #*    LARC  Initialization of Matrix Store and Operation Stores   #*
    #*##################################################################
    #* The routine initialize_larc() does the following:              #*
    #* * creates the matrix and op stores                             #*
    #* * preloads matrix store with: standard scalars and gates,      #*
    #*   and with all zero, identity, and (integer) Hadamard matrices #*
    #*   left to max matrix size                                      #*
    #* To make life easier we can read the parameters for the         #*
    #* the initialization from a prestored file.                      #*
    #*##################################################################

    #* read the parameter file into a python dictionary
    with open('../../InitParams/fft.init_params','r') as init_file:
        init_param = json.load(init_file)
        for p in init_param[computing_env]:
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

    #* warn if the commandline value for LARC_verbose differs from the parameter file value for p_verbose
    #*# MAKE OTHER CODE LIKE THIS
    if (LARC_verbose != p_verbose):
        print("WARNING: Using commandline (LARC_verbose = %d) rather than param file (p_verbose = %d)!" %(LARC_verbose,p_verbose))
        print("         The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    #* initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,LARC_verbose)

    #* if in CHATTY OR DEBUG mode start a reporting thread              
    if (verbose > 1):              
        mypy.create_report_thread(report_interval_seconds)

    #* Finished with initializing LARC
    if (verbose > 1):             
        if computing_env == 'desktop':
            print("We think we are running on a desktop")
        else:  #* Large memory for cs1, cs4, cs9
            print("We think we are running on a cycle server.")

        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("stopHogging check to see if program is to large to occur once every 10 minutes.\n")


    #*##########################################################
    #* if matrices are too large do not allow naive printing  #*            
    #*##########################################################
    if (level < 4):
        print_naive = 1
    else:
        print_naive = 0
    if (verbose > 1):             
        if print_naive:
            print("  The level is small enough that we can print files of naive matrices to the screen.")
        else: 
            print("  The level= %d, is too big to reasonable print naive formated matrices to the screen." %level)

    #*##############################################
    #* print and log the parameters
    #*##############################################
    print("MatExp: %d  OpExp: %d  max_level: % d  regionbitparm: %d"
          %(matrix_exponent,op_exponent,max_level,regionbitparam))


    #*##############################################
    #* print and log the size of the matrix store #*
    #*##############################################
    num_matrices = mypy.num_matrices_in_store()   
    print("After initialization there are %d matrices in the the store"
          %num_matrices)

    path_matrices = "MatricesDFT/Level"+str(level)
            
    vlen = level + 1
            
    #*##########################
    #* preload roots of unity #*
    #*##########################
    n = 2**vlen
    root_array = [0]*(n)
    for k in range(n):
        root_array[k]  = mypy.k_th_power_of_n_th_root_of_unity_pID(
                             k,n,LARC_verbose)

    #*##############################################
    #* print and log the size of the matrix store #*
    #*##############################################
    num_matrices = mypy.num_matrices_in_store()   
    print("After loading roots of unity there are %d matrices in the the store"
          %num_matrices)
    print("MatrixID of root_array[%d] is %d" %(n-1,root_array[n-1]))

        

    # #*############################################
    # #*  retrieve matrixIDs of identity matrices #*
    # #*############################################
    # Iden = [0]*vlen
    # for i in range(vlen):
    #     Iden[i] = mypy.get_identity_pID(i)
    #     print("The matrix ID of the identity of level %d is %d"
    #           %(i,Iden[i]))

        
    # #*##############################################
    # #* write out the Ism, inverse shuffle matrices #*
    # #*##############################################
    # Ism = [0]*vlen
    # for i in range(vlen):
    #     Ism[i] = mypy.create_invShufMat(i)
    #     # mypy.fprint_larcMatrixFile(Ism[i],path_matrices+"/ISM_"+str(i))
    # Ism_size = mypy.fprint_larcMatrixFile(
    #            Ism[level],path_matrices+"/ISM_"+str(level))
    # print("LARCsize of the largest inverse shuffle matrix is %d\n" %Ism_size)

    #*###########################################################
    #* If use_Pbar is not 1 then load all the individual       #*
    #* expand_Ism files.                                       #*
    #*###########################################################
    # expand_Ism = [0]*vlen
    # # expand_Ism[1] = Iden[level]
    # for i in range(2,vlen):
    #     expand_Ism[i]=mypy.kronecker_product(Iden[level-i],Ism[i])
    #     expand_Ism_size = mypy.fprint_larcMatrixFile(
    #            expand_Ism[i],path_matrices+"/expandISM_"+str(i))
    #     print("expand_Ism[%d] has LARCsize %d, matrixID %d"
    #       %(i,expand_Ism_size,expand_Ism[i]))
    # print("")

    # #*###############################################
    # #* if use_Pbar = 1  then create the product    #*
    # #*                  of the expand_Ism matrices #*
    # #*###############################################
    # if (use_Pbar == 1):
    #     prod_Ism = Iden[level]
    #     for i in range(2,vlen):
    #         prod_Ism = mypy.matrix_mult(prod_Ism,expand_Ism[i])
    #         prod_Ism_size = mypy.fprint_larcMatrixFile(
    #            prod_Ism,path_matrices+"/prodISM_"+str(i))
    #         print("prod_Ism #%d has LARCsize %d, matrixID is %d"
    #           %(i,prod_Ism_size,prod_Ism))
    #     print("")

    if (use_Pbar == 1): 
        prodISM_ID = mypy.read_larcMatrixFile(
            path_matrices+"/prodISM_"+str(level))
        if (verbose > 1):
            print("P_bar is set and the P_bar matrix is:")
            mypy.print_naive(prodISM_ID)
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading combined ISM the matrix store has size %d"
              %num_matrices)
    elif (use_Cbar == 1):
        expandISM_ID = [0]*vlen
        for i in range(2,vlen):
            expandISM_ID[i] = mypy.read_larcMatrixFile(
                            path_matrices+"/expandISM_"+str(i))
        if (verbose > 1):
            print("P_bar is 0 and the array of matrixIDs are:")
            print([expandISM_ID[i] for i in range(2,vlen)])
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading expanded ISMs the matrix store has size %d"
              %num_matrices)

        
    if ((use_Pbar == 2) and (use_Cbar == 2)): 
        prodISM_ID = mypy.read_larcMatrixFile(
            path_matrices+"/prodISM_"+str(level))
        if (verbose > 1):
            print("P_bar is set and the P_bar matrix is:")
            mypy.print_naive(prodISM_ID)
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading combined ISM the matrix store has size %d"
              %num_matrices)
        expandISM_ID = [0]*vlen
        for i in range(2,vlen):
            expandISM_ID[i] = mypy.read_larcMatrixFile(
                            path_matrices+"/expandISM_"+str(i))
        if (verbose > 1):
            print("P_bar is 0 and the array of matrixIDs are:")
            print([expandISM_ID[i] for i in range(2,vlen)])
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading expanded ISMs the matrix store has size %d"
              %num_matrices)

        

        
    # #*###############################
    # #* create C matrices in python #*
    # #*###############################
    # Cmat = [0]*vlen
    # for i in range(1,vlen):
    #     Cmat[i] = mypy.create_FFT_CMat(i)
    #     # mypy.fprint_larcMatrixFile(Cmat[i],path_matrices+"/Cmat_"+str(i))
    #     Cmat_size = mypy.fprint_larcMatrixFile(
    #         Cmat[i],path_matrices+"/CMAT_"+str(i))
    #     print("LARCsize of the C matrix of level %d is %d, matrixID is %d"
    #           %(i,Cmat_size,Cmat[i]))
    # print("\n")


    # #*###########################################################
    # #* Create tensor products of identities and Cmat matrices     #*
    # #*   that are full size matrices 2^level by 2^level        #*
    # #*   regardless of whether use_Cbar is 0 or 1              #*
    # #*###########################################################
    # expand_C = [0]*(vlen-1)
    # for i in range(vlen-1):
    #     expand_C[i]=mypy.kronecker_product(Iden[i],Cmat[level-i])
    #     expand_C_size = mypy.fprint_larcMatrixFile(
    #         expand_C[i],path_matrices+"/expandC_"+str(i))
    #     print("LARCsize of %i-th expand_C matrix is %d, matrixID is %d"
    #           %(i,expand_C_size,expand_C[i]))
    # print("\n")

    
    # #*###############################################
    # #* if use_Cbar = 1  then create the product    #*
    # #*                  of the expand_C matrices #*
    # #*###############################################
    # if (use_Cbar == 1):
    #     prod_C = Iden[level]
    #     for i in range(level):
    #         print("input to multiplication are matrices: %d %d"
    #               %(prod_C,expand_C[i]))
    #         prod_C = mypy.matrix_mult(prod_C,expand_C[i])
    #         prod_C_size = mypy.fprint_larcMatrixFile(
    #             prod_C,path_matrices+"/prodC_partial"+str(i))
    #         print("LARCsize of the partial prod_C matrix is %d, has matID %d"
    #               %(prod_C_size,prod_C))
    #     print("\n")

    if (use_Cbar == 1):
        prodC_ID = mypy.read_larcMatrixFile(
            path_matrices+"/prodC_"+str(level))
        if (verbose > 1):
            print("C_bar is set and the C_bar matrix is:")
            mypy.print_naive(prodC_ID)
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading combined C matrix the matrix store has size %d"
              %num_matrices)
    elif (use_Pbar == 1):
        expandC_ID = [0]*(vlen-1)
        for i in range(vlen-1):
            expandC_ID[i] = mypy.read_larcMatrixFile(
                            path_matrices+"/expandC_"+str(i))
        if (verbose > 1):
            print("C_bar is 0 and the array of matrixIDs are:")
            print(expandC_ID)
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading expanded C matrices the matrix store has size %d"
              %num_matrices)
        
    if ((use_Cbar == 2) and (use_Pbar == 2)):
        prodC_ID = mypy.read_larcMatrixFile(
            path_matrices+"/prodC_"+str(level))
        if (verbose > 1):
            print("C_bar is set and the C_bar matrix is:")
            mypy.print_naive(prodC_ID)
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading combined C matrix the matrix store has size %d"
              %num_matrices)
        expandC_ID = [0]*(vlen-1)
        for i in range(vlen-1):
            expandC_ID[i] = mypy.read_larcMatrixFile(
                            path_matrices+"/expandC_"+str(i))
        if (verbose > 1):
            print("C_bar is 0 and the array of matrixIDs are:")
            print(expandC_ID)
        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading expanded C matrices the matrix store has size %d"
              %num_matrices)
        
        

    # #*###############################################
    # #* Create a DFT matrix by making a loop over   #*
    # #* expand_C and expand_Ism matrices            #*
    # #*###############################################
    # DFT_stages = Iden[level]
    # for i in range(level):
    #     print("input to multiplication are matrices: %d %d"
    #           %(DFT_stages,expand_C[i]))
    #     DFT_stages = mypy.matrix_mult(DFT_stages,expand_C[i])
    #     DFT_stages_size = mypy.fprint_larcMatrixFile(
    #         DFT_stages,path_matrices+"/DFT_stages_Cportion"+str(i))
    #     print("partial DFT matrix C_stage %d has LARCsize %d, and matID %d"
    #           %(i,DFT_stages_size,DFT_stages))
    # print("\n")

    # for i in range(2,vlen):
    #     DFT_stages =mypy.matrix_mult(DFT_stages,expand_Ism[i])
    #     DFT_stages_size = mypy.fprint_larcMatrixFile(
    #         DFT_stages,path_matrices+"/DFT_stages_IsmPortion_"+str(i))
    #     print("partial DFT matrix Ism_stage %d has LARCsize %d, and matID %d"
    #           %(i,DFT_stages_size,DFT_stages))
    # print("")


    if ((use_Pbar == 0) and (use_Cbar == 0)) or ((use_Pbar == 2) and (use_Cbar == 2)):
        FFT_ID = mypy.read_larcMatrixFile(
            path_matrices+"/DFT_stages_IsmPortion_"+str(level))
        if (verbose > 1):
            print("Using the FFT matrix:")
            mypy.print_naive(FFT_ID)

        #*##############################################
        #* print and log the size of the matrix store #*
        #*##############################################
        num_matrices = mypy.num_matrices_in_store()   
        print("After loading FFT matrices the matrix store has size %d"
              %num_matrices)

        
  
    #*##########################################
    # create random vector of the full level   #
    #*##########################################

    dim_whole = 2**level  

    # IF YOU ARE USING ROW VECTORS THIS IS OLD CODE TO HELP
    # creating or finding the matrix associated with the array
    # row vector
    # orig_num_cols = dim_whole
    #v_ID = mypy.row_major_list_to_store(v_arr,0,level,orig_num_cols)
    # print("\nThe row vector has matrix ID: %d" %v_ID)

    if (vector_state == 1):
        # read the vector
        v_ID = mypy.read_larcMatrixFile(
            path_matrices+"/fixedRandomVector_"+str(i))
        print("Read previously stored random vector from file")
        mypy.print_naive(v_ID)
    else:
        # generate random 0-1 column vector
        randVals = [np.random.randint(0,2) for i in range(dim_whole)]
        # A_str = [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]
        v_arr = mypy.map_to_str(randVals,scalarTypeStr)
        orig_num_cols = 1
        v_ID = mypy.row_major_list_to_store(v_arr, level, 0,
                                                     orig_num_cols)
        print("\nThe column vector has matrix ID: %d" %v_ID)
        print("\nOur 0-1 random vector is:")
        print(randVals)

    if (vector_state == 2):
        LARCsize_vector = mypy.fprint_larcMatrixFile(
             v_ID,path_matrices+"/fixedRandomVector_"+str(level))
        print("Wrote random vector with LARCsize %d to file" %LARCsize_vector)
    

    #*##############################################
    #* print and log the size of the matrix store #*
    #*##############################################
    num_matrices = mypy.num_matrices_in_store()   
    print("After generating a random vector the matrix store has size %d"
          %num_matrices)


#    print("multiplying the FFT matrix by the random matrix")
#    u_ID = mypy.matrix_mult(F_ID,v_ID)
#    mypy.print_naive(u_ID)

    if ((use_Pbar == 0) and (use_Cbar == 0)) or ((use_Pbar == 2) and (use_Cbar == 2)):
        result_ID = mypy.matrix_mult(FFT_ID,v_ID)
        print("The matrixID when Pbar,Cbar 0 0 is %d" %result_ID)

    if ((use_Pbar == 1) and (use_Cbar == 1)) or ((use_Pbar == 2) and (use_Cbar == 2)):
        temp_ID = mypy.matrix_mult(prodISM_ID,v_ID)
        result_ID = mypy.matrix_mult(prodC_ID,temp_ID)
        print("The matrixID when Pbar,Cbar 1 1 is %d" %result_ID)

    if ((use_Pbar == 1) and (use_Cbar == 0)) or ((use_Pbar == 2) and (use_Cbar == 2)):
        temp_ID = mypy.matrix_mult(prodISM_ID,v_ID)
        for i in reversed(range(vlen-1)):  
            temp_ID = mypy.matrix_mult(expandC_ID[i],temp_ID)
        result_ID = temp_ID
        print("The matrixID when Pbar,Cbar 1 0 is %d" %result_ID)

    if ((use_Pbar == 0) and (use_Cbar == 1)) or ((use_Pbar == 2) and (use_Cbar == 2)):
        temp_ID = v_ID
        for i in reversed(range(2,vlen)):
            temp_ID = mypy.matrix_mult(expandISM_ID[i],temp_ID)
        result_ID = mypy.matrix_mult(prodC_ID,temp_ID)
        print("The matrixID when Pbar,Cbar 0 1 is %d" %result_ID)

        
    #*##############################################
    #* print and log the size of the matrix store #*
    #*##############################################
    num_matrices = mypy.num_matrices_in_store()   
    print("After taking FFT of the vector the matrix store has size %d"
          %num_matrices)
        

    print("You have reached the exit")
    #*  YOU ARE HERE
    sys.exit()
    
