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
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mypy
from ctypes import *
import json    #* for loading parameter files
import numpy as np
import random
from pathlib import Path


if __name__ == '__main__':


    #*###############################################################
    #* Commandline arguments are:                                  #*
    #*   level = (matrices are 2**level by 2**level),              #*
    #*   verbose = the local verbosity for this program, and       #*
    #*   LARC_verbose = verbosity passed to LARC initialization.   #*
    #* where: 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL          #*
    #*###############################################################

    if 6 == len(sys.argv):
        level = int(sys.argv[1]) 
        use_Cbar = int(sys.argv[2])
        use_Pbar = int(sys.argv[3])
        verbose = int(sys.argv[4])
        LARC_verbose = int(sys.argv[5])

    else:
        print("\nThis program requires five commandline integer inputs:")
        print("   level:  matrices will be 2**level by 2**level")
        print("   use_Cbar:  1 combine C_k terms in single matrix, 0 don't")
        print("   use_Pbar:  1 combine P_k terms in single matrix, 0 don't")
        print("   verbose:  the verbosity level for this program")
        print("   LARC_verbose:  the verbosity level for the LARC package")
        mypy.explain_verbosity()
        print("\nSample Usage:\n")
        print(" python store_fft_matrices.py 3 1 0 2 0")
        print("  level=3  final matrix will be 8x8,")
        print("  use_Cbar=1  group the C_k terms in a single matrix,")
        print("  use_Pbar=0  use the P_k terms individually,")
        print("  verbose=2  prints CHATTY comments in this python routine,")
        print("  LARC_verbose=0  prints only any errors from LARC package\n")
        print(" python store_fft_matrices.py 4 0 1 3 1")
        print("  level=4  final matrix will be 16x16,")
        print("  use_Cbar=0  use the C_k terms individually,")
        print("  use_Pbar=1  group the P_k terms in a single matrix,")
        print("  verbose=3  prints DEBUG comments in this python routine,")
        print("  LARC_verbose=1  prints warnings/errors from LARC package\n")
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
    #*  For the FFT application we have found that the following parameters
    #*  work reasonably well:
    #*    
    #*  SMALL STORES for working on desktop and max_level <= 8:    
    #*      matrix_exponent = 22
    #*      op_exponent = 19   
    #*    
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
        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")


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


    vlen = level + 1
            
    #*##########################
    #* preload roots of unity #*
    #*##########################
    n = 2**vlen
    root_array = [0]*(n)
    for k in range(n):
        root_array[k]  = mypy.k_th_power_of_n_th_root_of_unity_pID(
                             k,n,LARC_verbose)


    #*############################################
    #*  retrieve matrixIDs of identity matrices #*
    #*############################################
    Iden = [0]*vlen
    for i in range(vlen):
        Iden[i] = mypy.get_identity_pID(i)
        print("The matrix ID of the identity of level %d is %d"
              %(i,Iden[i]))

    #*################################################
    #*  Find name of subdirectory for storing files #*
    #*  and if it doesn't exist create it           #*
    #*################################################
    path_matrices = "MatricesDFT/Level"+str(level)
    Path(path_matrices).mkdir(parents=True,exist_ok=True)
        
    #*##############################################
    #* write out the Ism, inverse shuffle matrices #*
    #*##############################################
    Ism = [0]*vlen
    for i in range(vlen):
        Ism[i] = mypy.create_invShufMat(i)
    # Ism_size = mypy.fprint_larcMatrixFile(
    #           Ism[level],path_matrices+"/ISM_"+str(level))
    # print("LARCsize of the largest inverse shuffle matrix is %d\n" %Ism_size)

    #*###########################################################
    #* Create tensor products of identities and Ism matrices   #*
    #*   that are full size matrices 2^level by 2^level        #*
    #*   regardless of whether use_Pbar is 0 or 1              #*
    #*###########################################################
    expand_Ism = [0]*vlen
    # expand_Ism[1] = Iden[level]
    for i in range(2,vlen):
        expand_Ism[i]=mypy.kronecker_product(Iden[level-i],Ism[i])
        expand_Ism_size = mypy.fprint_larcMatrixFile(
               expand_Ism[i],path_matrices+"/expandISM_"+str(i))
        expand_Ism_size = mypy.fprint_larcMatrixFile(
               expand_Ism[i],path_matrices+"/expandISM_"+str(i))
        print("expand_Ism[%d] has LARCsize %d, matrixID %d"
          %(i,expand_Ism_size,expand_Ism[i]))
    print("");

    #*###############################################
    #* if use_Pbar = 1  then create the product    #*
    #*                  of the expand_Ism matrices #*
    #*###############################################
    if (use_Pbar == 1):
        prod_Ism = Iden[level]
        for i in range(2,vlen):
            prod_Ism = mypy.matrix_mult(prod_Ism,expand_Ism[i])
        prod_Ism_size = mypy.fprint_larcMatrixFile(
            prod_Ism,path_matrices+"/prodISM_"+str(level))
        print("prod_Ism_%d has LARCsize %d, matrixID is %d"
              %(level,prod_Ism_size,prod_Ism))
        print("");
        

    #*###############################
    #* create C matrices in python #*
    #*###############################
    Cmat = [0]*vlen
    for i in range(1,vlen):
        Cmat[i] = mypy.create_FFT_CMat(i)
#        Cmat_size = mypy.fprint_larcMatrixFile(
#            Cmat[i],path_matrices+"/CMAT_"+str(i))
#        print("LARCsize of the C matrix of level %d is %d, matrixID is %d"
#              %(i,Cmat_size,Cmat[i]))
#    print("\n")


    #*###########################################################
    #* Create tensor products of identities and Cmat matrices     #*
    #*   that are full size matrices 2^level by 2^level        #*
    #*   regardless of whether use_Cbar is 0 or 1              #*
    #*###########################################################
    expand_C = [0]*(vlen-1)
    for i in range(vlen-1):
        expand_C[i]=mypy.kronecker_product(Iden[i],Cmat[level-i])
        expand_C_size = mypy.fprint_larcMatrixFile(
            expand_C[i],path_matrices+"/expandC_"+str(i))
        print("LARCsize of %i-th expand_C matrix is %d, matrixID is %d"
              %(i,expand_C_size,expand_C[i]))
    print("\n")

    
    #*###############################################
    #* if use_Cbar = 1  then create the product    #*
    #*                  of the expand_C matrices #*
    #*###############################################
    if (use_Cbar == 1):
        prod_C = Iden[level]
        for i in range(level):
            print("input to multiplication are matrices: %d %d"
                  %(prod_C,expand_C[i]))
            prod_C = mypy.matrix_mult(prod_C,expand_C[i])
        prod_C_size = mypy.fprint_larcMatrixFile(
            prod_C,path_matrices+"/prodC_"+str(level))
        print("LARCsize of the prod_C_%d matrix is %d, has matID %d"
            %(level,prod_C_size,prod_C))
        print("\n")


    #*###############################################
    #* Create a DFT matrix by making a loop over   #*
    #* expand_C and expand_Ism matrices            #*
    #*###############################################
    DFT_stages = Iden[level]
    for i in range(level):
        print("input to multiplication are matrices: %d %d"
              %(DFT_stages,expand_C[i]))
        DFT_stages = mypy.matrix_mult(DFT_stages,expand_C[i])
        DFT_stages_size = mypy.fprint_larcMatrixFile(
            DFT_stages,path_matrices+"/DFT_stages_Cportion"+str(i))
        print("partial DFT matrix C_stage %d has LARCsize %d, and matID %d"
              %(i,DFT_stages_size,DFT_stages))
    print("\n")

    for i in range(2,vlen):
        DFT_stages =mypy.matrix_mult(DFT_stages,expand_Ism[i])
        DFT_stages_size = mypy.fprint_larcMatrixFile(
            DFT_stages,path_matrices+"/DFT_stages_IsmPortion_"+str(i))
        print("partial DFT matrix Ism_stage %d has LARCsize %d, and matID %d"
              %(i,DFT_stages_size,DFT_stages))
    print("");

  

# YOU ARE HERE
            
    #*#############################
    # write out Ism matrices  #
    # calculate the product 
    #*#############################
    # Ism_prod = [0]*vlen
    # Ism_prod[0] = Iden[0]
    # # mypy.fprint_larcMatrixFile(Ism_prod[0],path_matrices+"/ISM_prod0")
    # Ism_prod[1] = Ism[1]
    # # mypy.fprint_larcMatrixFile(Ism_prod[1],path_matrices+"/ISM_prod1")
    # for i in range(2,vlen):
    #     Ism_prod[i] = mypy.matrix_mult(
    #         mypy.kronecker_product(Iden[1],Ism_prod[i-1]), Ism[i])
    #     # mypy.fprint_larcMatrixFile(Ism_prod[i],path_matrices+"/ISM_prod"+str(i))

    #*#############################
    # write out lvl8 c matrices   #
    #*#############################
    # c_prod = [0]*vlen
    # c_prod[1] = c[1]
    # mypy.fprint_larcMatrixFile(c_prod[1],path_matrices+"/C_prod1")
    # for i in range(2,vlen):
    #     c_prod[i] = mypy.matrix_mult( c[i],
    #         mypy.kronecker_product(Iden[1],c_prod[i-1]) )
    #     mypy.fprint_larcMatrixFile(c_prod[i],path_matrices+"/C_prod"+str(i))

    #*#############################
    # write out fft matrices      #
    #*#############################
    # f = [0]*vlen
    # for i in range(1,vlen):
    #     f[i] = mypy.matrix_mult(c_prod[i],Ism_prod[i])
    #     mypy.fprint_larcMatrixFile(f[i],path_matrices+"/f"+str(i))


    #*##########################################
    # create random vector of the full level   #
    #*##########################################
    dim_whole = 2**level  
    randVals = [np.random.randint(0,2) for i in range(dim_whole)]
    # A_str = [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]
    v_arr = mypy.map_to_str(randVals,scalarTypeStr)

    # creating or finding the matrix associated with the array
    # row vector
    # orig_num_cols = dim_whole
    #v_ID = mypy.row_major_list_to_store(v_arr,0,level,orig_num_cols)

    # col vector
    orig_num_cols = 1
    v_ID = mypy.row_major_list_to_store(v_arr, level, 0, orig_num_cols)
    print("\nThe vector has matrix ID: %d" %v_ID)

    print("\nOur 0-1 random vector is:")
    print(randVals)

#    print("multiplying the FFT matrix by the random matrix")
#    u_ID = mypy.matrix_mult(F_ID,v_ID)
#    mypy.print_naive(u_ID)  


    print("You have reached the exit")
    # YOU ARE HERE
    sys.exit()
    
    # #*##############################
    # # inverse permutation matrices #
    # #*##############################
    # if (verbose == 3): # DEBUG mode (verbose = 3)
    #     print("\nISM_0 matrix is:")
    #     ISM_0 = mypy.create_invShufMat(0)
    #     if print_naive:
    #         mypy.print_naive(ISM_0)
            
    #     print("\nISM_1 matrix is:")
    #     ISM_1 = mypy.create_invShufMat(1)
    #     if print_naive:
    #         mypy.print_naive(ISM_1)
            
    #     print("\nISM_2 matrix is:")
    #     ISM_2 = mypy.create_invShufMat(2)
    #     if print_naive:
    #         mypy.print_naive(ISM_2)
                    
    #     print("\nISM_3 matrix is:")
    #     ISM_3 = mypy.create_invShufMat(3)
    #     if print_naive:
    #         mypy.print_naive(ISM_3)

    # ISM_name = "ISM_"+str(level)        
    # ISM_ID = mypy.create_invShufMat(level)
    # print("\nFormed %s the inverse shuffle matrix for a DFT: it has matrixID %d." %(ISM_name,ISM_ID))
    # filename = "Data/Out/%s_%s.json" %(ISM_name,scalarTypeStr) 
    # if (print_naive and (verbose > 1)):
    #     mypy.print_naive(ISM_ID)
    # if (verbose > 1):
    #     print("Printed %s into LARC compressed formated json file %s." %(ISM_name,filename))
    # ISM_larcSize = mypy.fprint_larcMatrixFile(ISM_ID,filename)
    # print("The %s matrix of level %d has larcSize %d" %(ISM_name,level,ISM_larcSize));


    # #*#######################
    # # print roots of unity  #
    # #*#######################
    # if (verbose == 3): # DEBUG mode (verbose = 3)
    #     print("\nTesting calculation for roots of unity:")
    #     mypy.print_pow2_roots_unity(1)
    #     mypy.print_pow2_roots_unity(2)
    #     mypy.print_pow2_roots_unity(3)


    # #*#############################
    # # create D matrices in python #
    # #*#############################
    # if (verbose == 3): # DEBUG mode (verbose = 3)
    #     print("\nD_0 matrix is:")
    #     D_0 = mypy.create_FFT_DMat(0)
    #     if print_naive:
    #         mypy.print_naive(D_0)

    #     print("\nD_1 matrix is:")
    #     D_1 = mypy.create_FFT_DMat(1)
    #     if print_naive:
    #         mypy.print_naive(D_1)

    #     print("\nD_2 matrix is:")
    #     D_2 = mypy.create_FFT_DMat(2)
    #     if print_naive:
    #         mypy.print_naive(D_2)

    #     print("\nD_3 matrix is:")
    #     D_3 = mypy.create_FFT_DMat(3)
    #     if print_naive:
    #         mypy.print_naive(D_3)

    # D_name = "D_"+str(level)        
    # D_ID = mypy.create_FFT_DMat(level)
    # print("\nFormed %s the D matrix for an DFT; it has matrixID %d." %(D_name,D_ID))
    # filename = "Data/Out/%s_%s.json" %(D_name,scalarTypeStr) 
    # if (print_naive and (verbose > 1)):
    #     mypy.print_naive(D_ID)
    # if (verbose > 1):
    #     print("Printed %s into LARC compressed formated json file %s." %(D_name,filename))
    # D_larcSize = mypy.fprint_larcMatrixFile(D_ID,filename)
    # print("The %s matrix of level %d has larcSize %d" %(D_name,level,D_larcSize));


    # #*#############################
    # # create C matrices in python #
    # #*#############################
    # if (verbose == 3): # DEBUG mode (verbose = 3)
    #     print("\nC_1 matrix is:")
    #     C_1 = mypy.create_FFT_CMat(1)
    #     if print_naive:
    #         mypy.print_naive(C_1)

    #     print("\nC_2 matrix is:")
    #     C_2 = mypy.create_FFT_CMat(2)
    #     if print_naive:
    #         mypy.print_naive(C_2)

    #     print("\nC_3 matrix is:")
    #     C_3 = mypy.create_FFT_CMat(3)
    #     if print_naive:
    #         mypy.print_naive(C_3)

    # C_name = "C_"+str(level)        
    # C_ID = mypy.create_FFT_CMat(level)
    # print("\nFormed %s the C matrix for a DFT; it has matrixID %d." %(C_name,C_ID))
    # filename = "Data/Out/%s_%s.json" %(C_name,scalarTypeStr) 
    # if (print_naive and (verbose > 1)):
    #     mypy.print_naive(C_ID)
    # if (verbose > 1):
    #     print("Printed %s into LARC compressed formated json file %s." %(C_name,filename))
    # C_larcSize = mypy.fprint_larcMatrixFile(C_ID,filename)
    # print("The %s matrix of level %d has larcSize %d" %(C_name,level,C_larcSize));

        
    # #*###############################
    # # create FFT matrices in python #
    # #*###############################
    # if (verbose == 3): # DEBUG mode (verbose = 3)
    #     print("\nF_1 matrix is:")
    #     F_1 = mypy.create_FFTMat(1)
    #     if print_naive:
    #         mypy.print_naive(F_1)

    #     print("\nF_2 matrix is:")
    #     F_2 = mypy.create_FFTMat(2)
    #     if print_naive:
    #         mypy.print_naive(F_2)

    #     print("\nF_3 matrix is:")
    #     F_3 = mypy.create_FFTMat(3)
    #     if print_naive:
    #         mypy.print_naive(F_3)

    # F_name = "F_"+str(level)        
    # F_ID = mypy.create_FFTMat(level)
    # print("\nFormed %s the DFT matrix; it has matrixID %d." %(F_name,F_ID))
    # filename = "Data/Out/%s_%s.json" %(F_name,scalarTypeStr) 
    # if (print_naive and (verbose > 1)):
    #     mypy.print_naive(F_ID)  
    
    # if (verbose > 1):
    #     print("Printed %s into LARC compressed formated json file %s." %(F_name,filename))
    # F_larcSize = mypy.fprint_larcMatrixFile(F_ID,filename)
    # print("The %s matrix of level %d has larcSize %d" %(F_name,level,F_larcSize));

    # #*##########################################
    # # create random vector of the full level   #
    # #*##########################################
    # dim_whole = 2**level  
    # randVals = [np.random.randint(0,2) for i in range(dim_whole)]
    # # A_str = [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]
    # v_arr = mypy.map_to_str(randVals,scalarTypeStr)

    # # creating or finding the matrix associated with the array
    # # row vector
    # # orig_num_cols = dim_whole
    # #v_ID = mypy.row_major_list_to_store(v_arr,0,level,orig_num_cols)

    # # col vector
    # orig_num_cols = 1
    # v_ID = mypy.row_major_list_to_store(v_arr, level, 0, orig_num_cols)
    # print("\nThe vector has matrix ID: %d" %v_ID)

    # print("\nOur 0-1 random vector is:")
    # print(randVals)

    # print("multiplying the FFT matrix by the random matrix")
    # u_ID = mypy.matrix_mult(F_ID,v_ID)
    # mypy.print_naive(u_ID)  

    # sys.exit(0)

    
    # # #*###################################################
    # # create random Toeplitz matrices with appropriate  #
    # #*###################################################
    # dim_whole = 2**level

    # randVals = [np.random.randint(0,1) for i in range(2*dim_whole-1)]

    # if scalarTypeStr in ('Real', 'MPReal', 'MPRational'):
    #     randVals = [ random.random() for i in range(2*dim_whole-1)]
    # elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
    #     randVals = [ np.complex(random.random(),random.random()) for i in range(2*dim_whole-1)]
    # elif scalarTypeStr in ('Integer', 'MPInteger'):
    #     randVals = [ random.randrange(0,10001) for i in range(2*dim_whole-1)]
    # else:
    #     raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # # create toeplitz matrix from the 2*dim_whole-1 random numbers
    # a = []
    # b = randVals
    # for i in range(dim_whole):
    #     a.append(list(b[dim_whole-i-1:2*dim_whole-i-1]))
    # amat = np.matrix(a)
    # alist = amat.reshape(-1).tolist()[0]
    # arr = mypy.map_to_str(alist, scalarTypeStr)
    # if ((verbose == 3) and print_naive):
    #     print("\nCreating random Toeplitz matrix")
    #     print("\n  The random values are: ", b)
    #     print("\n  The list that will be reshaped for the matrix is: ", a)
    #     print("\n  Which gives the shaped list of values:")
    #     print(alist)
    #     print("\n  And the string array:")
    #     print(mypy.str_scalarTypeArray(arr, len(alist)))

    # T_name = "T_"+str(level)
    # # creating or finding the matrix associated with the array
    # T_ID = mypy.row_major_list_to_store(arr, level, level, dim_whole)
    # filename = "Data/Out/%s_%s.json" %(T_name,scalarTypeStr)
    # if (print_naive and (verbose > 1)):
    #     mypy.print_naive(T_ID)
    # print("\nPrinted the Toeplitz matrix %s to file %s." %(T_name,filename))
    # T_larcSize = mypy.fprint_larcMatrixFile(T_ID,filename)
    # print("The %s matrix of level %d has larcSize %d" %(T_name,level,T_larcSize));

    
    # #*#########################################################
    # # apply inverse perm matrix PI_l to a Toeplitz matrix T_l #
    # # where l is the level.                                   #
    # #*#########################################################
    # print("\nNow multiplying the inverse shuffle matrix %s by the Toeplitz matrix %s" %(PI_name,T_name))
    # PIT_name = "PIT_"+str(level)
    # PIT_ID = mypy.matrix_mult(PI_ID,T_ID)
    # if (print_naive and (verbose > 1)):
    #     mypy.print_naive(PIT_ID)
    # filename = "Data/Out/%s_%s.json" %(PIT_name,scalarTypeStr)
    # PIT_larcSize = mypy.fprint_larcMatrixFile(PIT_ID,filename)
    # print("The %s matrix of level %d has larcSize %d" %(PIT_name,level,PIT_larcSize));

