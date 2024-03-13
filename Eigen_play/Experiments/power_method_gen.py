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
import json    # for loading parameter files
import numpy as np
import random



if __name__ == '__main__':


    #*##################################################################
    # Set the level (matrices are 2**level by 2**level                   #
    # and the verbosity (0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL) #
    #*####################################################################

    if 4 == len(sys.argv):
        level = int(sys.argv[1]) 
        verbose = int(sys.argv[2])
        LARC_verbose = int(sys.argv[3])

    else:
        print("\nThis program requires three commandline integer inputs:")
        print("   level:  matrices will be 2**level by 2**level")
        print("   verbose:  the verbosity level for this program")
        print("   LARC_verbose:  the verbosity level for the LARC package")
        mypy.explain_verbosity()
        print("\nSample Usage:")
        print("  python power_method.py 3 2 0")
        print("    level=3  means we will use 8x8 matrices,")
        print("    verbose=2  prints CHATTY comments in this python routine,")
        print("    LARC_verbose=0  prints only any errors from LARC package")
        print("  python power_method.py 4 3 1")
        print("    level=4  means we will use 16x16 matrices,")
        print("    verbose=3  prints DEBUG comments in this python routine,")
        print("    LARC_verbose=1  prints warnings/errors from LARC package\n")
        sys.exit()


    dim = 2**level
    if (verbose > 1):
       print("\nlevel=%d: Program will use %d by %d matrices." %(level,dim,dim))
    if ((verbose == 2) or (verbose == 4)):
        print("Verbosity:")
        print("\tThe local verbosity level is verbose=%d and the" %verbose)
        print("\tverbosity level for the LARC package is %d," %LARC_verbose)
        print("\twhere 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL.")


    #*####################################################
    #*  General Description of What this Program Does   ## 
    #*####################################################
    # input("\nHit return to continue")
    if (verbose > 1):
        print("\nThis code implements the power method for finding the largest")
        print("magnitude real eigenvalue of the starting matrix E.")
        print("The starting matrix E will be randomly generated from real")
        print("inputs, and the initial vector given to the power method will")
        print("be a real unit vector.")
        
    #*#############################
    #*  Generate a random matrix ##
    #*#############################
    E = np.random.normal(size=(dim,dim))
    print(E)

    #*######################
    #* LARC initalization ##
    #*######################
    userInput = input("\nEnter 'q'<RETURN> to quit, <RETURN> to initialize LARC and continue: ")
    if (userInput == 'q'):
        sys.exit(1)
    
        
    #*####################################################################
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*####################################################################
    scalarTypeStr = mypy.cvar.scalarTypeStr
#    if (scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
#        print("\nComplex scalarType chosen - you may proceed")
#    else:
#        print("\nYou have chosen a non-complex scalarType for finding")
#        print("the eigenvalues and engenvectors of a general real matrix")
#        print("This will not work.")
#        sys.exit(1)

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

    userInput= input("Press <Enter> to continue\n")
    

    #*#####################################
    #*    Print baseline usage report    ##
    #*#####################################
    if (verbose > 0):
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mypy.memory_and_time_report(0, "stdout")

    # computing_env = "debug"    
    # computing_env = "desktop"    
    # computing_env = "server"    
    # computing_env = "small"    
    # computing_env = "medium"    
    # computing_env = "large"    
    
    # read the parameter file into a python dictionary
    with open('../../InitParams/power_method.init_params','r') as init_file:
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

    # print("We are using regionbitparam %d  and zeroregionbitparam %d" %(regionbitparam,zeroregionbitparam))
    # input("STOP HERE")
    
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,LARC_verbose)

    if scalarTypeStr in ('Integer', 'MPInteger'):
        print("\nThis routine does not work with integer types!")
        sys.exit(1)
    if scalarTypeStr in ('MPRational','MPRatComplex'):
        print("\nThis routine does not give good results with this scalarType")
        print("Consider using either MPReal or MPComplex instead")

    # if verbosity higher than WARN start a reporting thread              
    if (verbose > 1):              
        mypy.create_report_thread(report_interval_seconds)

    # Finished with initializing LARC
    if (verbose > 1):             
        print("\nFinished creating LARC matrix and op stores and loading basic matrices.\n")
        print("stopHogging check to see if program is too large to occur once every 10 minutes.\n")


    #*##########################################################
    #* if matrices are too large do not allow naive printing  ##            
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


        
    #*##########################################
    #*  Take the numpy generated random       ##
    #*  orthogonal matrix M of dim=level**2   ##
    #*  and load it into LARC matrix.         ##
    #*##########################################

    userInput = input("\nEnter the matrix with known eigenstuff E into LARC(y/n)?")
    if (userInput == 'n'):
        sys.exit(1)

    # M = rvs(dim)
    # print(M)    
    # print("\nthe last row of this matrix m_last is:")
    # m_last = M[dim-1]
    # print(m_last)
    # print("\nconfirming orthonormality of M via M @ M.T: ")
    # # @ sign is matrix multiply
    # # note: not M*(M.T) - that is element-wise multiplication
    # shouldBeye = M @ M.T
    # print(shouldBeye)
    # print("\nconfirming orthonormality of M via M.T @ M: ")
    # alsoShouldBeye = M.T @ M
    # print(alsoShouldBeye)
    # print("\ngenerating random row vector v of dimension 1x%d" %dim)
    # v = np.random.normal(size=dim)
    # v = np.sort(v)
    # print(v)
    # print("\ngenerating diagonal matrix D with these values")
    # # replaces each element of diagonal with appropriate element of v
    # D = v*np.eye(dim)
    # print(D)
    # print("\ngenerating matrix E with known eigenvalues, eigenvectors")
    # E = (M.T) @ D @ M
    # print(E)
    # print("\nconfirming m_last is eigenvector with largest eigenvalue")
    # w = m_last @ E
    # print("\tw = m_last @ E, eigenvalue should be %g" %v[dim-1])
    # val = w / m_last        # element-wise division
    # for n in range(dim):
    #     print("n=%d: w[n]/m_last[n] = %g" %(n,val[n]))


    Py_matrix= E.tolist()
    if (verbose > 1):
        print("We used tolist to change float64 numpy matrix to float python")
        print(Py_matrix)
    
    userInput = input("\nKeep going(y/n)?")
    if (userInput == 'n'):
        sys.exit(1)

    # Making E into a single list
    if (verbose > 1):
        print("First we make a string of the elements of the matrix.")    
    E_list = []
    F_list = []
    for i in range (0,dim):
        F_list = list(Py_matrix[i])
        E_list = E_list + F_list
    # if (verbose > 1):
    if (verbose > 1):
        print("We make the matrix E into a list E_list")
        print(E_list)

    #a_arr = list(map(str,a_str))
    E_str_list = mypy.map_to_str(E_list,scalarTypeStr)

    # creating or finding the matrix associated with the array
    E_ID = mypy.row_major_list_to_store(E_str_list, level, level, dim)
    mypy.print_naive(E_ID)
    print("\n We have loaded the matrix E into LARC and it has matrixID %d\n"%E_ID)
    print("\n")

    # We will make a simple vector to start
    v = [0 for i in range (0,dim)]
    v[0] = 1
    if scalarTypeStr in ("Complex", "MPComplex", "MPRatComplex"):
        # adding this imaginary part to the input vector prevents
        # problems with the all-real case having a divide-by-zero in
        # recovering the angle of a complex eigenvalue
        v[1] = 1j

    print("\n We have created the vector v")
    print(v)
        
    #a_arr = list(map(str,a_str))
    v_str_list = mypy.map_to_str(v,scalarTypeStr)
    print("\n We tried to turn v into a list of strings")
    print(v_str_list)

    # creating or finding the matrix associated with the array
    v_ID = mypy.row_major_list_to_store(v_str_list, level, 0, 1)
    mypy.print_naive(v_ID)
    print("\n We have loaded the column vector v into LARC and it has matrixID %d\n"%v_ID)
    print("\n")

    v_old_ID = -1
    v_new_ID = v_ID
    # the maximum element in the unit vector is "1", which is pre-stored
    one_ID = mypy.get_identity_pID(0)
    maxEl_ID = one_ID
    counter = 0


    # TODO: Add query for which norm you want to use
    maxnorm = 1  # 1 means use the maxnorm, versus L2 norm
    if (maxnorm):
        chosen_norm = mypy.L_infty
    else:
        chosen_norm = mypy.L_2
           
    # eight_ID = mypy.get_valID_from_valString("8.0")

    # retrieve the next matrix ID to be assigned in LARC v_first_ID
    v_first_ID = mypy.num_matrices_created()
    # make exit condition that v_first_ID < v_new_ID < v_old_ID

    # MARK claims: v_norm decreasing is also a stopping condition

    loop_hash_table = set()
    

    print("\nStarting power method loop!")
    #while (v_old_ID != v_new_ID): # (v_old_ID = v_new_ID) converged  
    while (1):
        v_old_ID = v_new_ID
        old_maxEl_ID = maxEl_ID
        v_unnorm_ID = mypy.matrix_mult(E_ID,v_old_ID)
        maxEl_ID = mypy.matrix_element_with_maxNorm(v_unnorm_ID)
        # eventually, the eigenvalue will be given by the ratio
        # of the new maximum ID element and the old maximum ID element
        eval_ID = mypy.scalar_divide(maxEl_ID,old_maxEl_ID)
        if (maxnorm):
            # dividing by the element with maximal norm makes the
            # value of this element in the vector equal to 1.0
            # (+ 0.0i if complex); it remains maximal
            v_new_ID = mypy.scalar_divide(v_unnorm_ID, maxEl_ID)
            # when using l\infty norm, the vector is now normalized
            maxEl_ID = one_ID
        else:
            # dividing by the element with maximal norm makes the
            # value of this element in the vector equal to 1.0
            # (+ 0.0i if complex); it remains maximal
            v_temp_ID = mypy.scalar_divide(v_unnorm_ID, maxEl_ID)
            # when using l2 norm, we still need to normalize
            norm_ID = mypy.normID(v_temp_ID,mypy.L_2)
            v_new_ID = mypy.scalar_divide(v_temp_ID,norm_ID)
            maxEl_ID = mypy.scalar_divide(one_ID,norm_ID)
            
        # print("The normalized vector has ID %d\n" %normalizedID1)
        # mypy.print_naive(normalizedID1)
        if (v_new_ID == v_old_ID):
            evalue_str = mypy.get_scalar_value_string(eval_ID)
            print("Exiting: new and old eigenvectors are the same!")
            break
        if (v_new_ID < v_old_ID):
            if (v_new_ID in loop_hash_table):
                evalue_str = mypy.get_scalar_value_string(eval_ID)
                print("Exiting: detected loop convergence.")
                break
            else:
                print("\n\tnew is %d, old is %d" %(v_new_ID,v_old_ID))
                print("\tadding ID %d to hash table" %v_new_ID)
                loop_hash_table.add(v_new_ID)
        counter += 1
        if (1): #counter%100 == 0):
            evalue_str = mypy.get_scalar_value_string(eval_ID)
            print("The counter is %d, evalue_str %s"
                %(counter, evalue_str))
            # print("The IDs for v_old is %d and it is" %v_old_ID)
            # mypy.print_naive(v_old_ID)
            print("The IDs for v_new is %d and it is" %v_new_ID)
            mypy.print_naive(v_new_ID)

    print("After loop the counter is %d, the eigenvalue is %s"
          %(counter,evalue_str))
    print("and the eigenvector is")
    mypy.print_naive(v_new_ID)

    v_diff_ID = mypy.matrix_diff(v_new_ID,v_old_ID)
#    v_diff_norm_str = chosen_norm(v_diff_ID)
    v_diff_norm_ID = mypy.normID(v_diff_ID,chosen_norm)
    v_diff_norm_str = mypy.traceID(v_diff_norm_ID)
    print("The difference v_new and v_old has norm %s and the diff vector is " %v_diff_norm_str)
    mypy.print_naive(v_diff_ID)
    
    sys.exit(1)
    
    #*###########################################################
    # create random Toeplitz matrices with level given by input #
    #*###########################################################

    userInput = input("\nCreate a random Toeplitz matrix(y/n)?")
    if (userInput == 'n'):
        sys.exit(1)
        
    dim = 2**level

    if scalarTypeStr in ('Real', 'MPRational', 'MPReal'):
    # if mypy.cvar.scalarTypeDef in ('r', 'q'):
        randVals = [ random.random() for i in range(2*dim-1)]
    elif scalarTypeStr in ('Complex', 'MPRatComplex', 'MPComplex'):
    # elif mypy.cvar.scalarTypeDef in ('c', 'v'):
        randVals = [ np.complex(random.random(),random.random()) for i in range(2*dim-1)]
    elif scalarTypeStr in ('Integer', 'MPInteger'):
    # elif mypy.cvar.scalarTypeDef in ('i', 'z'):
        randVals = [ random.randrange(0,10001) for i in range(2*dim-1)]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # create toeplitz matrix from the 2*dim-1 random numbers
    a = []
    b = randVals
    for i in range(dim):
        a.append(list(b[dim-i-1:2*dim-i-1]))
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
    T_ID = mypy.row_major_list_to_store(arr, level, level, dim)
    filename = "Data/Out/%s_%s.json" %(T_name,scalarTypeStr)
    if (print_naive and (verbose > 1)):
        mypy.print_naive(T_ID)
    print("\nPrinted the Toeplitz matrix %s to file %s." %(T_name,filename))
    T_larcSize = mypy.fprint_larcMatrixFile(T_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(T_name,level,T_larcSize));

   
    
    #*############  THE FOLLOWING ONLY WORKS FOR COMPLEX TYPES #*############

    print("\nIf using a complex scalarType you can continue with FFT example")
    print("which calculates various matrices for block sparse DFT example")
    print("and then multiplies the FFT matrix by the random Toeplitz matrix.")
    userInput = input("Would you like to continue with the example (y/n)?")
    if (userInput == 'n'):
        sys.exit(1)

    
    #*##################################################################
    #*  For the FFT application we have found that the following parameters
    #*  work reasonably well:
    #*    
    #*  SMALL STORES for max_level <= 8:    
    #*      matrix_exponent = 22
    #*      op_exponent = 19   
    #*    
    #*    
    #*  LARGE STORES for max_level > 8:    
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

    # OLD NOTES
        #*###################################################
        #*  Sample failure values for LARC locality hash    #*
        #*  are to set         regionbitparam = 60            #*
        #*  and                zeroregionbitparam = 60      #*
        #*####################################################
        #*  Default values for LARC locality hash are       #*
        #*  both equal to DBL_MANT_DIG -2                   #*
        #*  regionbitparam = -1    # default is 53 bits       #*
        #*  zeroregionbitparam = -1 # OLD default is 1074 bits  #*
        #*  zeroregionbitparam = -1 # OLD default is 1074 bits  #*
        #*  NOTE:  testing shows -z 47 will work            #*
        #*####################################################
        # zeroregionbitparam = 52
        #*####################################################
        #*  TODO: find out the space in which this test fails!!!
        #*  DBL_MANT_DIG is the number of digits in FLT_MANT  #*
        #*  why aren't we using DBL_MANT_BITS  ??????? the number of bits
        #*  used in the mantissa
        #*####################################################

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

    #*##################################################
    #*   Print description sparse block algorithm     #*
    #*   from Cooley-Tukey Radix-2 Factorization      #*
    #*   see Van Loan, Computational Frameworks for   #*
    #*   the Fast Fourier Transform, p.21             #*
    #*##################################################
    if (verbose > 1):
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


    #*##############################
    # inverse permutation matrices #
    #*##############################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nPI_0 matrix is:")
        PI_0 = mypy.create_invShufMat(0)
        if print_naive:
            mypy.print_naive(PI_0)
            
        print("\nPI_1 matrix is:")
        PI_1 = mypy.create_invShufMat(1)
        if print_naive:
            mypy.print_naive(PI_1)
            
        print("\nPI_2 matrix is:")
        PI_2 = mypy.create_invShufMat(2)
        if print_naive:
            mypy.print_naive(PI_2)
                    
        print("\nPI_3 matrix is:")
        PI_3 = mypy.create_invShufMat(3)
        if print_naive:
            mypy.print_naive(PI_3)

    PI_name = "PI_"+str(level)        
    PI_ID = mypy.create_invShufMat(level)
    print("\nFormed %s the inverse shuffle matrix for a DFT: it has matrixID %d." %(PI_name,PI_ID))
    filename = "Data/Out/%s_%s.json" %(PI_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive(PI_ID)
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(PI_name,filename))
    PI_larcSize = mypy.fprint_larcMatrixFile(PI_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(PI_name,level,PI_larcSize));


    #*#######################
    # print roots of unity  #
    #*#######################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nTesting calculation for roots of unity:")
        mypy.print_pow2_roots_unity(1)
        mypy.print_pow2_roots_unity(2)
        mypy.print_pow2_roots_unity(3)


    #*#############################
    # create D matrices in python #
    #*#############################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nD_0 matrix is:")
        D_0 = mypy.create_FFT_DMat(0)
        if print_naive:
            mypy.print_naive(D_0)

        print("\nD_1 matrix is:")
        D_1 = mypy.create_FFT_DMat(1)
        if print_naive:
            mypy.print_naive(D_1)

        print("\nD_2 matrix is:")
        D_2 = mypy.create_FFT_DMat(2)
        if print_naive:
            mypy.print_naive(D_2)

        print("\nD_3 matrix is:")
        D_3 = mypy.create_FFT_DMat(3)
        if print_naive:
            mypy.print_naive(D_3)

    D_name = "D_"+str(level)        
    D_ID = mypy.create_FFT_DMat(level)
    print("\nFormed %s the D matrix for an DFT; it has matrixID %d." %(D_name,D_ID))
    filename = "Data/Out/%s_%s.json" %(D_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive(D_ID)
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(D_name,filename))
    D_larcSize = mypy.fprint_larcMatrixFile(D_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(D_name,level,D_larcSize));


    #*#############################
    # create C matrices in python #
    #*#############################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nC_1 matrix is:")
        C_1 = mypy.create_FFT_CMat(1)
        if print_naive:
            mypy.print_naive(C_1)

        print("\nC_2 matrix is:")
        C_2 = mypy.create_FFT_CMat(2)
        if print_naive:
            mypy.print_naive(C_2)

        print("\nC_3 matrix is:")
        C_3 = mypy.create_FFT_CMat(3)
        if print_naive:
            mypy.print_naive(C_3)

    C_name = "C_"+str(level)        
    C_ID = mypy.create_FFT_CMat(level)
    print("\nFormed %s the C matrix for a DFT; it has matrixID %d." %(C_name,C_ID))
    filename = "Data/Out/%s_%s.json" %(C_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive(C_ID)
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(C_name,filename))
    C_larcSize = mypy.fprint_larcMatrixFile(C_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(C_name,level,C_larcSize));

        
    #*###############################
    # create FFT matrices in python #
    #*###############################
    if (verbose == 3): # DEBUG mode (verbose = 3)
        print("\nF_1 matrix is:")
        F_1 = mypy.create_FFTMat(1)
        if print_naive:
            mypy.print_naive(F_1)

        print("\nF_2 matrix is:")
        F_2 = mypy.create_FFTMat(2)
        if print_naive:
            mypy.print_naive(F_2)

        print("\nF_3 matrix is:")
        F_3 = mypy.create_FFTMat(3)
        if print_naive:
            mypy.print_naive(F_3)

    F_name = "F_"+str(level)        
    F_ID = mypy.create_FFTMat(level)
    print("\nFormed %s the DFT matrix; it has matrixID %d." %(F_name,F_ID))
    filename = "Data/Out/%s_%s.json" %(F_name,scalarTypeStr) 
    if (print_naive and (verbose > 1)):
        mypy.print_naive(F_ID)  
    
    if (verbose > 1):
        print("Printed %s into LARC compressed formated json file %s." %(F_name,filename))
    F_larcSize = mypy.fprint_larcMatrixFile(F_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(F_name,level,F_larcSize));


    #*#########################################################
    # apply inverse perm matrix PI_l to a Toeplitz matrix T_l #
    # where l is the level.                                   #
    #*#########################################################
    print("\nNow multiplying the inverse shuffle matrix %s by the Toeplitz matrix %s" %(PI_name,T_name))
    PIT_name = "PIT_"+str(level)
    PIT_ID = mypy.matrix_mult(PI_ID,T_ID)
    if (print_naive and (verbose > 1)):
        mypy.print_naive(PIT_ID)
    filename = "Data/Out/%s_%s.json" %(PIT_name,scalarTypeStr)
    PIT_larcSize = mypy.fprint_larcMatrixFile(PIT_ID,filename)
    print("The %s matrix of level %d has larcSize %d" %(PIT_name,level,PIT_larcSize));

