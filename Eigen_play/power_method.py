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
 #*##############################################################*#

from __future__ import print_function
import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
from ctypes import *
import json    # for loading parameter files
import numpy as np
import random


## \file power_method.py
#
# \brief Uses the power method to find an extreme eigenvalue/vector
#
# This program uses the power method to find the eigenvalue of a matrix M with
# maximal norm values, along with its associated eigenvector. M is designed to
# be real symmetric, and the starting vector is real, so computation does not
# require complex math. M can be chosen to have known integer-valued
# eigenvalues, or to have random eigenvalues; in either case, the eigenvectors
# of M are random.

if __name__ == '__main__':


    #*#############################################################*#
    #* Commandline arguments are:                                  *#
    #*   level = (matrices are 2**level by 2**level),              *#
    #*   verbose = the local verbosity for this program, and       *#
    #*   LARC_verbose = verbosity passed to LARC initialization.   *#
    #* where: 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL          *#
    #*#############################################################*#

    if 5 == len(sys.argv):
        level = int(sys.argv[1]) 
        verbose = int(sys.argv[2])
        LARC_verbose = int(sys.argv[3])
        nonrandom = int(sys.argv[4])

    else:
        print("\nThis program requires four commandline integer inputs:")
        print("   level:  matrices will be 2**level by 2**level")
        print("   verbose:  the verbosity level for this program")
        print("   LARC_verbose:  the verbosity level for the LARC package")
        print("   nonrandom = 1 uses a matrix with known eigenvalues")
        print("             = 0 uses a matrix with random eigenvalues")
        mypy.explain_verbosity()
        print("\nSample Usage:")
        print(" python power_method.py 3 2 0 1")
        print("  level=3  means we will use 8x8 matrices,")
        print("  verbose=2  prints CHATTY comments in this python routine,")
        print("  LARC_verbose=0  prints only any errors from LARC package")
        print("  nonrandom=1  uses a matrix with known eigenvalues")
        print(" python power_method.py 4 3 1 0")
        print("  level=4  means we will use 16x16 matrices,")
        print("  verbose=3  prints DEBUG comments in this python routine,")
        print("  LARC_verbose=1  prints warnings/errors from LARC package")
        print("  nonrandom=0  uses a matrix with random eigenvalues\n")
        sys.exit()


    dim = 2**level

    if (verbose > 1):
       print("\nYou have chosen parameters:")
       print("\t%d = maximum level of matrices (size up to %d by %d)"
              %(level,dim,dim))
       print("\t%d = verbosity level for this routine" %verbose)
       print("\t%d = verbosity level for LARC code" %LARC_verbose)
       mypy.explain_verbosity()
        

    #*##################################################*#
    #*  General Description of What this Program Does   *# 
    #*##################################################*#
    # input("\nHit return to continue")
    if (verbose > 1):
        print("\nThis code implements the power method for finding the largest")
        print("magnitude eigenvalue of the starting matrix E")
        print("The starting matrix E will be randomly generated in such")
        print("a way that we know its eigenvalues and eigenvectors.")
        
        
    #*##################################################################*#
    #*  Generate a random orthonormal matrix using Householder method   *# 
    #*##################################################################*#
    verbose_section = 0    
    if (verbose > 1):
        print("\nOur first step toward producing a random starting matrix E")
        print("with known largest eigenvalue and associated eigenvector is to")
        print("create a random orthonormal matrix using the Householder method")
        print("\tEnter letter 'd' followed by <Return> to see the details,")
        print("\tor enter <Return> to skip to the results of this section.")
        user_input = input()
        if (user_input == 'd'):
            verbose_section = 1
        

    # The following code is from StackOverflow, and returns
    # a random orthonormal square matrix of given dimension dim
    if (verbose_section):
        print("\nGenerating dimension %d random orthonormal square matrix" %dim)
        print("using numpy with Householder transformation method.")
        print("This numpy code is from StackOverflow.\n")
        
    ## \brief creates a random orthonormal matrix of the input size
    #  \param dim The dimension of the matrix produced
    def rvs(dim=3):
        # random_state = np.random
        H = np.eye(dim)                 # identity of size dim,dim
        if (verbose_section):
            print("H matrix starts as an identity:")
            print(H)
            print("H will be modified until it is an orthonormal matrix")
        D = np.ones((dim,))             # all-ones vector of size dim (?)
        if (verbose_section):
            print("D vector starts as all ones:")
            print(D)
            print("")
        for n in range(1, dim):
            # x = random_state.normal(size=(dim-n+1,))
            x = np.random.normal(size=(dim-n+1,))
            if (verbose_section):
                print("Step n=%d:" %n)
                print("Sampling from normal distribution")
                print("select array x of length %d" %(dim-n+1))
                print(x)
                print("")
            D[n-1] = np.sign(x[0])
            if (verbose_section):
                print("D[%d] is changed to sign of x[0] so D is now:" %(n-1))
                print(D)
                print("")
            # np.sqrt(x*x).sum() is sqrt of sum of squares of elements of x
            x[0] -= D[n-1]*np.sqrt((x*x).sum())
            # this causes x[0] to change sign
            if (verbose_section):
                print("if x[0] was negative then add the norm of x to x[0]")
                print("other wise subtract the norm of x from x[0], giving:")
                print(x)
                print("")
            # Householder transformation
            # take identity subtract 2 times normalized outer product of x
            Hx = (np.eye(dim-n+1) - 2.*np.outer(x, x)/(x*x).sum())
            if (verbose_section):
                print("Hx=(np.eye(dim-n+1)-2.*np.outer(x, x)/(x*x).sum()) is:")
                print(Hx)
                print("")
            mat = np.eye(dim)
            # replace a subblock in bottom right of Identity with Hx
            mat[n-1:, n-1:,] = Hx
            if (verbose_section):
                print("mat=an identity matrix with one corner replaced by Hx:")
                print(mat)
                print("")
            # matrix multiply H times mat
            H = np.dot(H, mat)
            if (verbose_section):
                print("replace H with H times mat:")
                print(H)
                print("")
                input("hit return to continue loop for Householder method")
        # fix the last sign such that the determinant is 1
        # dim%2 is dim mod 2 and D.prod() is the product of all elements in D
        # D[-1] is just reverse indexing and is the last term of D
        D[-1] = (-1)**(1-(dim%2))*D.prod()
        if (verbose_section):
            print("\nAfter loop D is changed using numpy commmand")
            print("D[-1] = (-1)**(1-(dim%2))*D.prod() and becomes:")
            print(D)
            print("")
        # Equivalent to np.dot(np.diag(D), H) but faster, apparently
        # H.T is the transpose of H,
        # v*B is a matrix with row i of B replaced by v_i times that row.
        H = (D*H.T).T
        if (verbose_section):
            print("The final computation is H = (D*H.T).T)")
            print("H should be an orthnormal matrix.\n")
        return H

    H = rvs(dim)
    if (verbose > 1):
        print("The random orthonormal matrix H produced by the numpy code")
        print("implementing the Householder method is:")
        print(H)
        print("")
        input("hit return to continue")

        
    #*##################################################################*#
    #*  Use random orthogonormal matrix H, and diagonal matrix D with   *# 
    #*  selected eigenvectors to create  E =                            *# 
    #*##################################################################*#
    verbose_section = 0    
    if (verbose > 1):
        print("\nOur next step toward producing a random starting matrix E")
        print("with known largest eigenvalue and associated eigenvector is to")
        print("take a diagonal matrix D with our selected eigenvalues, and")
        print("modify it with our random orthonormal matrix H to create")
        print("The matrix E with known eigenvalues (the diagonal of D)")
        print("and known eigenvectors (the columns of H)")
        print("\tEnter letter 'd' followed by <Return> to see the details,")
        print("\tor enter <Return> to skip to the results of this section.")
        user_input = input()
        if (user_input == 'd'):
            verbose_section = 1
                  
    M = H

    # YOU ARE HERE!   
    if (verbose_section):
       print(M)
       print("\nthe last row of this matrix m_last is:")
    m_last = M[dim-1]
    if (verbose_section):
        print(m_last)
        print("\nconfirming orthonormality of M via M @ M.T: ")
        print("@ is numpy for matrix multiply.")
    # @ sign is matrix multiply
    # note: not M*(M.T) - that is element-wise multiplication
    shouldBeye = M @ M.T
    if (verbose_section):
        print(shouldBeye)
        print("\nconfirming orthonormality of M via M.T @ M: ")
    alsoShouldBeye = M.T @ M
    if (verbose_section):
        print(alsoShouldBeye)
    # nonrandom = 0
    if nonrandom:
        print("\ngenerating nonrandom row vector v of dimension 1x%d" %dim)
        print("\nThe vector entries will be the numbers 1 through %d" %dim)
        v = [i for i in range (1,dim+1)]
    else:
        print("\ngenerating random row vector v of dimension 1x%d" %dim)
        print("\nThe vector entries are:")
        v = np.random.normal(size=dim)
        v = np.sort(v)
        print(v)
    print("\nand will be the eigenvalues for the matrix")
    if (verbose_section):
        #print(v)
        print("\ngenerating diagonal matrix D with these values")
    # replaces each element of diagonal with appropriate element of v
    D = v*np.eye(dim)
    if (verbose_section):
        print(D)
    print("\ngenerating matrix E with known eigenvalues, eigenvectors")
    E = (M.T) @ D @ M
    print(E)
    if (verbose_section):
        print("\nconfirming m_last is eigenvector with largest eigenvalue")
        if (not nonrandom):
            print("(note that the smallest eigenvalue is %g)" %v[0])
        w = m_last @ E
        print("\tw = m_last @ E, eigenvalue should be %g" %v[dim-1])
        val = w / m_last        # element-wise division
        for n in range(dim):
            print("n=%d: w[n]/m_last[n] = %g" %(n,val[n]))

    #*####################*#
    #* LARC initalization *#
    #*####################*#
    userInput = input("\nEnter 'q'<RETURN> to quit, <RETURN> to initialize LARC and continue: ")
    if (userInput == 'q'):
        sys.exit(1)
    
        
    #*##################################################################*#
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*##################################################################*#
    scalarTypeStr = mypy.cvar.scalarTypeStr



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

    userInput= input("Press <Enter> to continue\n")
    

    #*###################################*#
    #*    Print baseline usage report    *#
    #*###################################*#
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
    with open('../InitParams/power_method.init_params','r') as init_file:
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

    # warn if the commandline value for LARC_verbose differs from the parameter file value for p_verbose        
    if (LARC_verbose > 0):
        if (LARC_verbose != p_verbose):
            print("NOTE: This program uses commandline (LARC_verbose = %d) " %LARC_verbose)
            print("      rather than the parameter file (p_verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    if (level > max_level):
        print("ERROR: max_level for machine is %d, but input level is %d" %(max_level, level))
        print("consider moving to a machine with more memory")
        sys.exit(0)

    # initialize LARC

    # print("We are using regionbitparam %d  and zeroregionbitparam %d" %(regionbitparam,zeroregionbitparam))
    # input("STOP HERE")
    
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,LARC_verbose)

    if scalarTypeStr in ('Integer', 'MPInteger'):
        print("\nThis routine does not work with integer types!")
        sys.exit(1)
    if scalarTypeStr in ('MPRational','MPRatComplex'):
        print("\nThis routine does not give good results with rational scalarTypes")
        print("Consider using either MPReal or MPComplex instead")
        input("\nHit return to continue, or ^C to quit")

    # if verbosity higher than WARN start a reporting thread              
    if (verbose > 1):              
        mypy.create_report_thread(report_interval_seconds)

    # Finished with initializing LARC
    if (verbose > 1):             
        print("\nFinished creating LARC matrix and op stores and loading basic matrices.\n")
        print("StopHogging check to see if program is too large to occur once every 10 minutes.\n")


    #*########################################################*#
    #* if matrices are too large do not allow naive printing  *#            
    #*########################################################*#
    if (level < 4):
        print_naive = 1
    else:
        print_naive = 0
    if (verbose > 1):             
        if print_naive:
            print("  The level is small enough that we can print files of naive matrices to the screen.")
        else: 
            print("  The level= %d, is too big to reasonable print naive formated matrices to the screen." %level)


        
    #*########################################*#
    #*  Take the numpy generated random       *#
    #*  orthogonal matrix M of dim=level**2   *#
    #*  and load it into LARC matrix.         *#
    #*########################################*#

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
    
