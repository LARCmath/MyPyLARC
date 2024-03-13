#!/usr/bin/env python3

 #*##############################################################*#
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
import mpmath


## \file basic_krylov.py
#
# \brief Uses a Krylov subspace method to find the dominant eigenvalues of
#        a matrix
#
# This program takes as input a vector y which has undergone several rounds of
# a power method where the vector is multiplied by a matrix M and then
# renormalized; the input to one round of this process is the output of the
# previous round. (A random input may be provided, but this will result in slow
# convergence of the Krylov subspace method.) Both the input vector y and the
# matrix M are expected to be in LARCMatrix compressed format.

# The program applies a number of rounds of multiplying y by the matrix M
# to generate y, My, M^2y, etc. Rows of each of these vectors are extracted to
# form varying-size matrices A and right-hand-side vectors z. For each size k,
# We solve for c in Ac=-z, and use these c coefficients in the characteristic
# polynomial equation of degree k. If the method works, the roots of these
# polynomials will include the dominant eigenvalues of our M, distinguished
# from spurious values by their consistent appearance as roots for all k-th
# degree polynomials with k greater than some k_0.
#
# If the input matrix M is not complex, the LARC computation involves only real
# numbers, and we recommend compiling LARC to use a multiprecision type,
# preferably MPREAL.
#
# The Python part of the method uses the packages NumPy and mpmath. NumPy is
# used to construct larger A matrices from smaller ones in an efficient manner
# with the hstack and vstack methods. mpmath is needed to maintain 256-bit
# precision within Python. Subroutines from this package used in this
# application include a linear equations solver using the LU decomposition and
# root-finding for the polynomials generated.
#

if __name__ == '__main__':


    #*#############################################################*#
    #* Commandline arguments are:                                  *#
    #*   level = (matrices are 2**level by 2**level),              *#
    #*   verbose = the local verbosity for this program, and       *#
    #*   LARC_verbose = verbosity passed to LARC initialization.   *#
    #*   LARC_verbose = verbosity passed to LARC initialization.   *#
    #* where: 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL          *#
    #*   inmat_name = name of LARCmatrix file providing matrix     *#
    #*   invec_name = name of LARCmatrix file providing initial    *#
    #*                starting point for power method iteration    *#
    #*#############################################################*#

    if 7 == len(sys.argv):
        level = int(sys.argv[1]) 
        verbose = int(sys.argv[2])
        LARC_verbose = int(sys.argv[3])
        inmat_name = sys.argv[4]
        invec_name = sys.argv[5]
        max_iterations = int(sys.argv[6])
    else:
        print("\nThis program requires the commandline inputs:")
        print("   level:  (int) matrices will be 2**level by 2**level")
        print("   verbose:  (int) the verbosity level for this program")
        print("   LARC_verbose: (int) the verbosity level for the LARC package")
        print("   invec_name: the name of the LARCmatrix file with the matrix")
        print("   invec_name: the name of the LARCmatrix file with the starting vector")
        print("   max_iterations: the maximum number of iterations")
        mypy.explain_verbosity()
        print("\nSample Usage:")
        print(" python basic_krylov.py 3 2 0 myMatrix.json startVec.json 50")
        print("  level=3  means we will use 8x8 matrices,")
        print("  verbose=2  prints CHATTY comments in this python routine,")
        print("  LARC_verbose=0  prints only any errors from LARC package")
        print("  inmat_name=myMatrix.json reads in the named file")
        print("  invec_name=startVec.json reads in the named file")
        print("  max_iterations=50  will stop after 50 iterations")
        sys.exit()

    dim = 2**level

    if (verbose > 1):
       print("\nYou have chosen parameters:")
       print("\t%d = maximum level of matrices (size up to %d by %d)"
              %(level,dim,dim))
       print("\t%d = verbosity level for this routine" %verbose)
       print("\t%d = verbosity level for LARC code" %LARC_verbose)
       print("\t%s = input matrix file name" %inmat_name)
       print("\t%s = input vector file name" %invec_name)
       print("\t%d = maximum number of iterations" %max_iterations)
       mypy.explain_verbosity()
        
    #*##################################################################*#
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*##################################################################*#
    scalarTypeStr = mypy.cvar.scalarTypeStr
    if scalarTypeStr in ('Integer', 'MPInteger'):
        print("\nThis routine does not work with integer types!")
        sys.exit(1)
    if scalarTypeStr in ('MPRational','MPRatComplex'):
        print("\nThis routine does not give good results with rational scalarTypes")
        print("Consider using either MPReal or MPComplex instead")
        input("\nHit return to continue, or ^C to quit")

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
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,LARC_verbose)

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
            print("  The level= %d, is too big to reasonably print naive formated matrices to the screen." %level)

###################################################3

    mpmath.mp.prec = 256
    dim = 2**level
    dim_ID = mypy.get_valID_from_valString(str(dim))

    # read in the matrix
    powMatrix = mypy.read_larcMatrixFile(inmat_name)
    mypy.set_hold_matrix(powMatrix)

    print("\n We have loaded the matrix into LARC")
    print("and it has matrixIDs %d\n" %powMatrix)

    invec = mypy.read_larcMatrixFile(invec_name)

    print("\n We have loaded the starting vector into LARC")
    print("and it has matrixID %d\n" %invec)

    # savedIDs is the maximum k in our k\times k matrix A, and also the max
    # degree of the polynomials generated and passed to the root finder.
    savedIDs = max_iterations
    r_index = [0]*savedIDs
    vIDvec = [-1]*savedIDs
    vIDvec[0] = invec

    normID = mypy.normID(invec,mypy.L_2)
    print('norm of invec is ',
        mypy.get_readableString_scalar_from_pID_and_coords(normID,0,0))

    # create first (1x1) A matrix

    # get random row index that has invec[0]!=0
    index_dim = dim
    index_values = [i for i in range(index_dim)]
    while (1):
        ival = random.randrange(0,index_dim)
        sID = mypy.get_scalarID_from_pID_and_coords(vIDvec[0],ival,0)
        if (mypy.matrix_is_zero(sID)==0): break

    # remove row index from further consideration
    r_index[0] = index_values.pop(ival)
    index_dim = index_dim - 1

    # make mpmath value from vector element, put into numpy matrix A
    a_str = mypy.get_readableString_scalar_from_pID_and_coords(sID,0,0)
    np_amat = np.matrix(mpmath.mpmathify(a_str))
#    print("np_amat = ",np_amat)

    with open('basic_krylov.out', 'w') as f:

        print("\nStarting power method loop!")
        for counter in range(1,savedIDs):
            print("counter = ",counter)
    
            # multiply by power method matrices to get next vector
            v_temp_ID = mypy.matrix_mult(powMatrix,vIDvec[counter-1])
            vIDvec[counter] = v_temp_ID
#            print("created vIDvec[",counter,"]")
            if mypy.matrix_is_zero(vIDvec[counter]):
                print("\tpower method produced zero vector!")
                sys.exit(0)

            # keep track of the norms of these vectors as they decrease
            normID = mypy.normID(vIDvec[counter],mypy.L_2)
#            print('norm of vIDvec[', counter, '] is ',
#                mypy.get_readableString_scalar_from_pID_and_coords(normID,0,0))

            # get numpy col vector of dimension counter-1 containing mpf data
            z_ID = [mypy.get_scalarID_from_pID_and_coords(vIDvec[counter],
                r_index[i],0) for i in range(counter)]
            z_str = [mypy.get_readableString_scalar_from_pID_and_coords(z_ID[i],
                0, 0) for i in range(counter)]
            np_cvec = np.matrix([[mpmath.mpmathify(z_str[i])]
                for i in range(counter)])
    #        print("np_cvec = ",np_cvec)
    #        print("np_amat = ",np_amat)

#           convert numpy matrix/vector to mpmath matrix/vector
            mp_cvec = mpmath.matrix(np_cvec)
            mp_amat = mpmath.matrix(np_amat)
    #        print("mp_cvec = ",mp_cvec)
    #        print("mp_amat = ",mp_amat)

#           solve for roots of characteristic polynomial
            if (counter>1):
#               find c in Ac=-z
                coeffs = mpmath.lu_solve(mp_amat,mp_cvec)
                coeffs = coeffs * -1
    #            print("coeffs = ",coeffs)

                # create vector describing characteristic polynomial
                polyvec = mpmath.matrix(1,counter+1)
                polyvec[0] = mpmath.mpmathify('1.0')
                for i in range(1,counter+1):
                    polyvec[i] = coeffs[counter-i]
#                print("polyvec = ",polyvec)

                # documentation says polyroots() implements the Durand-Kerner
                # method, which specifically fails for repeated roots. The
                # mpmath implementation avoids divide-by-zero somehow...
                rootvec = mpmath.polyroots(polyvec,maxsteps=1000,extraprec=100)
#                print("rootvec = ",rootvec)
                # since some roots may be complex, find and print absolute
                # values
                rootnorm = [mpmath.fabs(rootvec[i]) for i in range(counter)]
#                print("rootnorm = ",rootnorm)
                for i in range(counter):
                    print(counter,"\t\t", rootnorm[i], rootvec[i], file=f)
                    print(counter,"\t\t", rootnorm[i], rootvec[i])

                f.flush()
            # make A matrix rectangular by adding z onto the right side
            np_amat = np.hstack((np_amat,np_cvec))
    
            # choose new row randomly, making sure that the 
            # square matrix A that results is nonsingular
#            print("index_dim = ",index_dim)

            while (1):
                # pick random row index from list
                # but exclude any row for which vIDvec[counter] has a zero
                while (1):
                    ival = random.randrange(0,index_dim)
                    sID = mypy.get_scalarID_from_pID_and_coords(vIDvec[counter],
                        index_values[ival],0)
                    if (mypy.matrix_is_zero(sID)==0): break
                    elif (index_dim == 1):
                        print("last row rejected, exiting")
                        sys.exit(0)
                    else:
                        print("rejected zero row, index ",
                                index_values.pop(ival))
                        index_dim = index_dim - 1

                # create a row vector from all previous vectors
                # containing the elements at this index
                r_index[counter] = index_values[ival]
                r_ID = [mypy.get_scalarID_from_pID_and_coords(vIDvec[i],
                    r_index[counter],0) for i in range(counter+1)]
                r_str = [mypy.get_readableString_scalar_from_pID_and_coords(
                    r_ID[i], 0, 0) for i in range(counter+1)]

                # add this row to bottom of A, check determinant of
                # resulting square matrix
                np_rvec = [mpmath.mpmathify(r_str[i]) for i in range(counter+1)]
                new_amat = np.vstack((np_amat,np_rvec))
                det = mpmath.det(new_amat)
                print("counter = ",counter,", det = ",det)
                if (det==0.0):
                    print("det = 0, rejected row index ",index_values.pop(ival))
                    index_dim = index_dim - 1
                    if (index_dim > 0): continue # return to top of loop
                    print("no more rows to choose - exiting")
                    sys.exit(0)
                # remove row index from further consideration,
                # continue with new larger A matrix
                print("chose row index ",index_values.pop(ival))
                index_dim = index_dim - 1
                np_amat = new_amat
                # exit loop over random rows
                break
