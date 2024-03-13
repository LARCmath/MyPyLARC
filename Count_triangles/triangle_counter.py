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
import MyPyLARC as mypy
from ctypes import *
import json    # for loading parameter files
import numpy as np
import random
import time


if __name__ == '__main__':


    #*####################################################################
    # Set the level (matrices are 2**level by 2**level                   #
    # and the verbosity (0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL) #
    #*####################################################################

    startTime = time.process_time()
    currTime = startTime

    if 5 == len(sys.argv):
        max_level = int(sys.argv[1]) 
        verbose = int(sys.argv[2])
        LARC_verbose = int(sys.argv[3])
        data_path = sys.argv[4]

    else:
        print("\nThis program requires three commandline integer inputs")
        print("and a path to a matrix market data file in sparse format:")
        print("   max_level:  matrices will be up to (2**max_level)^2 size")
        print("   verbose:  the verbosity level for this program")
        print("   LARC_verbose:  the verbosity level for the LARC package")
        print("   data_path:  the path to the matrix market data file.")
        mypy.explain_verbosity()
        print("\nSample Usage:")
        print("  python triangle_counter.py 3 2 0 data/test1.mm")
        print("    max_level=3  means we will use 8x8 matrices,")
        print("    verbose=2  prints CHATTY comments in this python routine,")
        print("    LARC_verbose=0  prints only any errors from LARC package")
        print("    uses matrix market sparse file data/test1.mm")
        print("  python triangle_counter.py 4 3 1 data/test2.mm")
        print("    max_level=4  means we will use 16x16 matrices,")
        print("    verbose=3  prints DEBUG comments in this python routine,")
        print("    LARC_verbose=1  prints warnings/errors from LARC package\n")
        print("    uses matrix market sparse file data/test2.mm")
        sys.exit()


    if (verbose > 1):
        dim = 2**max_level
        print("\nmax_level=%d: Program will use %d by %d matrices."
              %(max_level,dim,dim))
    if ((verbose == 2) or (verbose == 4)):
        print("Verbosity:")
        print("\tThis routine's verbosity is verbose=%d and" %verbose)
        print("\tthe math package's is  LARC_verbose=%d," %LARC_verbose)
        print("\twhere 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL.")


    #*####################################################
    #*  General Description of What this Program Does   #* 
    #*####################################################
    # input("\nHit return to continue")
    if (verbose > 1):
        print("\nThis code counts the number of triangles in a graph")
        print("from the adjacency matrix, by taking the trace of")
        print("the third power of the matrix.")


    #*######################
    #* LARC initalization #*
    #*######################
    if (verbose > 1):
        userInput = input("\nEnter 'q'<RETURN> to quit, <RETURN> to initialize LARC and continue: ")
        if (userInput == 'q'):
            sys.exit(1)
    
        
    #*####################################################################
    #*Figure out the scalarType                                          #
    #*In the Makefile you can compile with different scalarType values   #
    #*Define string for using in formating filenames                     #
    #*####################################################################
    scalarTypeStr = mypy.cvar.scalarTypeStr


    #*####################################################
    #*   Find out if machine has a large amount of      #*
    #*   memory available so we can make bigger tables  #*
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
    #*    Print baseline usage report    #*
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
    with open('../InitParams/count_triangles.init_params','r') as init_file:
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
            p_max_level= p['max_level']
            regionbitparam = p['regionbitparam']
            zeroregionbitparam = p['zeroregionbitparam']
            report_interval_seconds = p['report_interval_seconds']
            min_memGiB_required = p['min_memGiB_required']
            p_verbose = p['verbose']

    # warn if the commandline value for LARC_verbose differs from p_verbose
    # warn if the commandline value for max_level differs from p_max_level
    if (LARC_verbose > 0):
        if (LARC_verbose != p_verbose):
            print("NOTE: This program uses commandline (LARC_verbose = %d) " %LARC_verbose)
            print("      rather than the parameter file (p_verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")
    if (max_level < p_max_level):
        max_level = p_max_level
        if (LARC_verbose > 0):
            print("NOTE: max_level was increased to parameter file value %d "
                  %max_level)
    if (max_level > p_max_level) and (LARC_verbose > 0):
        print("WARN: max_level is larger than parameter file value %d "
              %p_max_level)
        print("and could possibly cause problems if it is a lot bigger.")

    # initialize LARC

    # print("We are using regionbitparam %d  and zeroregionbitparam %d" %(regionbitparam,zeroregionbitparam))
    # input("STOP HERE")
    
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,LARC_verbose)

    if scalarTypeStr in ('Integer', 'MPInteger'):
        print("\nGood idea to use integer type!")
    else:
        print("\nYou are not using an integer, what does that do?")
#        sys.exit(1)
#    if scalarTypeStr in ('MPRational','MPRatComplex'):
#        print("\nThis routine does not work with rational types")
#         print("Consider using either MPReal or MPComplex instead")
#        sys.exit(1)

    # if verbosity higher than WARN start a reporting thread              
    if (verbose > 1):              
        mypy.create_report_thread(report_interval_seconds)

    # Finished with initializing LARC
    if (verbose > 1):             
        print("\nFinished creating LARC matrix and op stores and loading basic matrices.\n")
        print("StopHogging check to see if program is too large to occur once every 10 minutes.\n")


    #*##########################################################
    #* if matrices are too large do not allow naive printing  #*
    #*##########################################################
    if (max_level < 4):
        print_naive = 1
    else:
        print_naive = 0
    if (verbose > 1):             
        if print_naive:
            print("  The level is small enough that we can print files of naive matrices to the screen.")
        else: 
            print("  The level= %d, is too big to reasonable print naive formated matrices to the screen." %max_level)


# testbig.mmio.withheader
    nextTime = time.process_time()
    print("Time for initialization: %g seconds" %(nextTime-currTime))
    currTime = nextTime

    # This file has level 16 so if the initialization max_level is not
    # 2 it will crash:
    # print("WARNING: About to read level 16 matrix.")
    A_ID = mypy.read_matrixMarketExchange_file(data_path)
    print("Read %s and assigned the matrix MatrixID" %data_path)
    print(A_ID)
    size_A = mypy.fprint_larcMatrixFile(
        A_ID,"amat")
    print("The LARCsize of this matrix is ", size_A)

    nextTime = time.process_time()
    print("Time to read in matrix: %g seconds" %(nextTime-currTime))
    currTime = nextTime
    
    # A_str = [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]
    # A_arr = mypy.map_to_str(A_str,scalarTypeStr)

    # dim_whole = 2**level

    # creating or finding the matrix associated with the array
    # A_ID = mypy.row_major_list_to_store(A_arr, level, level, dim_whole)
    # print("The adjacency matrix is:\n")
    # mypy.print_naive(A_ID)
    # print("\n")

    B_ID = mypy.matrix_mult(A_ID,A_ID)
    print("The adjacency matrix squared is:\n")
    mypy.print_naive(B_ID)
    print("\n")
    size_B = mypy.fprint_larcMatrixFile(
        B_ID,"bmat")
    print("The LARCsize of this matrix is ", size_B)
    nextTime = time.process_time()
    print("Time for first matrix multiply: %g seconds" %(nextTime-currTime))
    currTime = nextTime

    C_ID = mypy.matrix_mult(A_ID,B_ID)
    print("The adjacency matrix to the third power is:\n")
    mypy.print_naive(C_ID)
    print("\n")
    size_C = mypy.fprint_larcMatrixFile(
        C_ID,"cmat")
    print("The LARCsize of this matrix is ", size_C)
    nextTime = time.process_time()
    print("Time for second matrix multiply: %g seconds" %(nextTime-currTime))
    currTime = nextTime

    traceStr = mypy.traceID(C_ID)
    print("The trace of the adjacency matrix to the third power is:\n")
    print(traceStr)
    print("\n")

    if (scalarTypeStr in ('Integer','Real','MPInteger','MPReal','MPRational')):
        num_triangles = int(float(traceStr)) // 6
        print("The number of triangles is %d\n" %num_triangles)
    else:
        print("Divide this number by 6 to get the number of triangles:\n")

    nextTime = time.process_time()
    print("Time for trace/division: %g seconds" %(nextTime-currTime))
    print("Total time for calculation: %g seconds" %(nextTime-startTime))
    

    sys.exit(1)

# OLD EXAMPLE WITH FEWER TRIANGLES

#    # This file has level 16 so if the initialization max_level is not
#    # 2 it will crash:
#    print("WARNING: About to read level 16 matrix.")
#    A_ID = mypy.read_matrixMarketExchange_file(data_path)
#    print("Read %s and assigned the matrix MatrixID" %data_path)
#    print(A_ID)
#    
#    # A_str = [0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]
#    # A_arr = mypy.map_to_str(A_str,scalarTypeStr)
#
#    # dim_whole = 2**level
#
#    # creating or finding the matrix associated with the array
#    # A_ID = mypy.row_major_list_to_store(A_arr, level, level, dim_whole)
#    # print("The adjacency matrix is:\n")
#    # mypy.print_naive(A_ID)
#    # print("\n")
#
#    B_ID = mypy.matrix_mult(A_ID,A_ID)
#    print("The adjacency matrix squared is:\n")
#    mypy.print_naive(B_ID)
#    print("\n")
#
#    C_ID = mypy.matrix_mult(A_ID,B_ID)
#    print("The adjacency matrix to the third power is:\n")
#    mypy.print_naive(C_ID)
#    print("\n")
#
#    traceStr = mypy.traceID(C_ID)
#    print("The trace of the adjacency matrix to the third power is:\n")
#    print(traceStr)
#    print("\n")
#    print("Divide this number by 6 to get the number of triangles:\n")
#    print("\n")

    

