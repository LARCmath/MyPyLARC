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

import numpy as np
import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
# sys.path.append(os.path.join(os.path.dirname(__file__),"../larc/src"))
# import pylarc
from ctypes import *
import json

if __name__ == '__main__':

    verbose = 1


    ####################################################
    ##   Print welcome message for LARC and MyPyLARC  ##
    ####################################################

    # look in How to Write Python
    # give a menu
    #     0 to see read two page overview of LARC and MyPyLARC
    #     1 getting started in LARC, initialization, io and creating matrices
    #     2 recursive operations in LARC, the matrix store, the operations store
    #     3 writing your own recursive function for LARC

    while True:
        print("**********************************************************************")
        print("*           WELCOME to the LARC Tutorial #1                          *")
        print("*   You can also see a list of contributors, introductory slides,    *")
        print("*   and a paper in larc/doc with more details on LARC and MyPyLARC.  *")
        print("*                                                                    *")
        print("*  Tutorial Menu:                                                    *")
        print("*    0. exit tutorial,                                               *")
        print("*    1. introduction to LARC ideas and MyPyLARC package,             *")
        print("*    2. inputing matrices into LARC and storing output in files      *")
        # print("*    3. LARC initialization from parameter files                     *")
        # print("*    4. recursive operations in LARC, matrix and operations stores   *")
        # print("*    5. writing your own recursive functions.                        *")  
        print("**********************************************************************")
        print("")
        user_input=input("Make a Tutorial selection by entering the item number: ")
        # if user_input in['0','1','2','3','4']:
        if user_input in['0','1','2']:
            break
        else:
            print("Please enter a number from the Tutorial Menu: ")

    if user_input=='0':
        print("\n\tExiting LARC and MyPyLARC Tutorial.\n\tPlease visit again!\n")
        sys.exit()

    if user_input=='1':
        print("")
        print("------------------------------------------------")
        print("LARC (Linear Algebra via Recursive Compression) is a software package")
        print("developed to store 2^r by 2^c matrices in a recursively compressed format")
        print("and to perform operations on the matrices without leaving that format.")
        print("")
        print("It was developed at the Center for Computing Sciences starting in 2013")
        print("and made generally available on GitHub in 2018.  The 2020 release of")
        print("LARC includes a BSD software license and the MyPyLARC demonstration")
        print("package that uses LARC for its vector and matrix calculations.")
        print("")
        print("The computational code for LARC is mostly in C, and there is a Python")
        print("wrapper which combines the functionality of the C package (by")
        print("utilizing SWIG) and some additional Python utilities.")
        print("")
        print("LARC uses quadrant submatrices to produce a recursive representation.")
        print("")
        print("             a  |  b                e  |  f                   A  |  B")
        print("Matrix  A =  -  -  -   Matrix B  =  -  -  -  ... Matrix  M =  -  -  -")
        print("             c  |  d                g  |  h                   C  |  D")
        print("")
        print("Matrix operations are performed recursively and are 'memoized' for reuse.")
        print("")        
        print("**********************************************************************")
        print("")
        yes_no=input("Would you like more details on recursive compression? (y/n): ")
        if yes_no =='y':
            print("")
            print("------------------------------------------------")
            print("LARC uses quadrant submatrices to produce a recursive representation.")
            print("(For similar representations of matrices and operation implementation")
            print("see Wise's quadtree representation for matrix multiplication on")
            print("parallel computers, [1], or the sparse block-recursive matrix FFT of")
            print("Cooley-Tukey, [2] p.21 of Van Loan's Computational Frameworks for the")
            print("FFT).")
            print("")
            print("             A  |  B                     M  |  N")
            print("Matrix  M =  -  -  -  ...   Matrix R  =  -  -  -   ...")
            print("             C  |  D                     P  |  Q")
            print("")
            print("Each new matrix is assigned a unique index (MatrixID) and is stored as")
            print("a list of the dimension exponents followed by the four MatrixID's of")
            print("its quadrant submatrices.  (LARC stores only matrices that are 2^r by")
            print("2^c dimensioned, where the exponent r is called the row_level and the")
            print("exponent c is the column_level.)")
            print("")
            # print("LARC uses several techniques to quickly recognize and")
            # print("retrieve the records of reused matrices and repeated")
            # print("matrix operations.")
            print("Storage and retrieval in the matrix store is via an efficient locality")
            print("hash system of the component MatrixID's with adjustable tolerances for")
            print("the equivalence of underlying scalar constants.")
            print("")
            print("Matrix operations are also defined recursively as functions of quadrant")
            print("submatrices. They are performed recursively and each is memoized into")
            print("a separate Operations store via hashing of the appropriately ordered")
            print("matrix indices with the operation name so that repeated work is")
            print("recognized and need not be performed more than once. (LARC can take")
            print("advantage of the commutativity of matrix addition by sorting the")
            print("matrixIDs of the summands before memoizing the operation.)")
            print("")
            print("All zero matrices and identity matrices have flags in their matrix")
            print("records so that operational identities can be used to short cut")
            print("calculations.  ")
            print("")
            print("LARC is designed to work efficiently on large power-of-two dimensioned")
            print("matrices (2^r x 2^c) which: have tensor structure; are created by block")
            print("algorithms; or have repeated submatrices, including those with a")
            print("limited number of distinct scalars (such as sparse matrices).")
            print("")
            print("LARC is able to impressively compress the matrices in applications of")
            print("interest (beyond the sparse representations available in standard")
            print("packages.)  However, the overhead for its recursive compression makes")
            print("it inefficient if used on large dense matrices with little repeated")
            print("quadtree structure.")
            print("")
            print("Matrices can be input/output in LARC compressed format to/from files")
            print("(using JSON - see below for a JSON minitutorial); this allows staged")
            print("algorithms and check-pointing.")
        print("")
        print("**********************************************************************")
        print("")
        yes_no=input("Would you like more details on LARCs locality preserving hash, finite\nprecision issues, and pseudo symbolic computation (y/n): ")
        if yes_no =='y':
            print("")
            print("------------------------------------------------")
            print("LARC incorporates a special CCS-developed retrieval method for scalars")
            print("called *neighborhood representative retrieval* that addresses finite")
            print("precision issues and allows numbers like 1.999999987 and 2.000000001")
            print("to be stored using the same MatrixID.")
            print("")
            print("LARC divides the space of base elements into small neighborhoods and")
            print("only stores a scalar if it is the first occurring representative of its")
            print("neighborhood.  By combining a neighborhood-defining function and a")
            print("hash function, a locality-preserving hash is created that allows LARC")
            print("to quickly check whether a scalar in a particular neighborhood has")
            print("already been stored.  When an attempt is made to store a scalar, if")
            print("LARC finds a stored record for an existing representative of that")
            print("scalar's neighborhood, then the new scalar is not stored and instead a")
            print("pointer to the existing record is returned.")
            print("")
            print("This significantly reduces (but does not eliminate) the spurious")
            print("spawning of nearly identical matrices.  It also allows LARC to mirror")
            print("and use important mathematical identities via selective pre-loading of")
            print("important scalars, such as the roots of unity for an FFT.")
            print("")
            print("The user specification of neighborhood tolerances allows for an")
            print("application-based balance between the cost of collapsing two critical")
            print("elements into a single element (from too loose a tolerance) and the")
            print("cost of having more than one representor for the same element (from")
            print("too tight a tolerance). (In some applications, having too tight a")
            print("tolerance may reduce the redundancy at the heart of the compression")
            print("scheme and result in overtaxing the system memory.)")
            print("")
            print("The existence of separated unique regional object identifiers")
            print("opens the door for the use of these representors in pseudo-symbolic")
            print("computation.")
            print("")
            print("References")
            print("[1] David Wise 1984")
            print("[2] Van Loan 1960")
            print("[3] Strassen 1969")
            print("")
            
    if user_input=='2':
        print("")
        print("------------------------------------------------")
        print("In this input/output demonstration we initialize the LARC")
        print("matrix store and operations store with some standard parameters,")
        print("then let you see how to get matrices into and out of LARC.")
        print("------------------------------------------------")
        print("")

    print("")
    print("**********************************************************************")
    print("")
    print("We recommend that the first time through this tutorial,")
    print("you respond 'n' to the following question.")
    print("When you want to learn about initialization details respond 'y'")
    print("or look at the code in tut1_larc_overview.py.")
    yes_no=input("Would you like to initialize in verbose mode?: (y/n) ")
    if yes_no =='n':
        verbose = 0
        # print("")
        # print("we don't care if you say no ...")
        # print("")
    elif yes_no == 'y':
        verbose = 1
    else:
        print("\nThe current value of verbose is %d\n" %verbose)

    print("")



    ######################################################
    ##   Find out if machine has a large amount of      ##
    ##   memory available so we can make bigger tables  ##
    ######################################################
    memory_available = mypy.memory_available_GiB()
    if verbose:
        print("\nThe memory available is %ld GiB" %memory_available)
        print("We will use this to select which computing_env to read from parameter file.")
        print("You could write code to select computing_env automatically.")

    if (memory_available > 200):
        if verbose:
            print("\nThis memory is more than 200 GiB\n")
        computing_env = 'large'
    else:    
        if (memory_available > 50):
            if verbose:
                print("\nThis memory is between 50 and 200 GiB\n")
            computing_env = 'medium'
        else:
            if verbose:
                print("\nThis memory is less than 50 GiB\n")
            computing_env = 'small'
            
    if verbose:
        print("This program believes the computing_environment is %s" %computing_env)

    userInput= input("Press <Enter> to continue\n")

    #######################################
    ##    Print baseline usage report    ##
    #######################################
    if verbose:
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mypy.rusage_report(0, "stdout")

        
    ##############################################
    ## read a parameter file into a dictionary  ##
    ##############################################
    if verbose:
        print("")
        print("The parameters echoed below can be set with python code:")    
        print("2^matrix_exponent is the size of the hash table storing matrices")
        print("2^op_exponent is the hash table size for remembered matrix operations")
        print("max_level is the deepest level of recursion allowed. Always a power of 2. RMB IS THIS TRUE")
        print("Values within 2^(-roundsigbits) of stored values are given the same representation")
        print("Values withing 2^(-trunc_to_zero_bits) of 0 are given the special 0 representation.")    
        print("report_interval_seconds provides a 'heartbeat' interval for screen reports.")
        print("verbose varies the amount of information reported")

   
    with open('../InitParams/tutorial.init_params','r') as init_file:
        init_param = json.load(init_file)
        if verbose:
            print("We have read these initialization parameters from a parameter file:")
            for p in init_param[computing_env]:
                print('  MatrixExponent: %d' %(p['matrix_exponent']))
                print('  OpExponent: %d' %(p['op_exponent']))
                print('  MaxLevel: %d' %(p['max_level']))
                print('  RoundSigBits: %d' %(p['rnd_sig_bits']))
                print('  TruncToZeroBits: %d' %(p['trunc_to_zero_bits']))
                print('  ReportIntervalSecs: %d' %(p['report_interval_seconds']))
                print('  MinMemRequiredGiB: %d' %(p['min_memGiB_required']))
                print('  parameter file verbose (ignored) is: %d' %(p['verbose']))
                print('')
        for p in init_param[computing_env]:
            matrix_exponent = p['matrix_exponent']
            op_exponent = p['op_exponent']
            max_level= p['max_level']
            rnd_sig_bits = p['rnd_sig_bits']
            trunc_to_zero_bits = p['trunc_to_zero_bits']
            report_interval_seconds = p['report_interval_seconds']
            min_memGiB_required = p['min_memGiB_required']
            p_verbose = p['verbose']
            # verbose = p['verbose']
            
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)
    if verbose:
        mypy.create_report_thread(report_interval_seconds)
    print_naive = 0
    print_nonzeros = 0

    if verbose:
        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("Kill check to see if program is too large to occur once every 10 minutes.\n")
        #print("LARC initialization is complete.\n")
        print("LARC initialization is complete.")
        print("")
        print("**********************************************************************")
        print("")

    if verbose:
        print("\nWe initialized in verbose mode.")
    else:
        print("\nWe initialized in non-verbose mode.")

        
    ##############################################################
    ##  In the Makefile you can compile with:                   ##
    ##   TYPE=INTEGER, TYPE=REAL, TYPE=COMPLEX,                 ##
    ##  or with multiprecision types:                           ##
    ##   TYPE=MPINTEGER, TYPE=MPRATIONAL, or TYPE=MPRatComplex  ##
    ##############################################################
    scalarTypeStr = mypy.cvar.scalarTypeStr
        
    ##  describe scalarTypes
    if verbose:
        print("\nAt compile time, you can specify which scalarType will be used")
        print("by using the TYPE modifier, e.g.:  make TYPE=COMPLEX.")
        print("The available types are:")
        print("\tINTEGER       (C int64_t)")
        print("\tREAL          (C long double, this is the default scalarType)")
        print("\tCOMPLEX       (C long complex)")
        print("\tMPINTEGER     (GMP multiprecision integer)")
        print("\tMPRATIONAL    (GMP multiprecision rational numbers)")
        print("\tMPRATCOMPLEX  (LARC structure with real and imag copies of GMP rational)")
        print("\nWhen writing python in which the algorithm depends on the")
        print("scalarType, one can query for the current scalarType using:")
        print("\t import MyPyLARC as mypy")
        print("\t mypy.cvar.scalarTypeStr")
        print("which returns a string specifying the scalarType, respectively:")
        print("Integer, Real, Complex, multi-precision (mp) MPInteger,")
        print("MPRational, and larc complex rational MPRatComplex.")
        

    ########################################
    # build vectors and matrices in python #
    ########################################        
    yes_no=input("Would you like to see some examples of matrix formation? (y/n): ")
    if yes_no =='y':
        print("")
        print("------------------------------------------------")
        print("A vector can be built using python 'list', 'map', 'str' commands along with")
        print("the LARC function mypy.row_major_list_to_store_matrixID.  LARC stores")
        print("only matrices that are 2^r by 2^c dimensioned, where the exponent r")
        print("is called the row_level and the exponent c is the column level.")
        print("For example, if one desires a column vector with 4 rows from [-1,0,1,2]")
        print("it is a 4 by 1 matrix so row_level = 2 and col_level = 0.")
        print("")
        print("(The .py source code of this tutorial has the following 5 lines)")
        print("A_arr = list(map(str,[-1,0,1,2]))")
        print("rowLevel = 2")
        print("colLevel = 0")
        print("length_row = 1 << colLevel   # length_row = 2^(colLevel)")
        print("A_mID = mypy.row_major_list_to_store_matrixID(A_arr,rowLevel,colLevel,length_row)")
        A_arr = list(map(str,[-1,0,1,2]))
        rowLevel = 2
        colLevel = 0
        length_row = 1 << colLevel   # length_row = 2^(colLevel)
        A_mID = mypy.row_major_list_to_store_matrixID(A_arr,rowLevel,colLevel,length_row)
        print("")
        print("Since this matrix is tiny one can print it out using")
        print("'mypy.print_naive_by_matID(A_mID)' : ")
        
        mypy.print_naive_by_matID(A_mID)
        print("")
        print("A 2 by 2 matrix can also be built with row_level = 1 and col_level = 1.")
        print("from the same input data [-1,0,1,2]:")
        print("")
        print("(The .py source code of this tutorial has the following 5 lines)")
        print("B_arr = list(map(str,[-1,0,1,2]))")
        print("rowLevel = 1")
        print("colLevel = 1")
        print("length_row = 1 << colLevel")
        print("B_mID = mypy.row_major_list_to_store_matrixID(B_arr,rowLevel,colLevel,length_row)")
        print("")
        B_arr = list(map(str,[-1,0,1,2]))
        rowLevel = 1
        colLevel = 1
        length_row = 1 << colLevel
        B_mID = mypy.row_major_list_to_store_matrixID(B_arr,rowLevel,colLevel,length_row)
        print("Since this matrix is tiny one can print it out using")
        print("'mypy.print_naive_by_matID(B_mID)' :")
        mypy.print_naive_by_matID(B_mID)
        
        
    print("")    
    print("**********************************************************************")
    print("")
    yes_no=input("Would you like to continue with examples of matrix file input and storage? (y/n): ")
    if yes_no =='y':
       # TESTING READING AND WRITING OF MATRICES
        filename_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr
        filename_naive = "Data/Out/sample.1.1.%s.naive" %scalarTypeStr
        filename_json = "Data/Out/sample.1.1.%s.json" %scalarTypeStr
        print("")   
        print("LARC can write matrices using three file formats, usually indicated")
        print("by the suffix:  '.rmm', '.naive', and '.json'.")
        print("\t'.rmm' refers to a row-major matrix format;")
        print("\t'.naive' to a standard matrix representation which lists all the")
        print("\t\telements of the matrix; and")   
        print("\t'.json' to a compressed LARC matrix format in a json container.")
        print("LARC can read matrices using three file formats:")
        print("\trow-major matrix format;")
        print("\tLARC compressed matrix format in json container; and")
        print("\tMatrix Market Exchange format.")
        print("")
        print("Examples:")
        print("")
        print('For instance, if scalarType = "Real" has been set then one can read ')
        print(" a 2 by 2 matrix of [[1, 0],[.400000000000222, .200000000000111]] from file: ")
        print(os.path.join(os.path.dirname(__file__),filename_rmm))
        print(" and write to the file: ")
        print(os.path.join(os.path.dirname(__file__),filename_naive))
        print("using the following four lines of python code:")
        print('   filename_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr ')
        print('   filename_naive = "Data/Out/sample.1.1.%s.naive" %scalarTypeStr ')        
        print("   sample_matrixID = mypy.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm))")
        print("   mypy.write_naive_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_naive))")
        print("")        
        print("Reading row major matrix format from file %s." %filename_rmm)
        sample_matrixID = mypy.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm))       
        print("")
        print("Writing the matrix in naive format to file %s." %filename_naive)
        print("(This prints every entry in the matrix so it should be used only for small matrices.)")
        mypy.write_naive_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_naive))
        print("")
        print("One can also print the sample matrix in naive format to the screen using ")
        print("   mypy.print_naive_by_matID(sample_matrixID)")
        print("")
        mypy.print_naive_by_matID(sample_matrixID)    
        print("")
        print("The following short JSON (JavaScript Object Notation) minitutorial")
        print("may be useful in understanding the formatting underlying LARC compressed files:")
        print("    Data is in name/value pairs")
        print("    	 (A name/value pair is a field name (in double quotes),")
        print("	 followed by a colon, followed by a value. eg: \"size\",10")
        print("    	    [JSON names require double quotes; JSON values can be: ")
        print("    	    number (integer or floating), string (double quoted)")
        print("    	    Boolean (true or false), array (in square brackets)")
        print("    	    object (in curly braces), null.]")
        print("    Data is separated by commas")
        print("    Curly braces hold objects")
        print("    Square brackets hold arrays")
        print("")
        print("Finally, to write the matrix in LARC recursively compressed format")
        print("to file %s , one can use the two python lines: " %filename_json)
        print('  filename_json = "Data/Out/sample.1.1.%s.json" %scalarTypeStr ')
        print("  mypy.write_larcMatrix_file_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_json))")
        print("")        
        print("(Several examples of various sizes can be seen in the directory Data/In.)")
        mypy.write_larcMatrix_file_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_json))
        
        print("")        
        print('Using  os.system("more %s" %filename_json)  , One can see that the')
        print("LARC recursively compressed file looks like this:")
        os.system("more %s" %filename_json)
        
        print("A 'json' container above holds the data associated to the matrix.")
        print("(Thus the first and last lines of the file are curly braces.)")
        print(" matrixID_max is the next assignable matrixID in the matrix store. ")
        print(" matid is the matrixID which was assigned to the matrix represented in this file")
        print("The table above uses one line for each submatrix defining the matrix,")
        print("Each line starts with the matrixID of the submatrix, then two numbers ")
        print(" representing the log base 2 of the row and column dimensions:")
        print(" ('row_level' and 'column_level').")
        print("If these levels are both 0, then the next value is a scalar.")
        print("Otherwise, there is a list of four matrixIDs which correspond to other submatrices")
        print(" which have already been recursively defined (i.e., have already been printed.)")
        print('The special final name-value pair "end":0 signifies the end of the table.')
        print("")
        print("**********************************************************************")
        print("")
        
        ##################################
        # read a row major formated file #
        ##################################
        
        
        
        #####################################
        # read a MatrixMarket formated file #
        #####################################
        
        


        #######################################################
        # Read a LARC matrix file in our own recursive format #
        #######################################################


        # print("*    3. LARC initialization from parameter files                     *")
        # print("*    4. recursive operations in LARC, matrix and operations stores   *")
        # print("*    5. writing your own recursive functions.                        *")  
 
