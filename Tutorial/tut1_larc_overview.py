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

##
# \file tut1_larc_overview.py
#
# \brief This program illustrates LARC initialization
# and some basic LARC operations.
#
# The specific operations shown are
# based on the user's response to various prompts.
# With each operation, this program provides an explanation
# of what is happening under the covers of the LARC package.

if __name__ == '__main__':

    verbose = 1


    #*##################################################
    #*   Print welcome message for LARC and MyPyLARC  ##
    #*##################################################

    # look in How to Write Python
    # give a menu
    #     0 to see read two page overview of LARC and MyPyLARC
    #     1 getting started in LARC, initialization, io and creating matrices
    #     2 recursive operations in LARC, the matrix store, the operations store
    #     3 writing your own recursive function for LARC

    while True:
        print("**********************************************************************")
        print("*           WELCOME to the MyPyLARC and the LARC Tutorial #1         *")
        print("**********************************************************************")
        print("*                                                                    *")
        print("*  Tutorial Menu:                                                    *")
        print("*    0. exit tutorial,                                               *")
        print("*    1. WORDS:  introduction to LARC ideas and MyPyLARC package,             *")
        print("*    2. DEMO:  initialization, matrix input, storage, operation, and output    *")
        print("*    3. RESOURCES: papers, slides, documentation, tutorial routines, etc.   *")
        print("**********************************************************************")
        print("")
        user_input=input("Make a Tutorial selection by entering the item number: ")
        # if user_input in['0','1','2','3','4']:
        if user_input in['0','1','2','3']:
            break
        else:
            print("Please enter a number from the Tutorial Menu: ")

    if user_input=='0':
        print("\n\tExiting LARC and MyPyLARC Tutorial.\n\tPlease visit again!\n")
        sys.exit()

    if user_input=='3':
        print("**********************************************************************")
        print("")
        mypy.list_explanatory_resources()
        print("")
        Userinput= input("Press <Enter> to continue\n")
        
        # Userinput= input("Press <Enter> to continue\n")
        # print("**********************************************************************")
        # print("*  RESOURCES:                                                               ")
        # print("*   A detailed explanatory paper on the LARC (Linear Algebra via   *")
        # print("*    Recursive Compression) package and the MyPyLARC Tutorial    *")
        # print("*   and Sample Applications package is available at:                       *")
        # print("*             MyPyLARC/doc/LARCandMyPyLARC.pdf                           *")
        # print("*                                                                    *")
        # print("*   Additional information can be found at the following locations:  *")
        # print("*     MyPyLARC/Tutorial/README:  short tutorial routine descriptions *")
        # print("*     MyPyLARC/Tutorial/newuser_instructions:                        *")
        # print("*                                suggested order to do tutorials     *")
        # print("*     MyPyLARC/README.md:        overview of MyPyLARC and LARC,      *")
        # print("*                                compiling, matrix operations list   *")
        # print("*     MyPyLARC/doc:              explanatory paper and poster        *")
        # print("*     MyPyLARC/larc/doc:         intro slides, contributors list     *")
        # print("*     MyPyLARC/html/index.html:  doxygen documentation to view       *")
        # print("*                                in browser                          *")
        # print("**********************************************************************")

    if user_input=='1':
        print("**********************************************************************")
        mypy. introduce_LARC_and_MyPyLARC()
        print("")
        Userinput= input("Press <Enter> to continue\n")
        print("**********************************************************************")
        mypy.explain_matrix_and_operation_storage()
        print("")
        Userinput= input("Press <Enter> to continue\n")
        print("**********************************************************************")
        mypy.explain_level()
        print("")
        Userinput= input("Press <Enter> to continue\n")
        print("**********************************************************************")
        mypy.explain_scalar_techniques()
        print("")
        Userinput= input("Press <Enter> to continue\n")
        print("**********************************************************************")

        # yes_no=input("LARC can also snap similar scalar values together. Learn more? (y/n): ")
        # if yes_no =='y':
        #     print("")
        #     print("------------------------------------------------")

 
  
        # print("")
        # print("**********************************************************************")
        # print("")
        # yes_no=input("Would you like more details on LARCs locality hashes, finite\nprecision issues, and pseudo symbolic computation (y/n): ")
        # if yes_no =='y':
        #     mypy.explain_collapsingScalars()
            
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



    #*####################################################
    #*   Find out if machine has a large amount of      ##
    #*   memory available so we can make bigger tables  ##
    #*####################################################
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

    #*#####################################
    #*    Print baseline usage report    ##
    #*#####################################
    if verbose:
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mypy.memory_and_time_report(0, "stdout")

        
    #*############################################
    #* read a parameter file into a dictionary  ##
    #*############################################
    if verbose:
        print("")
        print("The parameters echoed below can be set with python code:")    
        print("2^matrix_exponent is the size of the hash table storing matrices")
        print("2^op_exponent is the hash table size for remembered matrix operations")
        print("max_level is the deepest level of recursion allowed; it is also the log (base 2) of the dimension of the largest matrix allowed.")
        print("Values within 2^(-roundsigbits) of stored values are given the same representation")
        print("Values withing 2^(-zeroregionbitparam) of 0 are given the special 0 representation.")    
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
                print('  RegionBitParam: %d' %(p['regionbitparam']))
                print('  ZeroRegionBitParam: %d' %(p['zeroregionbitparam']))
                print('  ReportIntervalSecs: %d' %(p['report_interval_seconds']))
                print('  MinMemRequiredGiB: %d' %(p['min_memGiB_required']))
                print('  parameter file verbose (ignored) is: %d' %(p['verbose']))
                print('')
        for p in init_param[computing_env]:
            matrix_exponent = p['matrix_exponent']
            op_exponent = p['op_exponent']
            max_level= p['max_level']
            regionbitparam = p['regionbitparam']
            zeroregionbitparam = p['zeroregionbitparam']
            report_interval_seconds = p['report_interval_seconds']
            min_memGiB_required = p['min_memGiB_required']
            p_verbose = p['verbose']
            # verbose = p['verbose']
            
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)
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

    Userinput= input("Press <Enter> to continue\n")
    print("**********************************************************************")
    mypy.explain_verbosity()
    print("")

    if verbose:
        print("\nWe initialized in verbose mode.")
    else:
        print("\nWe initialized in non-verbose mode.")

    Userinput= input("Press <Enter> to continue\n")        
    print("**********************************************************************")
    # you can retrieve the scalarType you are using with this    
    scalarTypeStr = mypy.cvar.scalarTypeStr
    print("\nLARC is currently compiled with scalarType ")
    print(scalarTypeStr)
    print("")
    mypy.explain_scalarType()
    print("")

    Userinput= input("Press <Enter> to continue\n")
    print("**********************************************************************")
    mypy.print_larc_version();
    print("")

    #*######################################
    # build vectors and matrices in python #
    #*######################################        
    yes_no=input("Would you like to see some examples of matrix formation? (y/n): ")
    if yes_no =='y':
        print("")
        print("------------------------------------------------")
        print("A vector can be built using python 'list', 'map', 'str' commands along with")
        print("the LARC function mypy.row_major_list_to_store.  LARC stores")
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
        print("A_pID = mypy.row_major_list_to_store(A_arr,rowLevel,colLevel,length_row)")
        A_arr = list(map(str,[-1,0,1,2]))
        rowLevel = 2
        colLevel = 0
        length_row = 1 << colLevel   # length_row = 2^(colLevel)
        A_pID = mypy.row_major_list_to_store(A_arr,rowLevel,colLevel,length_row)
        print("")
        print("Since this matrix is tiny one can print it out using")
        print("'mypy.print_naive(A_pID)' : ")
        
        mypy.print_naive(A_pID)
        print("")
        print("A 2 by 2 matrix can also be built with row_level = 1 and col_level = 1.")
        print("from the same input data [-1,0,1,2]:")
        print("")
        print("(The .py source code of this tutorial has the following 5 lines)")
        print("B_arr = list(map(str,[-1,0,1,2]))")
        print("rowLevel = 1")
        print("colLevel = 1")
        print("length_row = 1 << colLevel")
        print("B_pID = mypy.row_major_list_to_store(B_arr,rowLevel,colLevel,length_row)")
        print("")
        B_arr = list(map(str,[-1,0,1,2]))
        rowLevel = 1
        colLevel = 1
        length_row = 1 << colLevel
        B_pID = mypy.row_major_list_to_store(B_arr,rowLevel,colLevel,length_row)
        print("Since this matrix is tiny one can print it out using")
        print("'mypy.print_naive(B_pID)' :")
        mypy.print_naive(B_pID)
        
        
    print("")    
    print("**********************************************************************")
    print("")
    yes_no=input("Would you like to continue with examples of matrix file input and storage? (y/n): ")
    if yes_no =='y':
       # TESTING READING AND WRITING OF MATRICES
        filename_in_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr
        filename_out_naive = "Data/Out/sample.1.1.%s.naive" %scalarTypeStr
        filename_out_json = "Data/Out/sample.1.1.%s.json" %scalarTypeStr
        filename_out_nonzeros = "Data/Out/sample.1.1.%s.nonzeros" %scalarTypeStr
        print("")   
        print("LARC can write matrices using three file formats, usually indicated")
        print("by the suffix:  '.naive', '.nonzeros', and '.json'.")
        print("\t'.naive' to a standard matrix representation which lists all the")
        print("\t\telements of the matrix;")   
        print("\t'.nonzeros' refers to an output of just the nonzero elements; and")
        print("\t'.json' to a compressed LARC matrix format in a json container.")
        print("LARC can read matrices using three file formats:")
        print("\trow-major matrix format;")
        print("\tLARC compressed matrix format in json container; and")
        print("\tMatrix Market Exchange format.")
        print("")
        print("Example:")
        print("")
        print('For instance, if scalarType = "Real" has been set then one can read ')
        print(" a 2 by 2 matrix of [[1, 0],[.400000000000222, .200000000000111]] from ")
        print(" the file \"Data/In/sample.1.1.Real.rmm\".")
        print("")
        print("Code Examples:")
        print("")
        print("To read in the file \"%s\" in row-major matrix format," %os.path.join(os.path.dirname(__file__),filename_in_rmm))
        print("use the following two lines of python code:")
        print('   filename_in_rmm = "%s"' %filename_in_rmm)
        print("   sample_packedID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_in_rmm))")
        print("Now performing the reading in of file \"%s\"." %filename_in_rmm)
        sample_packedID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_in_rmm))       
        print("")
        print("")
        print("To write out the file \"%s\" in naive format," %os.path.join(os.path.dirname(__file__),filename_out_naive))
        print("use the following two lines of python code:")
        print('   filename_out_naive = "%s"' %filename_out_naive)
        print("   mypy.fprint_naive(sample_packedID,os.path.join(os.path.dirname(__file__),filename_out_naive))")
        print("Now performing the writing out of file \"%s\"." %filename_out_naive)
        print("(This prints every entry in the matrix so it should be used only for small matrices.)")
        mypy.fprint_naive(sample_packedID,os.path.join(os.path.dirname(__file__),filename_out_naive))
        print("")
        print("One can also print the sample matrix in naive format to the screen using ")
        print("   mypy.print_naive(sample_packedID)")
        print("")
        mypy.print_naive(sample_packedID)
        print("")
        print("To write out the file \"%s\" in nonzeros format," %os.path.join(os.path.dirname(__file__),filename_out_nonzeros))
        print("use the following two lines of python code:")
        print('   filename_out_nonzeros = "%s"' %filename_out_nonzeros)
        print("   mypy.fprint_matrix_nonzeros(sample_packedID,os.path.join(os.path.dirname(__file__),filename_out_nonzeros))")
        print("Now performing the writing out of file \"%s\"." %filename_out_nonzeros)
        mypy.fprint_matrix_nonzeros(sample_packedID,os.path.join(os.path.dirname(__file__),filename_out_nonzeros))
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
        print("to file \"%s\" , one can use the two python lines: " %filename_out_json)
        print('  filename_out_json = "%s"' %filename_out_json)
        print("  mypy.fprint_larcMatrixFile(sample_packedID,os.path.join(os.path.dirname(__file__),filename_out_json))")
        print("")        
        print("(Several examples of various sizes can be seen in the directory Data/In.)")
        mypy.fprint_larcMatrixFile(sample_packedID,os.path.join(os.path.dirname(__file__),filename_out_json))
        
        print("")        
        print('Using  os.system("more %s")' %filename_out_json,' One can see that the')
        print("LARC recursively compressed file looks like this:")
        os.system("more %s" %filename_out_json)
        
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
        
        
        
        #*###################################
        # read a MatrixMarket formated file #
        #*###################################
        
        


        #*#####################################################
        # Read a LARC matrix file in our own recursive format #
        #*#####################################################


        # print("*    3. LARC initialization from parameter files                     *")
        # print("*    4. recursive operations in LARC, matrix and operations stores   *")
        # print("*    5. writing your own recursive functions.                        *")  
 
