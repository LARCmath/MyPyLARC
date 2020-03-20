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

 #Your_Project_template
 # Some example application "playgrounds" for eigen, fft, toeplitz
 # can be found under mypylarc

#Plan for checking whether A is near symmetric

#Calculate A from some internal LARC process
#Call transpose routine
#     B = transpose(A)
#Call matrix subtraction routine
#     C = (A,B,subtract)  or C = subtract(A,B)
#Call matrix_max_norm routine
#     d = matrix_max-norm(C)
#Convert LARC values for printing
#     e = convert_LARC2FLT d
#Print result
#     print("\n result is ",e,"\n")

 
 #Imports of utility python code
from __future__ import print_function

import os
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print("This code provides a template for using basic LARC routines \n")

    print("At this point in the code one specifies basic LARC parameters ")
    print("It would be good to have these as command line parameters ")
    print("and have them echo-ed out at the beginning of each run ")
    print("It also would be good to start with toy values during debugging ")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    rnd_sig_bits = 40   # (default value -1 =>62 bits )
    trunc_to_zero_bits = 20  # (default value -1 => 62 bits)
    verbose = 1



    mypy.create_report_thread(1800)
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)

    print("mat_store_exp = ",mat_store_exp)
    print("op_store_exp = ",op_store_exp)
    print("max_level = ",max_level)
    print("rnd_sig_bits = ",rnd_sig_bits,"; ( -1 implies a default 53 bits??)")
    print("trunc_to_zero_bits = ",trunc_to_zero_bits,"; ( -1 for  1074 bits??)")
    print("verbose = ",verbose)

    # scalarType is set at compile time (e.g. make TYPE=INTEGER) 
    # this command returns the scalar type in case you need to know it
    scalarTypeStr = mypy.cvar.scalarTypeStr;

    ## Calculate number of matrices created, then print part of matrix store
    #num_matrices_made = mypy.num_matrices_created()
    #print("\nDuring preload %d matrices have been created" %num_matrices_made)
    #print("The contents of the matrix store can be output to a file")
    #print("using the matrix_store_info_to_file function.")
    #end = num_matrices_made - 1
    #filename = "Data/Out/preload.%s.store" %scalarTypeStr
    #message_string = "After preload with parameters: %d, %d, %d" %(mat_store_exp,op_store_exp,max_level)
    #mypy.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename), message_string)


    #Create, Convert and Load auxiliary scalars
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C/Python subroutines:...
	     # ...

	     #SPECIFY IN C OR PYTHON or numpy
	     #const_pi = 3.141592653      #Real
	     #const_root2 = sqrt(2.)      #Real
	     #constant_eipi/4 = [1/const_root2,1/const_root2]  #Complex
	     #CONVERT AND LOAD INTO LARC
	     #...
	     #Or RECURSIVE GENERATION

    #Create, Convert and Load auxiliary row vectors
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C/Python subroutines:...
	     # ...
	     #SPECIFY IN C OR PYTHON or numpy
	     # List1 = [1., 2., 3., 4.,]
	     # V1Row = [1., 2., 3., 4.,]
	     #CONVERT AND LOAD INTO LARC
	     #...
	     # or RECURSIVE GENERATION

    #Create, Convert and Load auxiliary column vectors
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C/Python subroutines:...
	     # ...
	     # List1 = [1., 2., 3., 4.,]
	     # V1Col = [[1.], [2.], [3.], [4.],]
	     #...
	     # or RECURSIVE GENERATION

    #Create, Convert and load auxiliary matrices
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     #...
	     # List1 = [1., 2., 3., 4.,]
	     # M = [ [1. , 2.] , [3. , 4.] ]
	     # ...
	     #CONVERT AND LOAD INTO LARC
	     #...
	     # or RECURSIVE GENERATION

    ###############################################
    # build array in C from Python list of scalars
    ###############################################
    print("\nUsing row_major_list_to_store on data entered from python\n")

    # create entries for matrix
    if scalarTypeStr in ('Integer', 'MPInteger'):
        a_str = [1, 3, 5, 6,8, 6, 3, 1,-9, 11, 13, 15,16, 13, 12, 10]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a_str = [1+2j, 3+4j, 5+6j, 7+8j,
                 8+7j, 6+5j, 3+4j, 1+2j,
                 9+10j, 11+12j, 13+14j, 15+16j,
                 16+15j, 14+13j, 12+11j, 10+9j]
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational'):
        a_str = [1, 3, 5, 6,8, 6, 3, 1,-9, 11, 13, 15,16, 13, 12, 10]
    else:
        raise Exception('Do not know how to build matrix for type %s.'%(mypy.cvar.scalarTypeStr,))

    #a_arr = list(map(str,a_str))
    a_arr = mypy.map_to_str(a_str,scalarTypeStr)
    
    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    a_ID = mypy.row_major_list_to_store_matrixID(a_arr, level, level, dim_whole)

    ################################################################
    # Produce the naive print of a Matrix identified by its Mat_ID #
    ################################################################

    mypy.print_naive_by_matID(a_ID)
    print("\n")

    #Perform LARC Linear Operations (or LARC/PYTHON hybrid operations???)
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     # ...



    ###########################################################################
    # Return the Mat_ID of the transpose of a Matrix identified by its Mat_ID #
    ###########################################################################
    
    userC_arr_ID = mypy.matrix_adjoint_matrixID(a_ID)
    
    ################################################################
    # Produce the naive print of a Matrix identified by its Mat_ID #
    ################################################################
    mypy.print_naive_by_matID(userC_arr_ID)
    print("\n")

    ####################################################################
    # Return the Mat_ID of the diff of Matrices identified by Mat_ID s #
    ####################################################################
    userD_arr_ID = mypy.matrix_diff_matrixID(a_ID,userC_arr_ID)

    ################################################################
    # Produce the naive print of a Matrix identified by its Mat_ID #
    ################################################################
    mypy.print_naive_by_matID(userD_arr_ID)
    print("\n")

    ####################################################################
    # Produce the maximum element of a Matrix identified by its Mat_ID #
    ####################################################################

    maxabselement = mypy.matrix_maxnorm_matrixID(userD_arr_ID)

    print("max absolute element of diff of A and A_transpose ",maxabselement)




    #Evaluate Properties of LARC objects (norms, is_symmetric, is_toeplitz, ...
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     # ...

    #Output Results
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     # ...

print(" If ", maxabselement, " is very very small, then the matrix is near symmetric. Otherwise not")
  

