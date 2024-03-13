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
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *


##
# \file tut4_matrix_build_and_io.py
#
# \brief This code tests some basic matrix reading and building routines.
#
# 


if __name__ == '__main__':
	# This version references matrices by packedID instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    verbose = 3      # 0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL
    mypy.create_report_thread(1800)
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)


    #*############################################################
    #*  In the Makefile you can compile with:                   ##
    #*     TYPE=INTEGER, TYPE=REAL, TYPE=COMPLEX,               ##
    #*  or with multiprecision types:                           ##
    #*     TYPE=MPINTEGER, TYPE=MPRATIONAL,                     ##
    #*     TYPE=MPREAL, TYPE=MPCOMPLEX, or TYPE=MPRATCOMPLEX    ##
    #*############################################################

    #*  Define string for using in formating filenames
    scalarTypeStr = mypy.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\nDuring preload %d matrices have been created" %num_matrices_made)
    print("The contents of the matrix store can be output to a file")
    print("using the fprint_store_info_for_matrixID_range function.")
    end = num_matrices_made - 1
    filename = "Data/Out/preload.%s.store" %scalarTypeStr
    message_string = "After preload with parameters: %d, %d, %d" %(mat_store_exp,
                                                                   op_store_exp,
                                                                   max_level)
    mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename), message_string)
    # "After preload with parameters: 26, 24, 10.")

    #*#############################################
    # build array in C from Python list of scalars
    #*#############################################
    print("\nUsing row_major_list_to_store on data entered from python\n")

    # create entries for matrix
    if scalarTypeStr in ('Integer', 'MPInteger', 'Clifford'):
        a_str = [1, 3, 5, 6,8, 6, 3, 1,-9, 11, 13, 15,16, 13, 12, 10]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a_str = [1+2j, 3+4j, 5+6j, 7+8j,
                 8+7j, 6+5j, 3+4j, 1+2j,
                 9+10j, 11+12j, 13+14j, 15+16j,
                 16+15j, 14+13j, 12+11j, 10+9j]
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational'):
        a_str = [1, 3, 5, 6,8, 6, 3, 1,-9, 11, 13, 15,16, 13, 12, 10]
    elif scalarTypeStr in ('Boolean', 'Upper', 'Lower'):
        a_str = [1, 0, 1, 1,0, 0, 0, 1, 1,  1,  0,  1, 0,  0,  1,  1]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    #a_arr = list(map(str,a_str))
    a_arr = mypy.map_to_str(a_str,scalarTypeStr)
    
    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    a_ID = mypy.row_major_list_to_store(a_arr, level, level, dim_whole)
    mypy.print_naive(a_ID)
    print("\n")

    # Make a parent matrix from four copies of the a_ID matrix
    print("Creating matrix from get_pID_from_four_sub_pIDs on panel input and writing LARCMatrix file\n")
    panel = [a_ID]*4   # alternatively panel=[a_ID,a_ID,a_ID,a_ID]
    a_ID_parent = mypy.get_pID_from_four_sub_pIDs(a_ID,a_ID,a_ID,a_ID,3,3)
    mypy.print_naive(a_ID_parent)
    filename = "Data/Out/testfile.%s.json" %scalarTypeStr
    mypy.fprint_larcMatrixFile(a_ID_parent,os.path.join(os.path.dirname(__file__), filename))
    
    #  PLAYING WITH PREWRITTEN REVERSIBLE NAND LARCMatrix FILE
    print("About to test read LARCMatrix file\n")
    filename = "Data/In/nand.%s.json" %scalarTypeStr
    nand_packedID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix nand file\n")
    mypy.print_naive(nand_packedID)
    print("\n")

    # TESTING READING AND WRITING OF MATRICES
    print("Testing reading row major matrix format and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr
    filename_naive = "Data/Out/sample.1.1.%s.naive" %scalarTypeStr
    filename_json = "Data/Out/sample.1.1.%s.json" %scalarTypeStr

    print(os.path.join(os.path.dirname(__file__),filename_rmm))
    
    sample_packedID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    mypy.fprint_naive(sample_packedID,os.path.join(os.path.dirname(__file__),filename_naive))
    mypy.fprint_larcMatrixFile(sample_packedID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the rrm sample matrix in naive format to screen\n")
    mypy.print_naive(sample_packedID)

    print("Testing reading row major nonsquare matrices and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "Data/In/sample.1.2.%s.rmm" %scalarTypeStr
    filename_naive = "Data/Out/sample.1.2.%s.naive" %scalarTypeStr
    filename_json = "Data/Out/sample.1.2.%s.json" %scalarTypeStr

    sample_packedID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm))
    mypy.fprint_naive(sample_packedID,os.path.join(os.path.dirname(__file__),filename_naive))
    mypy.fprint_larcMatrixFile(sample_packedID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the nonsquare rrm sample matrix\n")
    mypy.print_naive(sample_packedID)

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_packedID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))

    print("We read in the LARCMatrix file\n")
    mypy.print_naive(samp_packedID)
    print("\n")

    print("Testing printing nonzeros to file.\n")
    filename_rmm = "Data/In/sample.1.3.%s.rmm" %scalarTypeStr
    filename_nonzeros = "Data/Out/sample.1.3.%s.nonzeros" %scalarTypeStr
    filename_json = "Data/Out/sample.1.3.%s.json" %scalarTypeStr
    
    sample_packedID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__), filename_rmm))
    mypy.fprint_matrix_nonzeros(sample_packedID,os.path.join(os.path.dirname(__file__),filename_nonzeros))
    mypy.fprint_larcMatrixFile(sample_packedID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Here is the matrix we are testing for printing out nonzero values")
    mypy.print_naive(sample_packedID)
    print("\n")

    # give example of printing the list of unique scalars
    print("\n*******************************************************")
    print("\nOne can print out a text file containing a list (one per line) of all the")
    print("unique scalars for a particular matrix. For example the scalars ...\n")
    filename = "Data/Out/testfile.%s.scalars" %scalarTypeStr
    our_path = os.path.join(os.path.dirname(__file__), filename)
    print("Printing scalar file of matrix with ID %d to file %s" %(a_ID_parent,our_path))
    mypy.fprint_uniqScalar_file(a_ID_parent,our_path)
    print("*******************************************************\n")


    print("Testing reading in Matrix Market Exchange format.\n")
    filename_mm = "Data/In/matrixMarketExchange1"
    sample_packedID = mypy.read_matrixMarketExchange_file(os.path.join(os.path.dirname(__file__), filename_mm))
    print("Here is the matrix we are testing for Matrix Market Exchange format")
    mypy.print_naive(sample_packedID)
    print("There are more examples of Matrix Market Exchange format files in the Count_trianges directory.")
    
    print("\n")

    # make CNOT
    print("\nHere is the CNOT (reversible XOR) matrix\n")

    #  CNOT_arr = mypy.buildArray([1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0])
    CNOT_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0]))
    CNOT_packedID = mypy.row_major_list_to_store(CNOT_arr,level,level,dim_whole)
    mypy.print_naive(CNOT_packedID)

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    start = end + 1
    end = num_matrices_made - 1
    filename = "Data/Out/cnot.%s.store" %scalarTypeStr
    mypy.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")

    # build Zero matrices
    print("\nHere is the level 2 zero matrix\n")
    # Z2_arr = mypy.buildArray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    Z2_arr = list(map(str,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
    Z2_packedID = mypy.row_major_list_to_store(Z2_arr,level,level,dim_whole)
    mypy.print_naive(Z2_packedID)

    # build Identity matrices
    print("\nHere is the level 2 identity matrix\n")
    # I2_arr = mypy.buildArray([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1])
    I2_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]))
    I2_packedID = mypy.row_major_list_to_store(I2_arr,level,level,dim_whole)
    mypy.print_naive(I2_packedID)

    # build a doubly-controlled NOT (base unit for reversible computing)
    print("\nThe 3 bit reversible AND (Toffoli) matrix with target 3rd.\n")
    # get_pID_from_four_sub_pIDs is under construction
    TOFFOLI_packedID= mypy.get_pID_from_four_sub_pIDs(I2_packedID,Z2_packedID,Z2_packedID,CNOT_packedID,3,3)
    mypy.print_naive(TOFFOLI_packedID)
    filename = "Data/Out/toffoli.%s.naive" %scalarTypeStr
    mypy.fprint_naive(TOFFOLI_packedID,os.path.join(os.path.dirname(__file__),filename))

    # use CCNOT to build a NAND
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
    # build 8x8 matrix I tensor I tensor X
    NOT_packedID = mypy.cvar.packedID_NOT;
    I1_packedID = mypy.get_identity_pID(1);
    not3_packedID = mypy.kronecker_product(I1_packedID,mypy.kronecker_product(I1_packedID, NOT_packedID));
    nand_from_Toff_packedID = mypy.matrix_mult(not3_packedID,TOFFOLI_packedID);
    mypy.print_naive(nand_from_Toff_packedID)
    filename = "Data/Out/nandfromtoff.%s.naive" %scalarTypeStr
    mypy.fprint_naive(nand_from_Toff_packedID,os.path.join(os.path.dirname(__file__),filename))
    # build from panels
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd, constructed a different way.\n")
    Z1_packedID = mypy.get_zero_pID(1,1);
    topleft_packedID = mypy.get_pID_from_four_sub_pIDs(NOT_packedID,Z1_packedID,Z1_packedID,NOT_packedID,2,2);
    botright_packedID = mypy.get_pID_from_four_sub_pIDs(NOT_packedID,Z1_packedID,Z1_packedID,I1_packedID,2,2);
    nand_from_panels_packedID = mypy.get_pID_from_four_sub_pIDs(topleft_packedID,Z2_packedID,Z2_packedID,botright_packedID,3,3);
    mypy.print_naive(nand_from_panels_packedID)
    filename = "Data/Out/nandfrompanels.%s.naive" %scalarTypeStr
    mypy.fprint_naive(nand_from_panels_packedID,os.path.join(os.path.dirname(__file__),filename))
    # compare with other nands
    if (nand_from_Toff_packedID != nand_packedID):
        print("problem with building NAND from Toffoli")
    if (nand_from_panels_packedID != nand_packedID):
        print("problem with building NAND from panels")


    # read the (1,2) entry from the nand matrix        
    scalarID = mypy.get_scalarID_from_pID_and_coords(nand_packedID, 1, 2);
    print("matrix ID of scalar in position (1,2) of NAND is %d"
        %mypy.matrixID_from_packedID(scalarID))
    print("The value of scalar in position (1,2) of NAND is")
    mypy.print_naive(scalarID)

    # read the (0,1) entry from the nand matrix        
    scalarID = mypy.get_scalarID_from_pID_and_coords(nand_packedID, 0, 1);
    print("matrix ID of scalar in position (0,1) of NAND is %d"
        %mypy.matrixID_from_packedID(scalarID))
    print("The value of scalar in position (0,1) of NAND is")
    mypy.print_naive(scalarID)
    print("The same value accessed directly using indices and get_readableString_scalar_from_pID_and_coords")
    myStr = mypy.get_readableString_scalar_from_pID_and_coords(nand_packedID,0,1);
    print("%s" %myStr);

    print("\nOne can use a commandline to check whether two LARC") 
    print("compressed matrices represent the same matrix, by using:")
    print("python ../larc/src/python/larcMatrix_files_matrix_equal.py  file1 file2")
