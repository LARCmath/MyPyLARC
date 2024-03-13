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
import random # for use in toeplitz
import numpy as np
from ctypes import *


## 
# \file test_hashstats.py
#
# \brief This program demonstrates the hash statistics output.
#
#  This lets the user see which statistics are collected when the LARC
#  src/larc.h has been modified to have the \#define HASHSTATS uncommented.
#  The program repeats various matrix operations in LARC and then gives
#  the report on the hashchain lengths and other statistics.
#

if __name__ == '__main__':
	# This version references matrices by packedIDs instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 14
    op_store_exp = 14
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)

    #*   START OF CODE STOLEN FROM test_cleaning.py

    # Define string for using in formating filenames
    scalarTypeStr = mypy.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1
    filename = "Data/Out/preload.%s.store" %scalarTypeStr
    mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"After preload with parameters: 26, 24, 10.")

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_pID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))

    print("We read in the LARCMatrix file\n")
    mypy.print_naive(samp_pID)
    print("\n")

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    if scalarTypeStr in ('Integer', 'MPInteger'):
        a_str = [1, 3, 5, 6,
                 8, 6, 3, 1,
                 -9, 11, 13, 15,
                 16, 13, 12, 10]
    elif scalarTypeStr == 'Boolean':
        a_str = [1, 0, 0, 0,
                 0, 0, 0, 1,
                 0, 1, 1, 1,
                 1, 1, 1, 0]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a_str = [1+2j, 3+4j, 5+6j, 7+8j,
                 8+7j, 6+5j, 3+4j, 1+2j,
                 9+10j, 11+12j, 13+14j, 15+16j,
                 16+15j, 14+13j, 12+11j, 10+9j]
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational', 'Clifford'):
        a_str = [1, 3, 5, 6,
                 8, 6, 3, 1,
                 -9, 11, 13, 15,
                 16, 13, 12, 10]
    elif scalarTypeStr in ('Upper', 'Lower'):
        a_str = [0.1, 0.3, 0.5, 0.6,
                 0.8, 0.6, 0.3, 0.1,
                 0, 1, 0.013, 0.15,
                 0.16, 0.13, 0.12, 0.10]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    arr = mypy.map_to_str(a_str, scalarTypeStr)
    print('arr:', arr)

    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    pID = mypy.row_major_list_to_store(arr, level, level, dim_whole)
    mypy.print_naive(pID)
    print("\n")

    # make a parent matrix from four copies of the packedID matrix
    print("Creating matrix from get_pID_from_four_sub_pIDs on panel input and writing LARCMatrix file\n")
    panel = [pID]*4   # alternatively panel=[pID,pID,pID,pID]
    pID_parent = mypy.get_pID_from_four_sub_pIDs(pID,pID,pID,pID,3,3)
    mypy.print_naive(pID_parent)
    filename = "Data/Out/testfile.%s.json" %scalarTypeStr
    mypy.fprint_larcMatrixFile(pID_parent,os.path.join(os.path.dirname(__file__), filename))

    
    #  PLAYING WITH PREWRITTEN REVERSIBLE NAND LARCMatrix FILE
    print("About to test read LARCMatrix file\n")
    filename = "Data/In/nand.%s.json" %scalarTypeStr
    nand_pID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix NAND file\n")
    mypy.print_naive(nand_pID)
    print("\n")

    # TESTING READING AND WRITING OF MATRICES
    print("Testing reading row major matrix format and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr
    filename_naive = "Data/Out/sample.1.1.%s.naive" %scalarTypeStr
    filename_json = "Data/Out/sample.1.1.%s.json" %scalarTypeStr
    
    sample_pID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    mypy.fprint_naive(sample_pID,os.path.join(os.path.dirname(__file__),filename_naive))
    mypy.fprint_larcMatrixFile(sample_pID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the rrm sample matrix in naive format to screen\n")
    mypy.print_naive(sample_pID)


    print("Testing reading row major nonsquare matrices and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "Data/In/sample.1.2.%s.rmm" %scalarTypeStr
    filename_naive = "Data/Out/sample.1.2.%s.naive" %scalarTypeStr
    filename_json = "Data/Out/sample.1.2.%s.json" %scalarTypeStr
    
    sample_pID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    mypy.fprint_naive(sample_pID,os.path.join(os.path.dirname(__file__),filename_naive))
    mypy.fprint_larcMatrixFile(sample_pID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the nonsquare rrm sample matrix\n")
    mypy.print_naive(sample_pID)

    print("Testing printing nonzeros to file.\n")
    filename_rmm = "Data/In/sample.1.3.%s.rmm" %scalarTypeStr
    filename_nonzeros = "Data/Out/sample.1.3.%s.nonzeros" %scalarTypeStr
    filename_json = "Data/Out/sample.1.3.%s.json" %scalarTypeStr
    
    sample_pID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__), filename_rmm))
    mypy.fprint_matrix_nonzeros(sample_pID,os.path.join(os.path.dirname(__file__),filename_nonzeros))
    mypy.fprint_larcMatrixFile(sample_pID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Here is the matrix we are testing for printing out nonzero values")
    mypy.print_naive(sample_pID)
    
    print("\n")


    # make CNOT
    print("\nHere is the CNOT (reversible XOR) matrix\n")
    CNOT_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0]))
    CNOT_pID = mypy.row_major_list_to_store(CNOT_arr,level,level,dim_whole)
    mypy.print_naive(CNOT_pID)

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    start = end + 1
    end = num_matrices_made - 1
    filename = "Data/Out/cnot.%s.store" %scalarTypeStr
    mypy.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")

    # build Zero matrices
    print("\nHere is the level 2 zero matrix\n")
    Z2_arr = list(map(str,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
    Z2_pID = mypy.row_major_list_to_store(Z2_arr,level,level,dim_whole)
    mypy.print_naive(Z2_pID)

    # build Identity matrices
    print("\nHere is the level 2 identity matrix\n")
    I2_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]))
    I2_pID = mypy.row_major_list_to_store(I2_arr,level,level,dim_whole)
    mypy.print_naive(I2_pID)

    # build a doubly-controlled NOT (base unit for reversible computing)
    print("\nThe 3 bit reversible AND (Toffoli) matrix with target 3rd.\n")
    TOFFOLI_pID= mypy.get_pID_from_four_sub_pIDs(I2_pID,Z2_pID,Z2_pID,CNOT_pID,3,3)
    mypy.print_naive(TOFFOLI_pID)
    filename = "Data/Out/toffoli.%s.naive" %scalarTypeStr
    mypy.fprint_naive(TOFFOLI_pID,os.path.join(os.path.dirname(__file__),filename))

    # use CCNOT to build an NAND
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
    NOT_packedID = mypy.cvar.packedID_NOT;
    I1_packedID = mypy.get_identity_pID(1);
    not3_packedID = mypy.kronecker_product(I1_packedID,mypy.kronecker_product(I1_packedID, NOT_packedID));
    nand_from_Toff_packedID = mypy.matrix_mult(not3_packedID,TOFFOLI_pID);
    mypy.print_naive(nand_from_Toff_packedID)
    filename = "Data/Out/nandfromtoff.%s.naive" %scalarTypeStr
    mypy.fprint_naive(nand_from_Toff_packedID,os.path.join(os.path.dirname(__file__),filename))


    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_pID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the json file\n")
    mypy.print_naive(samp_pID)
    
    print("does scalarM1_val print?")
    scalarM1_val = '-1'
    scalarM1_pID = mypy.get_valID_from_valString(scalarM1_val)
    mypy.print_naive(scalarM1_pID)
    
    print("testing scalar_mult:")
    samp2_pID = mypy.scalar_mult(scalarM1_pID,samp_pID)
    mypy.print_naive(samp2_pID)
    
    print("testing addition:")
    samp3_pID = mypy.matrix_add(samp_pID,samp2_pID)
    mypy.print_naive(samp3_pID)
    
    # save input packedIDs for testing op store hash chains later
    in1_test_sum_pID = samp_pID
    in2_test_sum_pID = samp2_pID
    
    print("testing adjoint:")
    samp3_pID = mypy.adjoint(samp_pID)
    mypy.print_naive(samp3_pID)
    adj_pID = samp3_pID
    
    print("testing non-square matrix mult:")
    samp4_pID = mypy.matrix_mult(samp3_pID,samp_pID)
    mypy.print_naive(samp4_pID)
    print("")
    samp4_pID = mypy.matrix_mult(samp_pID,samp3_pID)
    mypy.print_naive(samp4_pID)
    print("testing kron product:")
    samp4_pID = mypy.kronecker_product(samp_pID,samp_pID)
    mypy.print_naive(samp4_pID)
    

    print("testing join:")
    samp4_pID = mypy.join(samp_pID,samp_pID)
    mypy.print_naive(samp4_pID)
    print("testing stack:")
    samp4_pID = mypy.stack(samp_pID,samp_pID)
    mypy.print_naive(samp4_pID)


    #*  TESTING DELETION
    print("\nPreparing to delete a matrix from the store.\n")
    filename = "Data/Out/temp.%s.json" %scalarTypeStr
    mypy.fprint_larcMatrixFile(nand_pID, os.path.join(os.path.dirname(__file__),filename))
    mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    print("Previous range printed ended with matrixID %d\n" %end)
    if (end == num_matrices_made-1) :
        print("Nothing new since last matrix store print\n")
    else :
        start = end + 1
        end = num_matrices_made - 1
        filename = "Data/Out/nand.%s.store" %scalarTypeStr
        mypy.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__),filename),"Loaded NAND")

    # get the hashID and print the hash chain corresponding to a matrix we are about to delete
    hashID = mypy.hash_pID(nand_pID)
    comment = "hash chain before removal"
    filename = "Data/Out/hashChain.beforeMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    mypy.fprint_nonscalar_hash_chain_info(hashID, out_path, comment)

    
    # Test deletion of a matrix
    print("Testing removal of matrix from the matrix store\n")
    num_matrices_made =  mypy.num_matrices_created()
    end = num_matrices_made - 1
    filename = "Data/Out/nandYES.%s.store" %scalarTypeStr
    mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"Before Removed NAND")

    mypy.remove_matrix_from_store(nand_pID)
    filename = "Data/Out/nandNO.%s.store" %scalarTypeStr
    filename_json = "Data/Out/temp.%s.json" %scalarTypeStr
    print("\nDeleting the NAND matrix with matrixID", mypy.matrixID_from_packedID(nand_pID),"from store, which had been read from %s\n"  %filename_json)
    mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND")

    comment = "hash chain after removal"
    filename = "Data/Out/hashChain.afterMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    mypy.fprint_nonscalar_hash_chain_info(hashID, out_path, comment)


    mypy.list_op_names()
    
    # Test op store hash chains after deletion
	
    # print ops store report before deletion of a matrix
    mypy.op_store_report("stdout")
    
    # get hash for a operation record and print the hash chain
    op_type = mypy.get_op_type_from_string_name("SUM")
    sum_hashID = mypy.hash_from_op(in1_test_sum_pID,in2_test_sum_pID,op_type)
    if (sum_hashID == -1):
        print("invalid matrixID requested for op hash chain")
    else:
        mypy.op_hash_chain_info_to_screen(sum_hashID, "op hash chain before deletion")
        
    
    # delete the first input matrix
    mypy.remove_matrix_from_store(in1_test_sum_pID)
    print("deleting matrix with matrixID", mypy.matrixID_from_packedID(in1_test_sum_pID))
    
    # traverse the op_hash_chain by trying some new sums
    # test1_pID = mypy.matrix_add(samp4_pID,samp4_pID)
    # test2_pID = mypy.matrix_add(test1_pID,samp4_pID)
    # test3_pID = mypy.matrix_add(test2_pID,test1_pID)
    
    # set a hold on a matrix by packedID to see if it is immune to cleaning
    print("The matrixID of matrix to be held is", mypy.matrixID_from_packedID(adj_pID))
    mypy.set_hold_matrix(adj_pID)
    
    # clean the matrix store and print it again
    mypy.clean_matrix_storage()
    filename = "Data/Out/nandNOcleaned.%s.store" %scalarTypeStr
    mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND and cleaned matrix store")

    # clean the op store 
    for hash in range(1<<op_store_exp):
        mypy.clean_op_hash_chain(hash) 

    # print ops store report after deletion of a matrix
    mypy.op_store_report("stdout")

    # print same op hash chain again after deleting a matrix, holding a matrix and cleaning
    if (sum_hashID != -1):
        mypy.op_hash_chain_info_to_screen(sum_hashID, "op hash chain after deletion and cleaning")
    else:
        print("invalid matrixID requested for op hash chain")

    #*   END OF CODE STOLEN FROM test_cleaning.py


    #*  START OF CODE STOLEN FROM toeplitz.py

    # build array in C from Python list of scalars
    # print("Using row_major_list_to_store on data entered from python\n")

    # parameters for entering the python array into the store
    level = 8
    dim_whole = 2**level

    if scalarTypeStr in ('Real', 'MPReal', 'MPRational', 'Clifford', 'Upper', 'Lower'):
        randVals = [ random.random() for i in range(2*dim_whole-1)]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        randVals = [ np.complex(random.random(),random.random()) for i in range(2*dim_whole-1)]
    elif scalarTypeStr in ('Integer', 'MPInteger'):
        randVals = [ random.randrange(0,10001) for i in range(2*dim_whole-1)]
    elif scalarTypeStr == 'Boolean':
        randVals = [ random.randrange(0,2) for i in range(2*dim_whole-1)]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # create toeplitz matrix from the 2*dim_whole-1 random numbers
    a = []
    b = randVals
    # print("b: ", b)
    for i in range(dim_whole):
        a.append(list(b[dim_whole-i-1:2*dim_whole-i-1]))
    # print("a: ", a)
    amat = np.matrix(a)
    serial = mypy.add_numpy_matrix_to_matrix_store(amat)
    # filename_json = "Data/Out/toeplitz.lev%d.%s.json" %(level,scalarTypeStr)
    # mypy.print_naive(serial)
    # print("\n")
    # mypy.fprint_larcMatrixFile(serial,os.path.join(os.path.dirname(__file__),filename_json))



    #*  END OF CODE STOLEN FROM toeplitz.py


    #* The following lines only work if the C code has been compiled with
    #* #define HASHSTATS
    print("\nWARNING The last portion of test_hashstats.c only reports if ")
    print("the C code has been combiled with the #define HASHSTATS, otherwise")
    print("there is an error that there is no attribute matrix_hashstats.\n")

    filename_a = os.path.join(os.path.dirname(__file__),"Data/Out/hashaccesses.mat")
    filename_n = os.path.join(os.path.dirname(__file__),"Data/Out/hashnodes.mat")
    filename_r = os.path.join(os.path.dirname(__file__),"Data/Out/hashreport.mat")
    mypy.matrix_hashstats(filename_a,filename_n,filename_r)

    filename_a = os.path.join(os.path.dirname(__file__),"Data/Out/hashaccesses.op")
    filename_n = os.path.join(os.path.dirname(__file__),"Data/Out/hashnodes.op")
    filename_r = os.path.join(os.path.dirname(__file__),"Data/Out/hashreport.op")
    mypy.op_hashstats(filename_a,filename_n,filename_r)

