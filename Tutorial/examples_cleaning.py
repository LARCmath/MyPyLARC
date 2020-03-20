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

import os
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixIDs instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 15
    op_store_exp = 5
    max_level = 10
    rnd_sig_bits = -1   # default value
    trunc_to_zero_bits = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)

    scalarTypeStr = mypy.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1
    filename = "Data/Out/preload.%s.store" %scalarTypeStr
    iS_string = "After preload with parameters: " + str(mat_store_exp) + ", ";
    iS_string = iS_string + str(op_store_exp) + ", " + str(max_level) + "."
    mypy.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),iS_string)

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_mID = mypy.read_larcMatrix_file_return_matID(os.path.join(os.path.dirname(__file__),filename))

    print("We read in the LARCMatrix file\n")
    mypy.print_naive_by_matID(samp_mID)
    print("\n")

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    #############################
    # create a matrix in python #
    #############################
    if scalarTypeStr in ('Integer', 'MPInteger'):
        a_str = [1, 3, 5, 6,
                 8, 6, 3, 1,
                 -9, 11, 13, 15,
                 16, 13, 12, 10]
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a_str = [1+2j, 3+4j, 5+6j, 7+8j,
                 8+7j, 6+5j, 3+4j, 1+2j,
                 9+10j, 11+12j, 13+14j, 15+16j,
                 16+15j, 14+13j, 12+11j, 10+9j]
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational'):
        a_str = [1, 3, 5, 6,
                 8, 6, 3, 1,
                 -9, 11, 13, 15,
                 16, 13, 12, 10]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # turn the string into an array
    #a_arr = list(map(str,a_str))
    a_arr = mypy.map_to_str(a_str,scalarTypeStr)
    
    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    mID = mypy.row_major_list_to_store_matrixID(a_arr, level, level, dim_whole)
    mypy.print_naive_by_matID(mID)
    print("\n")

    # make a parent matrix from four copies of the matrixID matrix
    print("Creating matrix from get_matID_from_four_subMatIDs on panel input and writing LARCMatrix file\n")
    panel = [mID]*4   # alternatively panel=[mID,mID,mID,mID]
    mID_parent = mypy.get_matID_from_four_subMatIDs(mID,mID,mID,mID,3,3)
    mypy.print_naive_by_matID(mID_parent)
    filename = "Data/Out/testfile.%s.json" %scalarTypeStr
    mypy.write_larcMatrix_file_by_matID(mID_parent,os.path.join(os.path.dirname(__file__), filename))
    
    #  PLAYING WITH PREWRITTEN REVERSIBLE NAND JSON FILE
    print("About to test read LARCMatrix file\n")
    filename = "Data/In/nand.%s.json" %scalarTypeStr
    nand_mID = mypy.read_larcMatrix_file_return_matID(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix nand file\n")
    mypy.print_naive_by_matID(nand_mID)
    print("\n")

    # TESTING READING AND WRITING OF MATRICES
    print("Testing reading row major matrix format and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr
    filename_naive = "Data/Out/sample.1.1.%s.naive" %scalarTypeStr
    filename_json = "Data/Out/sample.1.1.%s.json" %scalarTypeStr
    
    sample_mID = mypy.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    mypy.write_naive_by_matID(sample_mID,os.path.join(os.path.dirname(__file__),filename_naive))
    mypy.write_larcMatrix_file_by_matID(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the rrm sample matrix in naive format to screen\n")
    mypy.print_naive_by_matID(sample_mID)

    print("Testing reading row major nonsquare matrices and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "Data/In/sample.1.2.%s.rmm" %scalarTypeStr
    filename_naive = "Data/Out/sample.1.2.%s.naive" %scalarTypeStr
    filename_json = "Data/Out/sample.1.2.%s.json" %scalarTypeStr
    
    sample_mID = mypy.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    mypy.write_naive_by_matID(sample_mID,os.path.join(os.path.dirname(__file__),filename_naive))
    mypy.write_larcMatrix_file_by_matID(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the nonsquare rrm sample matrix\n")
    mypy.print_naive_by_matID(sample_mID)

    print("Testing printing nonzeros to file.\n")
    filename_rmm = "Data/In/sample.1.3.%s.rmm" %scalarTypeStr
    filename_nonzeros = "Data/Out/sample.1.3.%s.nonzeros" %scalarTypeStr
    filename_json = "Data/Out/sample.1.3.%s.json" %scalarTypeStr
    
    sample_mID = mypy.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__), filename_rmm))
    mypy.write_matrix_nonzeros_by_matID(sample_mID,os.path.join(os.path.dirname(__file__),filename_nonzeros))
    mypy.write_larcMatrix_file_by_matID(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Here is the matrix we are testing for printing out nonzero values")
    mypy.print_naive_by_matID(sample_mID)
    
    print("\n")

    # make CNOT
    print("\nHere is the CNOT (reversible XOR) matrix\n")
    CNOT_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0]))
    CNOT_mID = mypy.row_major_list_to_store_matrixID(CNOT_arr,level,level,dim_whole)
    mypy.print_naive_by_matID(CNOT_mID)

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = mypy.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    start = end + 1
    end = num_matrices_made - 1
    filename = "Data/Out/cnot.%s.store" %scalarTypeStr
    mypy.matrix_store_info_to_file(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")

    # build Zero matrices
    print("\nHere is the level 2 zero matrix\n")
    Z2_arr = list(map(str,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
    Z2_mID = mypy.row_major_list_to_store_matrixID(Z2_arr,level,level,dim_whole)
    mypy.print_naive_by_matID(Z2_mID)

    # build Identity matrices
    print("\nHere is the level 2 identity matrix\n")
    I2_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]))
    I2_mID = mypy.row_major_list_to_store_matrixID(I2_arr,level,level,dim_whole)
    mypy.print_naive_by_matID(I2_mID)

    # build a Toffoli (base unit for reversible computing)
    print("\nThe 3 bit reversible AND (Toffoli) matrix with target 3rd.\n")
    TOFFOLI_mID= mypy.get_matID_from_four_subMatIDs(I2_mID,Z2_mID,Z2_mID,CNOT_mID,3,3)
    mypy.print_naive_by_matID(TOFFOLI_mID)
    filename = "Data/Out/toffoli.%s.naive" %scalarTypeStr
    mypy.write_naive_by_matID(TOFFOLI_mID,os.path.join(os.path.dirname(__file__),filename))

    # use Toffoli to build an NAND
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
    NOT_matrixID = mypy.cvar.matID_NOT;
    I1_matrixID = mypy.get_identity_matrixID(1);
    not3_matrixID = mypy.kronecker_product_matrixID(I1_matrixID,mypy.kronecker_product_matrixID(I1_matrixID, NOT_matrixID));
    nand_from_Toff_matrixID = mypy.matrix_mult_matrixID(not3_matrixID,TOFFOLI_mID);
    mypy.print_naive_by_matID(nand_from_Toff_matrixID)
    filename = "Data/Out/nandfromtoff.%s.naive" %scalarTypeStr
    mypy.write_naive_by_matID(nand_from_Toff_matrixID,os.path.join(os.path.dirname(__file__),filename))

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_mID = mypy.read_larcMatrix_file_return_matID(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix file\n")
    mypy.print_naive_by_matID(samp_mID)
    
    print("does scalarM1_val print?")
    scalarM1_val = '-1'
    scalarM1_mID = mypy.get_valID_from_valString(scalarM1_val)
    mypy.print_naive_by_matID(scalarM1_mID)
    
    print("testing scalar_mult:")
    samp2_mID = mypy.scalar_mult_matrixID(scalarM1_mID,samp_mID)
    mypy.print_naive_by_matID(samp2_mID)
    
    print("testing addition:")
    samp3_mID = mypy.matrix_add_matrixID(samp_mID,samp2_mID)
    mypy.print_naive_by_matID(samp3_mID)
    
    # save input matrixIDs for testing op store hash chains later
    in1_test_sum_mID = samp_mID
    in2_test_sum_mID = samp2_mID
    
    print("testing adjoint:")
    samp3_mID = mypy.matrix_adjoint_matrixID(samp_mID)
    mypy.print_naive_by_matID(samp3_mID)
    adj_mID = samp3_mID
    
    print("testing non-square matrix mult:")
    samp4_mID = mypy.matrix_mult_matrixID(samp_mID,samp3_mID)
    mypy.print_naive_by_matID(samp4_mID)
    # print("")
    # samp4_mID = mypy.matrix_mult_matrixID(samp3_mID,samp_mID)
    # mypy.print_naive_by_matID(samp4_mID)

    print("testing kron product:")
    samp4_mID = mypy.kronecker_product_matrixID(samp_mID,samp_mID)
    mypy.print_naive_by_matID(samp4_mID)
    

    print("testing join:")
    samp4_mID = mypy.join_matrixID(samp_mID,samp_mID)
    mypy.print_naive_by_matID(samp4_mID)
    print("testing stack:")
    samp4_mID = mypy.stack_matrixID(samp_mID,samp_mID)
    mypy.print_naive_by_matID(samp4_mID)


    ##  TESTING DELETION
    print("\nPreparing to delete a matrix from the store.\n")
    filename = "Data/Out/temp.%s.json" %scalarTypeStr
    mypy.write_larcMatrix_file_by_matID(nand_mID, os.path.join(os.path.dirname(__file__),filename))
    mypy.read_larcMatrix_file_return_matID(os.path.join(os.path.dirname(__file__),filename))

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
        mypy.matrix_store_info_to_file(start,end,os.path.join(os.path.dirname(__file__),filename),"Loaded NAND")

    # get the hashID and print the hash chain corresponding to a matrix we are about to delete
    hashID = mypy.matrix_hashID_from_matrixID(nand_mID)
    comment = "hash chain before removal"
    filename = "Data/Out/hashChain.beforeMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    # os.path.dirname(__file__) =  /.ccs/u01/jszito/LARC/tests/python
    # os.path.join("/.ccs/u01/jszito/LARC/tests/python","Data/Out/hashChain.beforeMatrixRemove")
    # out_path = "/.ccs/u01/jszito/LARC/tests/dat/out/hashChain.beforeMatrixRemove"
    mypy.matrix_hash_chain_info_to_file(hashID, out_path, comment)
    
    print("Now that we've loaded several matrices into the matrix store")
    print("(both as preloads and manually), we'll remove one of them.")
    user_input = input("Press any key to continue.")
    
    # Test deletion of a matrix
    print("\nTesting removal of matrix from the matrix store\n")
    num_matrices_made =  mypy.num_matrices_created()
    end = num_matrices_made - 1
    filename = "Data/Out/nandYES.%s.store" %scalarTypeStr
    mypy.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"Before Removed NAND")

    mypy.remove_matrix_from_mat_store_by_matrixID(nand_mID)
	
    filename = "Data/Out/nandNO.%s.store" %scalarTypeStr
    filename_json = "Data/Out/temp.%s.json" %scalarTypeStr
    print("\nDeleting the NAND matrix with matrixID", nand_mID,"from store, which had been read from %s\n"  %filename_json)
    mypy.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND")

    comment = "hash chain after removal"
    filename = "Data/Out/hashChain.afterMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    mypy.matrix_hash_chain_info_to_file(hashID, out_path, comment)

    print("We can also check the operations store and remove any operation")
    print("that has the removed matrix as one of its inputs.")
    user_input = input("Press any key to continue.")
    
    mypy.list_op_names()
    
    # Test op store hash chains after deletion
	
    # print ops store report before deletion of a matrix
    mypy.op_store_report("stdout")
    
    # print op store hash chain for "SUM"
    op_name = "SUM"
    
    # get hash for a operation record and print the hash chain
    sum_hashID = mypy.op_hashID_by_matrixIDs(in1_test_sum_mID,in2_test_sum_mID,op_name)
    if (sum_hashID == -1):
        print("invalid matrixID requested for op hash chain")
    else:
        mypy.op_hash_chain_info_to_screen(sum_hashID, "op hash chain before deletion")
        
    # delete the first input matrix
    mypy.remove_matrix_from_mat_store_by_matrixID(in1_test_sum_mID)
    print("deleting matrix with matrixID", in1_test_sum_mID)
    
    # traverse the op_hash_chain by trying some new sums
    # test1_mID = mypy.matrix_add_matrixID(samp4_mID,samp4_mID)
    # test2_mID = mypy.matrix_add_matrixID(test1_mID,samp4_mID)
    # test3_mID = mypy.matrix_add_matrixID(test2_mID,test1_mID)
    
    # set a hold on a matrix by matrixID to see if it is immune to cleaning
    print("The matrixID of matrix to be held is", adj_mID)
    mypy.set_hold_matrix_from_matrixID(adj_mID)
    
    # clean the matrix store and print it again
    mypy.clean_matrix_store()
    filename = "Data/Out/nandNOcleaned.%s.store" %scalarTypeStr
    mypy.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND and cleaned matrix store")

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
    
