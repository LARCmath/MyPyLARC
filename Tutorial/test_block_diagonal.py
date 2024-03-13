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
# \file test_block_diagonal.py
#
# \brief Demonstrates the list_block_diagonal_matrixID() subroutine.
#
#  This program creates a matrix and returns the matrixIDs of the block
#  diagonal submatrices of that matrix.  It does testing for various
#  error conditions.
#


## \brief Returns a list of matrixIDs for the diagonal blocks of a matrix
#  \param topMatID The MatrixID of the large matrix
#  \param blockLevel The level of the diagonal blocks
#  \param verbose Controls the information printed by the function
#  \return A list of matrixIDs, of size 2**(topMat_level-blockLevel)
#
def list_block_diagonal_matrixIDs(topMatID, blockLevel, verbose):
    # NOTE: verbose = 0 is run quiet
    #       verbose = 1 is warning level, for example let you know if
    #                   off diagonals are not zero
    #       verbose = 2 is print chatty information and warnings
    """
    This function returns the list of matrixIDs for the diagonal
    subblocks (with level blockLevel k) of a given matrix M (of level n).
    For more details see the comments in the recursive function
    recurse_list_block_diagonal_matrixIDs()
    """

    # Retrieve the level of topMatID and exception check
    topLevel = mypy.matrix_row_level(topMatID)
    if (topLevel != mypy.matrix_col_level(topMatID)):
        print("ERROR: in list_block_diagonal_matrixIDs:")
        print("\targument must be square matrix")
        return []

    if (topLevel<blockLevel):
        print("ERROR: in list_block_diagonal_matrixIDs:")
        print("\tblocklevel is greater than the level of the passed matrix")
        return []
    elif (topLevel==blockLevel):
        return [topMatID]

    # Create an array to contain all the matrixIDs for blocks
    array = [0]*(2**(topLevel-blockLevel))
    if (verbose>1):
        print("The length of the array is %d" %len(array))

    # Set initial parameters for recursive call
    startArrayIndex = 0
    arrayLengthExponent = topLevel-blockLevel
    currentLevel = topLevel
    array[0] = topMatID
    recurse_list_block_diagonal_matrixIDs(array,
                                          startArrayIndex,
                                          arrayLengthExponent,
                                          currentLevel,
                                          blockLevel,
                                          verbose)        
    return array
    

##
# \ingroup larc
# \brief This function is the recursive call for list_block_diagonal_matrixIDs
# \param array Eventually, a list of matrixIDs of diagonal subblocks
# \param startArrayIndex The start of the portion of the array worked on
# \param arrayLengthExponent The length of the portion of the array worked on
# \param currentLevel The size of the matrices in this part of the recursion
# \param blockLevel The size of the matrices for which we want the MatrixIDs
# \param verbose Determines how much info is printed to screen
#
def recurse_list_block_diagonal_matrixIDs(array,
                                          startArrayIndex,
                                          arrayLengthExponent,
                                          currentLevel,
                                          blockLevel,
                                          verbose):
    """
    This recursive function returns the list of matrixIDs for the diagonal
    subblocks (with level blockLevel==k) of a given matrix M (of level==n).
    The list "array" must pre-exist and be of size 2**(n-k) 
          with array[0] = matrixID of M.
    The parameters startArrayIndex, arrayLengthExponent, and currentLevel
    are modified before passing to the next recursive calls. Their initial
    values on first call to the routine must be:
        startArrayIndex: 0
        arrayLengthExponent: n-k (log base 2 of the size of array)
        currentLevel: n
    """
    # The matrixID of this given matrix will initially be passed as entry 0 in
    # the array. When the algorithm is finished array will contain the list of
    # matrixIDs of the blocks in the order from top left to bottom right.
    #     array[startArrayIndex] contains the matrixID that will be quartered.
    #     2**arrayLengthExponent is the length of the array currently being
    #       modified
    #     currentLevel is initially n and will be decremented during the
    #       recursion.
    #     blockLevel is the size of the blocks that we want and stays = k.
    # 
    # NOTE: verbose = 0 is run quiet
    #       verbose = 1 is warning level, for example let you know if
    #                   off diagonals are not zero
    #       verbose = 2 is print chatty information and warnings

    # Branch: divide the matrix into quadrants,
    #    then run recursively on the top left and bottom right
    # currentLevel will be greater than the blockLevel.
    if (verbose>1):
       print("the matrixID that we are quartering is %d" %(
           mypy.matrixID_from_packedID(array[startArrayIndex])))

    # Find matrixID for the topLeftQuadrant and bottomRightQuadrant of the current matrix
    mat0pID = mypy.get_pID_of_indexed_submatrix(array[startArrayIndex],0)
    mat3pID = mypy.get_pID_of_indexed_submatrix(array[startArrayIndex],3)

    # Check to see if off diagonals are zero matrices, if not WARN
    if verbose:
        mat1pID = mypy.get_pID_of_indexed_submatrix(array[startArrayIndex],1)
        mat2pID = mypy.get_pID_of_indexed_submatrix(array[startArrayIndex],2)
        print("matrix IDs for off diagonals are %d and %d" %(
           mypy.matrixID_from_packedID(mat1pID),
           mypy.matrixID_from_packedID(mat2pID),))
        if (mypy.matrix_is_zero(mat1pID)==0):
            print("WARNING: top right quadrant of current submatrix is nonzero")
            print("\tcurrentLevel is %d, startArrayIndex is %d, and arrayLengthExponent is %d"
                  %(currentLevel,startArrayIndex,arrayLengthExponent))
        if (mypy.matrix_is_zero(mat2pID)==0):
            print("WARNING: bottom left quadrant of current submatrix is nonzero")
            print("\tcurrentLevel is %d, startArrayIndex is %d, and arrayLengthExponent is %d"
                  %(currentLevel,startArrayIndex,arrayLengthExponent))
            
    # Load these matrixIDs into the array at the start and halfway through
    offsetArrayIndex = startArrayIndex + 2**(arrayLengthExponent-1)
    array[startArrayIndex] = mat0pID
    array[offsetArrayIndex] = mat3pID
    if (verbose>1):
        print("the top diagonal quadrant has matrixID %d" %(mypy.matrixID_from_packedID(array[startArrayIndex])))
        print("the bottom diagonal quadrant has matrixID %d" %(mypy.matrixID_from_packedID(array[offsetArrayIndex])))
        
    # check to see if we are at the bottom of the recursion
    if (currentLevel-1==blockLevel):
        if (verbose>1):
            print("At bottom level of recursion") 
        return
    
    # Recurse
    if (verbose>1):
        print("recursing on top diagonal quadrant with array index %d" %(startArrayIndex))
    recurse_list_block_diagonal_matrixIDs(array,
                                          startArrayIndex,
                                          arrayLengthExponent-1,
                                          currentLevel-1,
                                          blockLevel,
                                          verbose)       
    if (verbose>1):
        print("recursing on bottom diagonal quadrant with array index %d" %(offsetArrayIndex))
    recurse_list_block_diagonal_matrixIDs(array,
                                          offsetArrayIndex,
                                          arrayLengthExponent-1,
                                          currentLevel-1,
                                          blockLevel,
                                          verbose)              

if __name__ == '__main__':
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
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
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1
    filename = "Data/Out/preload.%s.store" %scalarTypeStr
    mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"After preload with parameters: 26, 24, 10.")

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    # create a matrix in python
    if scalarTypeStr in ('Integer', 'MPInteger'):
        a = np.matrix([[1, 3, 0, 0],
                       [8, 6, 0, 0],
                       [0, 0, 13, 15],
                       [0, 0, 12, 10]])
    elif scalarTypeStr == 'Boolean':
        a = np.matrix([[1, 1, 0, 0],
                       [1, 1, 0, 0],
                       [0, 0, 1, 1],
                       [0, 0, 1, 0]])
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a = np.matrix([[1+2j, 3+4j, 5+6j, 7+8j],
                       [8+7j, 6+5j, 3+4j, 1+2j],
                       [9+10j, 11+12j, 13+14j, 15+16j],
                       [16+15j, 14+13j, 12+11j, 10+9j]])
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational', 'Clifford'):
        a = np.matrix([[1, 3, .5, 6],
                       [8, 6, 3, .1],
                       [-9, 11, 13, 1.5],
                       [16, 13, 12, 10]])
    elif scalarTypeStr in ('Upper', 'Lower'):
        a = np.matrix([[0.1, 0.3, 0.5, 0.6],
                       [0.8, 0.6, 0.3, 0.1],
                       [0, 1, 0.013, 0.15],
                       [0.16, 0.13, 0.12, 0.10]])
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # test our new function matrix_is_zero(z_pID))
    if scalarTypeStr in ('Integer', 'MPInteger'):
        z = np.matrix([[0, 0],[0, 0]])
        z_pID = mypy.add_numpy_matrix_to_matrix_store(z)
        mypy.print_naive(z_pID)
        print("\nThe fucntion matrix_is_zero returns:")
        print(mypy.matrix_is_zero(z_pID))
        print("\n")
    
    # turn the matrix into an array by reading off each row in turn (row major format)
    serial = mypy.add_numpy_matrix_to_matrix_store(a)
    mypy.print_naive(serial)
    print("\n")

    level = 2
    dim_whole = 2**level

    # Make a parent matrix from four copies of the serial matrix
    print("Creating matrix from get_pID_from_four_sub_pIDs on panel input and writing LARCMatrix file\n")

    zero2 = mypy.get_zero_pID(2,2)
    
    # panel = [serial,zero2,zero2,serial]
    top_pID = mypy.get_pID_from_four_sub_pIDs(serial,zero2,zero2,serial,3,3)
    mypy.print_naive(top_pID)
    filename = "Data/Out/blocktest.%s.json" %scalarTypeStr
    mypy.fprint_larcMatrixFile(top_pID,os.path.join(os.path.dirname(__file__), filename))


    verbose = 1    # 0 run quiet, 1 warn if not block diagonal, 2 be chatty

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 4
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = list_block_diagonal_matrixIDs(top_pID,blockLevel,verbose)
    print(list(map(mypy.matrixID_from_packedID,block_array)))

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 3
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = list_block_diagonal_matrixIDs(top_pID,blockLevel,verbose)
    print(list(map(mypy.matrixID_from_packedID,block_array)))

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 2
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = list_block_diagonal_matrixIDs(top_pID,blockLevel,verbose)
    print(list(map(mypy.matrixID_from_packedID,block_array)))
    
    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 1
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = list_block_diagonal_matrixIDs(top_pID,blockLevel,verbose)
    print(list(map(mypy.matrixID_from_packedID,block_array)))

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 0
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = list_block_diagonal_matrixIDs(top_pID,blockLevel,verbose)
    print(list(map(mypy.matrixID_from_packedID,block_array)))


