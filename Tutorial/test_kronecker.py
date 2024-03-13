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
# \file test_kronecker.py
#
# \brief This program creates kronecker products of various combinations of 
# matrices, row vectors, and column vectors.
#


if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print("This code tests the kronecker product routine")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)


    # In the Makefile you can compile with different scalarType values
    # Define string for using in formating filenames
    scalarTypeStr = mypy.cvar.scalarTypeStr

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    # create a matrix in python
    if scalarTypeStr in ('Integer', 'MPInteger'):
        a = np.matrix([[1, 2, 3, 4],[5, 6, 7, 8]])
        b = np.matrix([[9, 10, 11, 12]])
    elif scalarTypeStr == 'Boolean':
        a = np.matrix([[1, 0, 1, 0],[0, 1, 1, 1]])
        b = np.matrix([[0, 1, 1, 0]])
    elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
        a = np.matrix([[1+2j, 2+3j, 3+4j, 4+5j],[5+1j, 6+2j, 7+3j, 8+4j]])
        b = np.matrix([[9-1j, 10-2j, 11-3j, 12-4j]])
    elif scalarTypeStr in ('Real', 'MPReal', 'MPRational', 'Clifford'):
        a = np.matrix([[1, 2, 3, 4],[5, 6, 7, 8]])
        b = np.matrix([[9, 10, 11, 12]])
    elif scalarTypeStr in ('Upper', 'Lower'):
        a = np.matrix([[0.1, 0.2, 0.3, 0.4],[0.5, 0.6, 0.7, 0.8]])
        b = np.matrix([[0.9, 0.10, 0.11, 0.12]])
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    print("COLUMN_COLUMN CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 3, 0, 1)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 2, 0, 1)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product A \otimes B:")
    serial_c = mypy.kronecker_product(serial_a,serial_b)
    mypy.print_naive(serial_c)
    print("\n")

    print("kronecker product B \otimes A:")
    serial_c = mypy.kronecker_product(serial_b,serial_a)
    mypy.print_naive(serial_c)
    print("\n")

    print("ROW_COLUMN CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 0, 3, 8)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 2, 0, 1)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product A \otimes B:")
    serial_c = mypy.kronecker_product(serial_a,serial_b)
    mypy.print_naive(serial_c)
    print("\n")

    print("COLUMN_ROW CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 3, 0, 1)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 0, 2, 4)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product A \otimes B:")
    serial_c = mypy.kronecker_product(serial_a,serial_b)
    mypy.print_naive(serial_c)
    print("\n")

    print("ROW_ROW CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 0, 3, 8)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 0, 2, 4)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product A \otimes B:")
    serial_c = mypy.kronecker_product(serial_a,serial_b)
    mypy.print_naive(serial_c)
    print("\n")

    print("kronecker product B \otimes A:")
    serial_c = mypy.kronecker_product(serial_b,serial_a)
    mypy.print_naive(serial_c)
    print("\n")

    print("MATRIX_ROW CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 1, 2, 4)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 0, 2, 4)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product A \otimes B:")
    serial_c = mypy.kronecker_product(serial_a,serial_b)
    mypy.print_naive(serial_c)
    print("\n")

    print("MATRIX_COLUMN CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 1, 2, 4)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 2, 0, 1)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product A \otimes B:")
    serial_c = mypy.kronecker_product(serial_a,serial_b)
    mypy.print_naive(serial_c)
    print("\n")

    print("ROW_MATRIX CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 1, 2, 4)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 0, 2, 4)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product B \otimes A:")
    serial_c = mypy.kronecker_product(serial_b,serial_a)
    mypy.print_naive(serial_c)
    print("\n")

    print("COLUMN_MATRIX CASE")
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = mypy.map_to_str(alist, scalarTypeStr)
    # print('array A:', mypy.str_scalarTypeArray(aarr, len(alist)))
    print('array A:', aarr)

    # creating or finding the matrix associated with the array
    serial_a = mypy.row_major_list_to_store(aarr, 1, 2, 4)
    mypy.print_naive(serial_a)
    print("\n")

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = mypy.map_to_str(blist, scalarTypeStr)
    # print('array B:', mypy.str_scalarTypeArray(barr, len(blist)))
    print('array B:', barr)

    # creating or finding the matrix associated with the array
    serial_b = mypy.row_major_list_to_store(barr, 2, 0, 1)
    mypy.print_naive(serial_b)
    print("\n")

    print("kronecker product B \otimes A:")
    serial_c = mypy.kronecker_product(serial_b,serial_a)
    mypy.print_naive(serial_c)
    print("\n")
