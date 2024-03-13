#!/usr/bin/env python3

# NOTE: to run on command line
#   python3 -m unittest -v test_unittest_matmath
# To run individual test classes from the module, do (for example) 
#   python3 -m unittest -v test_unittest_matmath.TestMatrixMaxnorm

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

 # This test set looks at routines in matmath.c and compares
 # their output with an output calculated in python.

from __future__ import print_function

import os 
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
from ctypes import *
import random
import unittest
import math


##
# \file test_unittest_matmath.py
#
# \brief This is a collection of unit tests for the various matrix functions
# in larc/src/matmath.c.
#
# See also working_example_unittest.py


# Python routine to compute the complex conjugate of a number:
def conj(val):
    if type(val) == type(1j):
        return val.conjugate()
    elif type(val) == type(1.0):
        return val
    elif type(val) == type(1):
        return val
    else:
        assert False, "Unknown value type {0} of {1}.".format(type(val), val)


# Python routine to take matrix entry list and square each entry
# then multiply each entry by scale_factor
def squareEntries(entries, scale_factor):
    return [scale_factor*(conj(entry)*entry) for entry in entries]


# tests the matmath.c routine  matrix_entrySquared
class TestMatrixEntrySquared(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with mypy.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 2
        self.col_level = 2
        self.val_range= [-100, 100]
        self.sparsity = 0.3

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_entrysquared_random_zeroscale(self):
        """
        0 * entrysquared(A) = 0 should work for random A.
        """
        scale_factor = 0
        n = (2**self.row_level)*(2**self.col_level)

        vals = mypy.matrix_random_entries(self.scalarTypeStr, n, self.val_range[0], self.val_range[1], self.sparsity)
        squared_vals = [0 for val in vals]

        vals_pID = mypy.row_major_list_to_store(mypy.map_to_str(vals, self.scalarTypeStr), self.row_level, self.col_level, 2**self.col_level)
        squared_vals_pID = mypy.row_major_list_to_store(mypy.map_to_str(squared_vals, self.scalarTypeStr), self.row_level, self.col_level, 2**self.col_level)
        self.assertEqual(mypy.matrix_entrySquared(vals_pID, str(scale_factor)), squared_vals_pID)

    @unittest.skipIf(mypy.cvar.scalarTypeStr in ("Real", "Complex", "MPReal", "MPComplex"), "ScalarType must not be Real, Complex, MPReal or MPComplex")
    def test_entrysquared_random_noscale(self):
        """
        1 * entrysquared(A) should work for random A.
        """
        scale_factor = 1
        n = (2**self.row_level)*(2**self.col_level)

        vals = mypy.matrix_random_entries(self.scalarTypeStr, n, self.val_range[0], self.val_range[1], self.sparsity)
        squared_vals = squareEntries(vals, 1)

        vals_pID = mypy.row_major_list_to_store(mypy.map_to_str(vals, self.scalarTypeStr), self.row_level, self.col_level, 2**self.col_level)
        squared_vals_pID = mypy.row_major_list_to_store(mypy.map_to_str(squared_vals, self.scalarTypeStr), self.row_level, self.col_level, 2**self.col_level)
        squared_pID = mypy.matrix_entrySquared(vals_pID, "1") 
        failure_msg = "matrices displayed above"
        if squared_pID != squared_vals_pID:
            print("Failure in test_entrysquared_random_noscale - printing relevant matrices.")
            print("Original matrix:")
            mypy.print_naive(vals_pID)
            print("Squared (by python) matrix:")
            mypy.print_naive(squared_vals_pID)
            print("Squared (by LARC) matrix:")
            mypy.print_naive(squared_pID)
        self.assertEqual(squared_pID, squared_vals_pID, failure_msg)

    @unittest.skipIf(mypy.cvar.scalarTypeStr in ("Real", "Complex", "MPReal", "MPComplex"), "ScalarType must not be Real, Complex, MPReal or MPComplex")
    def test_entrysquared_random_scale(self):
        """
        s * entrysquared(A) should work for random A.
        """
        scale_factor = mypy.matrix_random_entry(self.scalarTypeStr, -10, 10)
        n = (2**self.row_level)*(2**self.col_level)

        vals = mypy.matrix_random_entries(self.scalarTypeStr, n, self.val_range[0], self.val_range[1], self.sparsity)
        squared_vals = squareEntries(vals, scale_factor)

        vals_pID = mypy.row_major_list_to_store(mypy.map_to_str(vals, self.scalarTypeStr), self.row_level, self.col_level, 2**self.col_level)
        squared_vals_pID = mypy.row_major_list_to_store(mypy.map_to_str(squared_vals, self.scalarTypeStr), self.row_level, self.col_level, 2**self.col_level)
        squared_pID = mypy.matrix_entrySquared(vals_pID, mypy.value_to_string(scale_factor, self.scalarTypeStr))
        failure_msg = "matrices displayed above"
        if squared_pID != squared_vals_pID:
            print("Failure in test_entrysquared_random_scale - printing relevant matrices.")
            print("Scale factor: {0} ({1})".format(scale_factor, scale_factor.hex()))
            print("Original matrix:")
            mypy.print_naive(vals_pID)
            print("Squared (by python) matrix:")
            mypy.print_naive(squared_vals_pID)
            print("Squared (by LARC) matrix:")
            mypy.print_naive(squared_pID)
        self.assertEqual(squared_pID, squared_vals_pID, failure_msg)


class TestMatrixMaxnorm(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with mypy.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = 0.3

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_maxnorm_zero(self):
        """
        maxnorm(0) = 0
        """
        #* grab Zero matrix
        zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        maxnormID = mypy.normID(zeroMatID,0)
        maxnorm = mypy.traceID(maxnormID)
        if self.scalarTypeStr in ("Complex", "MPComplex", "MPRatComplex"):
            maxnormList = list(map(float, maxnorm.split("+I*")))
            if (len(maxnormList)==2):
                self.assertEqual(maxnormList[1],0)
            maxnorm = maxnormList[0]
        else:
            maxnorm = float(maxnorm)
        self.assertEqual(maxnorm, 0)

    def test_maxnorm_id(self):
        """
        maxnorm(id) = 1
        """
        idMatID = mypy.get_identity_pID(self.row_level)
        maxnormID = mypy.normID(idMatID,0)
        maxnorm = mypy.traceID(maxnormID)
        if self.scalarTypeStr in ("Complex", "MPComplex", "MPRatComplex"):
            maxnormList = list(map(float, maxnorm.split("+I*")))
            if (len(maxnormList)==2):
                self.assertEqual(maxnormList[1],0.0)
            maxnorm = maxnormList[0]
        else:
            maxnorm = float(maxnorm)
        self.assertEqual(maxnorm, 1)

    def test_maxnorm_negated(self):
        """
        maxnorm(-A) = maxnorm(A)
        """
        # build -1
        neg1 = mypy.get_valID_from_valString("-1")
        # generate random matrix
        randMatAID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        negAMatID = mypy.scalar_mult(neg1, randMatAID)
        self.assertEqual(mypy.normID(randMatAID,0), mypy.normID(negAMatID,0))

    def test_maxnorm_random(self):
        """
        placing != entry 5 in zero matrix to make A
        maxnorm(A) = 5
        """
        # get a zero matrix and place a 5 some where in the matrix
        zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        # TO BE VERY CONFUSING we need an (i,j) index in the matrix we use:
        # for i: we know row_level <= 2^row_level so the index i=self.row_level is a valid row index
        # for j: we know col_level <= 2^col_level so the index j=self.col_level is a valid col index
        matAID = mypy.replace_scalar_in_matrix_by_string_and_coords(zeroMatID, self.row_level, self.col_level, "5");
        maxnormID = mypy.normID(matAID,0)
        maxnorm = mypy.traceID(maxnormID)
        if self.scalarTypeStr in ("Complex", "MPRatComplex", "MPComplex"):
#            maxnorm = sum(map(float, maxnorm.split("+I*")))
            maxnormList = list(map(float, maxnorm.split("+I*")))
            if (len(maxnormList)==2):
                self.assertEqual(maxnormList[1],0.0)
            maxnorm = maxnormList[0]
        else:
            maxnorm = float(maxnorm)
        self.assertEqual(maxnorm, 5)


class TestMatrixL2norm(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = 0.3

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_L2norm(self):
        """
        l2Norm(0) = 0
        """
        #* grab Zero matrix
        zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
#        maxnormID = normID(zeroMatID,0)
#        maxnorm = traceID(maxnormID)
        l2NormID = mypy.normID(zeroMatID,2)
        l2Norm = mypy.traceID(l2NormID)
        if self.scalarTypeStr in ("Complex", "MPRatComplex", "MPComplex"):
#            l2Norm = sum(map(float, l2Norm.split("+I*")))
            l2NormList = list(map(float, l2Norm.split("+I*")))
            if (len(l2NormList)==2):
                self.assertEqual(l2NormList[1],0)
            l2Norm = l2NormList[0]
        else:
            l2Norm = float(l2Norm)
        self.assertEqual(l2Norm, 0)

#    @unittest.skipIf(cvar.scalarTypeStr in ("Integer", "MPInteger", "MPRational","MPRatComplex"), "ScalarType must not be Integer, MPInteger, MPRational or MPRatComplex")
    @unittest.skipIf(mypy.cvar.scalarTypeStr in ("Integer", "MPInteger" ), "ScalarType must not be Integer or MPInteger")
    def test_l2Norm_id(self):
        """
        l2Norm(id) = sqrt(8.0)
        """
        idMatID = mypy.get_identity_pID(self.row_level)
        l2NormID = mypy.normID(idMatID,2)
        l2Norm = mypy.traceID(l2NormID)
        if self.scalarTypeStr in ("Complex", "MPComplex"):
            # l2Norm = sum(map(float, l2Norm.split("+I*")))
            l2NormList = list(map(float, l2Norm.split("+I*")))
            # norm should be real
            if (len(l2NormList) == 2):
                self.assertEqual(l2NormList[1],0.0)
            l2Norm = l2NormList[0]
        elif self.scalarTypeStr in ("MPRational", "MPRatComplex"):
            if (self.scalarTypeStr == "MPRatComplex"):
                #    need to split complex into real and imaginary
                #    and confirm imaginary part is zero
                l2NormList = l2Norm.split("+I*")
                if (len(l2NormList) == 2):
                    l2ImagList = l2NormList[1].split("/") 
                    if (len(l2ImagList) == 2):
                        l2Imag = float(l2ImagList[0])/float(l2ImagList[1])
                    else:
                        l2Imag = float(l2ImagList[0])
                    self.assertEqual(l2Imag,0.0)
                l2Norm = l2NormList[0]
            #    need to split real rational into numerator and denominator
            l2NormList = l2Norm.split("/")
            l2Norm = float(l2NormList[0])/float(l2NormList[1])
        else:
            l2Norm = float(l2Norm)
        self.assertEqual(l2Norm, math.sqrt(8.0))

    def test_l2Norm_negated(self):
        """
        maxnorm(-A) = maxnorm(A)
        """
        # build -1
        neg1 = mypy.get_valID_from_valString("-1")
        # generate random matrix
        randMatAID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        negAMatID = mypy.scalar_mult(neg1, randMatAID)
        self.assertEqual(mypy.normID(randMatAID,2), mypy.normID(negAMatID,2))

    def test_l2Norm_random(self):
        """
        placing != entry 5 in zero matrix to make A
        maxnorm(A) = 5
        """
        # get a zero matrix and place a 5 some where in the matrix
        zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        # TO BE VERY CONFUSING we need an (i,j) index in the matrix we use:
        # for i: we know row_level <= 2^row_level so the index i=self.row_level is a valid row index
        # for j: we know col_level <= 2^col_level so the index j=self.col_level is a valid col index
        matAID = mypy.replace_scalar_in_matrix_by_string_and_coords(zeroMatID, self.row_level, self.col_level, "5");
        l2NormID = mypy.normID(matAID,2)
        l2Norm = mypy.traceID(l2NormID)
        if self.scalarTypeStr in ("Complex", "MPRatComplex", "MPComplex"):
#            l2Norm = sum(map(float, l2Norm.split("+I*")))
            l2NormList = list(map(float, l2Norm.split("+I*")))
            if (len(l2NormList)==2):
                self.assertEqual(l2NormList[1],0)
            l2Norm = l2NormList[0]
        else:
            l2Norm = float(l2Norm)
        self.assertEqual(l2Norm, 5)



class TestMatrixSparsity(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with mypy.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (mypy.cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_sparse_val(self):
        """
        Verify sparsity(A) = sparsity value for generating random matrix A. 
        """
        sparsity = 0
        n = (2**self.row_level)*(2**self.col_level)
        # get a nonzero random sparsity
        while sparsity == 0:
            sparsity = random.random()
        randMatID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], sparsity)
        # convert to number of zeros
        # requested sparsity may be inaccurate by less than one zero 
        self.assertEqual(int(round(sparsity * n)), int(mypy.matrix_count_entries(randMatID, "0")))

    def test_sparse_0(self):
        """
        Matrix generated with sparcity value 0 has sparsity 0.
        """
        onesMatID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, 1, 2, 0)
        self.assertEqual(int(mypy.matrix_count_entries(onesMatID, "0")), 0)

    def test_sparse_id(self):
        """
        Identity matrix of order n has sparsity (n-1)/n.
        """
        idMatID = mypy.get_identity_pID(self.row_level)
        self.assertEqual(int(mypy.matrix_count_entries(idMatID, "0")), (2**self.row_level-1)*(2**self.row_level))

    def test_sparse_0_mat(self):
        """
        Zero matrix has sparsity 0.
        """
        zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        self.assertEqual(int(mypy.matrix_count_entries(zeroMatID, "0")), (2**self.row_level) * (2**self.col_level))

class TestMatrixSum(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with mypy.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (mypy.cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3
        # get Zero matrix
        self.zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        # generate two random matrices of same size
        self.randMatAID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = mypy.get_valID_from_valString("-1")

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_add_additive_inverse(self):
        """
        A + (-A) = 0
        """
        a = mypy.scalar_mult(self.neg1, self.randMatAID)
        self.assertEqual(mypy.matrix_add(self.randMatAID, a), self.zeroMatID)

    def test_add_add_zero(self):
        """
        A + 0 = A
        """
        self.assertEqual(mypy.matrix_add(self.randMatAID, self.zeroMatID), self.randMatAID)

    def test_add_double(self):
        """
        A + A = 2A
        """
        scalMatID = mypy.get_valID_from_valString("2")
        prodMatID = mypy.scalar_mult(scalMatID, self.randMatAID)
        self.assertEqual(mypy.matrix_add(self.randMatAID, self.randMatAID), prodMatID)

    def test_add_commutativity(self):
        """
        A + B = B + A
        """
        a = mypy.matrix_add(self.randMatAID, self.randMatBID)
        b = mypy.matrix_add(self.randMatBID, self.randMatAID)
        self.assertEqual(a,b)

    def test_add_commutativity_op_shortcut(self):
        """
        If A+B in op store, calculating B+A shouldn't require additional calculations.
        """
        a = mypy.matrix_add(self.randMatAID, self.randMatBID)
        num_matrices_made = mypy.num_matrices_created()
        b = mypy.matrix_add(self.randMatBID, self.randMatAID)
        num_matrices_made2 = mypy.num_matrices_created()
        self.assertEqual(num_matrices_made, num_matrices_made2)


class TestMatrixDifference(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with mypy.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (mypy.cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3
        # get a Zero matrix
        self.zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        # generate two random matrices of same size
        self.randMatAID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = mypy.get_valID_from_valString("-1")

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_diff_additive_inverse(self):
        """
        A - A = 0
        """
        self.assertEqual(self.zeroMatID, mypy.matrix_diff(self.randMatAID, self.randMatAID))

    def test_diff_subtract_zero(self):
        """
        A - 0 = A
        """
        self.assertEqual(self.randMatAID, mypy.matrix_diff(self.randMatAID, self.zeroMatID))

    def test_diff_uni_invert(self):
        """
        0 - A = -A
        """
        self.assertEqual(mypy.matrix_diff(self.zeroMatID, self.randMatAID), mypy.scalar_mult(self.neg1, self.randMatAID))

    def test_diff_double(self):
        """
        A - (-A) = A + A
        """
        a = mypy.matrix_diff(self.randMatAID, mypy.scalar_mult(self.neg1, self.randMatAID))
        b = mypy.matrix_add(self.randMatAID, self.randMatAID)
        self.assertEqual(a,b)

    def test_diff_bin_invert(self):
        """
        A - B = -(B - A)
        """
        a = mypy.matrix_diff(self.randMatAID, self.randMatBID)
        b = mypy.matrix_diff(self.randMatBID, self.randMatAID)
        self.assertEqual(a, mypy.scalar_mult(self.neg1, b))

    @unittest.skip("not supported") # *.skip is python 2.7 and above
    def test_diff_size_mismatch(self):
        """
        Fail if size mismatch.
        """
        # generate matrix with one less row
        randMatCID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level-1, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # generate matrix with one less col
        randMatCID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level-1, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # could try using self.assertRaises() to show an exception is raised
        self.assertEqual(self.randMatAID, mypy.matrix_diff(self.randMatAID, randMatCID))


class TestMatrixValueCount(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with mypy.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        mypy.initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (mypy.cvar is in pylarc)
        self.scalarTypeStr = mypy.cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3

    def tearDown(self):
        # clean the matrix store after every test
        mypy.clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_count0_sparsity_rand(self):
        """
        Verify count0(A) = sparsity*size value for generating random matrix A. 
        """
        sparsity = 0
        n = (2**self.row_level)*(2**self.col_level)
        # get a nonzero random sparsity
        while sparsity == 0:
            sparsity = random.random()
        randMatID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, self.val_range[0], self.val_range[1], sparsity)
        # convert to number of zeros
        # requested sparsity may be inaccurate by less than one zero 
        self.assertEqual(int(round(sparsity * n)), int(mypy.matrix_count_entries(randMatID, "0")))

    def test_count0_sparse_0(self):
        """
        Matrix generated with sparcity value 0 has 0 zero entries.
        """
        onesMatID = mypy.matrix_random_matrixID(self.scalarTypeStr, self.row_level, self.col_level, 1, 2, 0)
        self.assertEqual(0, int(mypy.matrix_count_entries(onesMatID, "0")))

    def test_count0_id_mat(self):
        """
        nxn identiy matrix has n^2 - n zeros.
        """
        idMatID = mypy.get_identity_pID(self.row_level)
        self.assertEqual(int(mypy.matrix_count_entries(idMatID, "0")), ((2**self.row_level) - 1)*(2**self.row_level))

    def test_count1_id_mat(self):
        """
        nxn identiy matrix n ones.
        """
        idMatID = mypy.get_identity_pID(self.row_level)
        self.assertEqual(int(mypy.matrix_count_entries(idMatID, "1")), 2**self.row_level)

    def test_count0_0_mat(self):
        """
        nxn zero matrix has n^2 zeros.
        """
        n = (2**self.row_level)*(2**self.col_level);
        zeroMatID = mypy.get_zero_pID(self.row_level, self.col_level)
        self.assertEqual(n, int(mypy.matrix_count_entries(zeroMatID, "0")))


#needed for python 2.6 (without unittest2)
if __name__ == "__main__":
    unittest.main()


