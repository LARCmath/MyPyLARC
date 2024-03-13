#!/usr/bin/env python3

# NOTE: to run on command line
#   python3 -m unittest -v example_unittest_scalars.py
# To run individual test classes from the module, do (for example) 
#   python3 -m unittest -v example_unittest_scalars.TestMatrixMaxnorm

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
import MyPyLARC as myp
# from pylarc import *
from ctypes import *
import math
import unittest


##
# \file example_unittest_scalars.py
#
# \brief This contains a sample set of unittests for scalars in LARC.
#
#  To use the python unittest module and run this code type:
#  python3 -m unittest example_unittest_scalars.py
#

class TestMatrixMaxnorm(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarTypeStr = myp.cvar.scalarTypeStr

    def test_twobytwo_doubling(self):
        level = 1
        dim_whole = 2**level

        # with myp.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        myp.initialize_larc(26,24,10,-1,-1,0)   # run quiet

        if self.scalarTypeStr in ("Integer", "MPInteger"):
            arr_a = [4, 0, 0, 3]
            arr_b = [8, 0, 0, 6]
        else:
            arr_a = [.75, 0, 0, .125]
            arr_b = [1.5, 0, 0, .25]
        a_ID = myp.row_major_list_to_store(myp.map_to_str(arr_a, self.scalarTypeStr), level, level, dim_whole)
        b_ID = myp.row_major_list_to_store(myp.map_to_str(arr_b, self.scalarTypeStr), level, level, dim_whole)
        self.assertEqual(myp.matrix_add(a_ID, a_ID), b_ID)

    def test_twobytwo_subtraction(self):
        level = 1
        dim_whole = 2**level

        # with myp.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        myp.initialize_larc(26,24,10,-1,-1,0)    # run quiet

        if self.scalarTypeStr in ("Integer", "MPInteger"):
            arr_a = [4, 0, 0, 3]
            arr_b = [8, 0, 0, 6]
        else:
            arr_a = [.75, 0, 0, .125]
            arr_b = [1.5, 0, 0, .25]
        a_ID = myp.row_major_list_to_store(myp.map_to_str(arr_a, self.scalarTypeStr), level, level, dim_whole)
        b_ID = myp.row_major_list_to_store(myp.map_to_str(arr_b, self.scalarTypeStr), level, level, dim_whole)
        self.assertEqual(myp.matrix_diff(b_ID, a_ID), a_ID)

    def test_twobytwo_negation(self):
        level = 1
        dim_whole = 2**level

        # with myp.stdout_redirected():
        #     initialize_larc(26,24,10,-1,-1,1)
        myp.initialize_larc(26,24,10,-1,-1,0)    # run quiet

        if self.scalarTypeStr in ("Integer", "MPInteger"):
            arr_a = [4, 0, 0, 3]
            arr_b = [8, 0, 0, 6]
        else:
            arr_a = [.75, 0, 0, .125]
            arr_b = [1.5, 0, 0, .25]
        a_ID = myp.row_major_list_to_store(myp.map_to_str(arr_a, self.scalarTypeStr), level, level, dim_whole)
        zero_ID = myp.get_zero_pID(level, level)
        self.assertEqual(myp.matrix_diff(a_ID, a_ID), zero_ID)

    def test_subtract_self(self):
        # with myp.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        myp.initialize_larc(26,24,10,-1,-1,0)   # run quiet

        if self.scalarTypeStr in ("Integer", "MPInteger"):
            pID = myp.get_valID_from_valString("40")
        else:
            pID = myp.get_valID_from_valString("0.4")
        scalarM1 = myp.get_valID_from_valString("-1")
        neg_pID = myp.matrix_mult(scalarM1, pID)
        zero = myp.get_zero_pID(0, 0)
        self.assertEqual(myp.matrix_add(neg_pID, pID), zero)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex", "MPReal", "MPComplex", "MPRational", "MPRatComplex"), "ScalarType must be real or complex ")
    def test_zeroBitThresh_read_below(self):
        """
        when zeroBitThresh = z any scalar y with abs(y) < 2^(-z-1) will be
        in the zero region in particular, y = 2^(-z-2)
        """
        zeroBitThresh = 3

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, -1, zeroBitThresh,1)
        myp.initialize_larc(26, 24, 10, -1, zeroBitThresh,0)  # run quiet

        small = 1.0/(2**(zeroBitThresh+2))
        smallID = myp.get_valID_from_valString(myp.value_to_string(small, self.scalarTypeStr))
        zero = myp.get_zero_pID(0, 0)
        self.assertEqual(smallID, zero)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_zeroBitThresh_read_equal(self):
        """
        when zeroBitThresh = t anything y with y>=2^t will be nonzero.
        """
        zeroBitThresh = 3

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, -1, zeroBitThresh,1)
        myp.initialize_larc(26, 24, 10, -1, zeroBitThresh,0)  # run quiet

        notSmall = 1.0/(2**zeroBitThresh)
        notSmallID = myp.get_valID_from_valString(myp.value_to_string(notSmall, self.scalarTypeStr))
        zero = myp.get_zero_pID(0, 0)
        self.assertNotEqual(notSmallID, zero)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_zeroBitThresh_read_above(self):
        """
        when zeroBitThresh = t anything y with y>=2^t will be nonzero.
        """
        zeroBitThresh = 3

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, -1, zeroBitThresh,1)
        myp.initialize_larc(26, 24, 10, -1, zeroBitThresh,0)  # run quiet

        notSmall = 3.0/(2**(zeroBitThresh+1))
        notSmallID = myp.get_valID_from_valString(myp.value_to_string(notSmall, self.scalarTypeStr))
        zero = myp.get_zero_pID(0, 0)
        self.assertNotEqual(notSmallID, zero)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_sighash_read_roundDown(self):
        """
        When sighash = t, real number x 
            = (-1)^sign (1 + \sum_1^{52} b_{52-i}2^{-i}) * 2^{e-1023}
        is rounded down to 
            = (-1)^sign (1 + \sum_1^{t-1} b_{52-i}2^{-i}) * 2^{e-1023}
        or rounded up to 
            = (-1)^sign (1 + \sum_1^{t-1} (b_{52-i}+1)2^{-i}) * 2^{e-1023}
        if b_{52-t} == 1.

        In particular, the largest power of two we can add to 1 and not change it
        is 2^-(t+1). Here we test that we can add it. 
        """
        sighash = 4

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, sighash,-1,1)
        myp.initialize_larc(26, 24, 10, sighash,-1,0)  # run quiet

        oneID = myp.get_valID_from_valString("1")
        m = 1 + 1.0/(2**(sighash+2))
        pID = myp.get_valID_from_valString(myp.value_to_string(m, self.scalarTypeStr))

        self.assertEqual(pID, oneID)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_sighash_read_roundUp(self):
        """
        When sighash = t, real number x 
            = (-1)^sign (1 + \sum_1^{52} b_{52-i}2^{-i}) * 2^{e-1023}
        is rounded down to 
            = (-1)^sign (1 + \sum_1^{t-1} b_{52-i}2^{-i}) * 2^{e-1023}
        or rounded up to 
            = (-1)^sign (1 + \sum_1^{t-1} (b_{52-i}+1)2^{-i}) * 2^{e-1023}
        if b_{52-t} == 1.

        In particular, the largest power of two we can add to 1 and not change it
        is 2^-(t+1). Here we test that a 2^-t doesn't work (because it rounds
        up). 
        """
        sighash = 4

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, sighash,-1,1)
        myp.initialize_larc(26, 24, 10, sighash,-1,0)  # run quiet

        oneID = myp.get_valID_from_valString("1")
        m = 1 + 1.0/(2**sighash)
        pID = myp.get_valID_from_valString(myp.value_to_string(m, self.scalarTypeStr))

        self.assertNotEqual(pID, oneID)

    @unittest.skipUnless(myp.cvar.scalarTypeDef in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_sighash_read_noRound(self):
        """
        when sighash = t, real number x 
            = (-1)^sign (1 + \sum_1^{52} b_{52-i}2^{-i}) * 2^{e-1023}
        is rounded down to 
            = (-1)^sign (1 + \sum_1^{t-1} b_{52-i}2^{-i}) * 2^{e-1023}
        or rounded up to 
            = (-1)^sign (1 + \sum_1^{t-1} (b_{52-i}+1)2^{-i}) * 2^{e-1023}
        if b_{52-t} == 1.

        In particular, the largest power of two we can subtract from 1 and 
        not change it is 2^-(t+1). Here we show that it works. 
        """
        sighash = 4

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, sighash,-1,1)
        myp.initialize_larc(26, 24, 10, sighash,-1,0)  # run quiet

        oneID = myp.get_valID_from_valString("1")
        m = 1 - 1.0/(2**(sighash+1))
        pID = myp.get_valID_from_valString(myp.value_to_string(m, self.scalarTypeStr))

        self.assertEqual(pID, oneID)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_sighash_read_round(self):
        """
        when sighash = t, real number x 
            = (-1)^sign (1 + \sum_1^{52} b_{52-i}2^{-i}) * 2^{e-1023}
        is rounded down to 
            = (-1)^sign (1 + \sum_1^{t-1} b_{52-i}2^{-i}) * 2^{e-1023}
        or rounded up to 
            = (-1)^sign (1 + \sum_1^{t-1} (b_{52-i}+1)2^{-i}) * 2^{e-1023}
        if b_{52-t} == 1.

        In particular, the largest power of two we can subtract from 1 and 
        not change it is 2^-(t+1). Here we show subtracting 2^-t does not
        work. (This one is weird: remember to shift by power of two so that 
        the necessary 'implicit' bit is there before rounding). 
        """
        sighash = 4

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, sighash,-1,1)
        myp.initialize_larc(26, 24, 10, sighash,-1,0)  # run quiet

        oneID = myp.get_valID_from_valString("1")
        m = 1 - 1.0/(2**(sighash))
        pID = myp.get_valID_from_valString(myp.value_to_string(m, self.scalarTypeStr))

        self.assertNotEqual(pID, oneID)

    # @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    @unittest.skipIf(True, "Not sure why sqrt2*sqrt2 = 2.0 in floating point")
    def test_sqrt2_squared_is_2(self):
        """
        sqrt(2) squared is the same as 2.
        """

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, -1,-1,1)
        myp.initialize_larc(26, 24, 10, -1,-1,0)  # run quiet

        sqrt2ID = myp.get_valID_from_valString(myp.value_to_string(math.sqrt(2), self.scalarTypeStr))
        twoID = myp.get_valID_from_valString("2")

        self.assertEqual(myp.matrix_mult(sqrt2ID, sqrt2ID), twoID)

    @unittest.skipIf(myp.cvar.scalarTypeStr in ("Integer", "MPInteger"), "ScalarType must not be Integer or MPInteger")
    def test_sqrt2_formula_works(self):
        """
        (sqrt2+sqrt2)/2 is the same as sqrt2.
        """

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, -1,-1,1)
        myp.initialize_larc(26, 24, 10, -1,-1,0)  # run quiet

        halfID = myp.get_valID_from_valString(myp.value_to_string(0.5, self.scalarTypeStr))
        sqrt2ID = myp.get_valID_from_valString(myp.value_to_string(math.sqrt(2), self.scalarTypeStr))
        doubleID = myp.matrix_add(sqrt2ID, sqrt2ID)
        resultID = myp.matrix_mult(halfID, doubleID)
        failure_msg = "matrices displayed above"
        if resultID != sqrt2ID:
            print("Failure in test_sqrt2_formula_works - printing relevant matrices.")
            print("Half matrix:")
            myp.print_naive(halfID)
            print("Sqrt2 matrix:")
            myp.print_naive(sqrt2ID)
            print("Double Sqrt2 matrix:")
            myp.print_naive(doubleID)
            print("Result matrix:")
            myp.print_naive(resultID)
        self.assertEqual(resultID, sqrt2ID, failure_msg)

    @unittest.skipUnless(myp.cvar.scalarTypeStr in ("Real", "Complex"), "ScalarType must be Real or Complex")
    def test_inv_sqrt2_formula_works(self):
        """
        ((1/sqrt2)^2) is the same as 1/2.
        """

        # with myp.stdout_redirected():
        #    initialize_larc(26, 24, 10, -1,-1,1)
        myp.initialize_larc(26, 24, 10, -1,-1,0)  # run quiet

        halfID = myp.get_valID_from_valString(myp.value_to_string(0.5, self.scalarTypeStr))
        resultID = myp.matrix_mult(myp.cvar.packedID_inv_sqrt_2, myp.cvar.packedID_inv_sqrt_2)
        self.assertEqual(resultID, halfID)


