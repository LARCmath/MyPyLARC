#!/usr/bin/env python3

 #*##############################################################*#
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
 #*##############################################################*#

from __future__ import print_function

import os
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mypy
import random
import math
from ctypes import *

## \file potential.py
#  \brief Contains utility routines for calculating potential energy surfaces.

## \brief makes a diagonal square matrix with the values of a row vector on its diagonal
def put_row_vector_on_diagonal(vectorID,level):
    '''
    This routine takes as input the matrixID for a LARC row vector,
    and returns the matrixID of a diagonal square matrix with the row
    vector elements on the diagonal.
    '''
    # row vector is [aID,bID,-1, -1], so we want submatrices 0 and 1
    leftID = mypy.get_pID_of_indexed_submatrix(vectorID,0)
    rightID = mypy.get_pID_of_indexed_submatrix(vectorID,1)
    #print("leftID=%d, rightID=%d" %(leftID,rightID))
    if (level==1):
        # two scalars are just put on the diagonal of a 2x2 matrix
        zeroID = mypy.get_zero_pID(0,0)
        #mypy.print_naive(leftID)
        #mypy.print_naive(rightID)
        #print("constructing 2x2 matrix")
        lev1ID = mypy.get_pID_from_four_sub_pIDs(
               leftID,zeroID,zeroID,rightID,1,1)
        #print("matrix is:")
        #mypy.print_naive(lev1ID)
        return lev1ID

    # recursive calls returns diagonal matrices of one lower level
    leftdiagID = put_row_vector_on_diagonal(leftID,level-1)
    # print("leftdiagID is", leftdiagID)
    rightdiagID = put_row_vector_on_diagonal(rightID,level-1)
    # print("rightdiagID is", rightdiagID)

    # make the two matrices diagonal blocks of a larger matrix
    zeroID = mypy.get_zero_pID(level-1,level-1)
    return mypy.get_pID_from_four_sub_pIDs(
           leftdiagID,zeroID,zeroID,rightdiagID,level,level)

## \brief makes a diagonal square matrix with the values of a column vector on its diagonal
def put_col_vector_on_diagonal(vectorID,level):
    '''
    This routine takes as input the matrixID for a LARC column vector,
    and returns the matrixID of a diagonal square matrix with the column
    vector elements on the diagonal.
    '''
    # column vector is [aID,-1,bID,-1], so we want submatrices 0 and 2
    topID = mypy.get_pID_of_indexed_submatrix(vectorID,0)
    bottomID = mypy.get_pID_of_indexed_submatrix(vectorID,2)
    #print("topID=%d, bottomID=%d" %(topID,bottomID))
    if (level==1):
        # two scalars are just put on the diagonal of a 2x2 matrix
        zeroID = mypy.get_zero_pID(0,0)
        #mypy.print_naive(topID)
        #mypy.print_naive(bottomID)
        #print("constructing 2x2 matrix")
        lev1ID = mypy.get_pID_from_four_sub_pIDs(
               topID,zeroID,zeroID,bottomID,1,1)
        #print("matrix is:")
        #mypy.print_naive(lev1ID)
        return lev1ID

    # recursive calls returns diagonal matrices of one lower level
    topdiagID = put_col_vector_on_diagonal(topID,level-1)
    #print("topdiagID is", topdiagID)
    #mypy.print_naive(topdiagID)
    bottomdiagID = put_col_vector_on_diagonal(bottomID,level-1)
    #print("bottomdiagID is", bottomdiagID)
    #mypy.print_naive(bottomdiagID)

    # make the two matrices diagonal blocks of a larger matrix
    zeroID = mypy.get_zero_pID(level-1,level-1)
    return mypy.get_pID_from_four_sub_pIDs(
           topdiagID,zeroID,zeroID,bottomdiagID,level,level)

## \brief Produces a 1D potential energy surface for a harmonic oscillator
def harmonic_potential(x, ctimes, e_min):
    # x is list, return list 
    # the harmonic oscillator potential is just a parabola
    # the parameter ctimes affects the slope of the parabola
    #    and the spacing of the energy eigenvalues
    # the potential has minimum energy e_min
    # e_min = -100000.0
    # ctimes = 1
    y = [(ctimes*i*i + e_min) for i in x]
    return y

## \brief Produces a 1D potential energy surface for the Morse potential
def morse_potential(r, alpha, e_min, r_min):
    # r is list, return list
    # alpha is a scaling factor that determines how fast the
    #     exponentials go to zero/infinity
    # the minimum energy e_min is achieved when r=r_min
    # energy approaches zero asymptotically as r increases
    pow1 = [-alpha*(t-r_min) for t in r]
    y = [ e_min*(math.exp(2*t)-2*math.exp(t)) for t in pow1]
    return y

if __name__ == '__main__':

    if 2 == len(sys.argv):
        level = int(sys.argv[1])
    else:
        level = 2

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp, op_store_exp,
                  max_level,regionbitparam,zeroregionbitparam,verbose)

    scalarTypeStr = mypy.cvar.scalarTypeStr
    dim = 2**level

    # build random array
    if scalarTypeStr in ("Real", "MPRational", "MPReal"):
        randVals = [ random.random() for i in range(dim)]
    elif scalarTypeStr in ("Complex", "MPRatComplex", "MPComplex"):
        randVals = [ np.complex(random.random(),random.random()) for i in range(dim)]
    elif scalarTypeStr in ("Integer", "MPInteger"):
        randVals = [ random.randrange(0,10001) for i in range(dim)]
    else:
        raise Exception('Do not know how to build matrix for type %s.'%scalarTypeStr)
    randVals_strs = mypy.map_to_str(randVals, scalarTypeStr)
#    print(randVals)
#    print(randVals_strs)

#    potVecID = mypy.row_major_list_to_store(randVals_strs, 
#                 level, 0, 1)
#    mypy.print_naive(potVecID)
#    print("")
#    potMatID = put_col_vector_on_diagonal(potVecID,level)
#    mypy.print_naive(potMatID)
    potVecID = mypy.row_major_list_to_store(randVals_strs,
                 0, level, dim)
    mypy.print_naive(potVecID)
    print("")
    potMatID = put_row_vector_on_diagonal(potVecID,level)
    mypy.print_naive(potMatID)

