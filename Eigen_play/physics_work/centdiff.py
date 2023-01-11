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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mypy
from ctypes import *
import json    # for loading parameter files

## \file centdiff.py
#  \brief Routines for calculating second derivatives via matrix multiplication and the central difference formula

## \brief builds a tridiagonal matrix with 1 -2 1 on diagonals
def recursive_build_central_diff_matrix(level):
    """
    This routine returns a tuple consisting of the diagonal and
    off-diagonal blocks of a tridiagonal matrix with -2 on the
    diagonal and 1 on the off-diagonals. This is a matrix instantiation
    of the central difference formula for the second derivative, so
    multiplication of a vector by this matrix produces a result that
    is approximately the second derivative of that vactor. We have left
    out the 1/(\Delta x)^2 scaling factor, where \Delta x is the (fixed)
    distance between the points in the discretization.
    """

    if (level==1):
	# diagonal is [-2,1;1,-2]; right off-diagonal is [0,0;1,0]; 
	# left off-diagonal is [0,1;0,0]; 
        minus2ID = mypy.get_valID_from_valString("-2")
        oneID = mypy.get_identity_pID(0);
        zeroID = mypy.get_zero_pID(0,0)
        rodID = mypy.get_pID_from_four_sub_pIDs(zeroID,
                    zeroID,oneID,zeroID,1,1)
        lodID = mypy.get_pID_from_four_sub_pIDs(zeroID,
                    oneID,zeroID,zeroID,1,1)
        diagID = mypy.get_pID_from_four_sub_pIDs(minus2ID,
                    oneID,oneID,minus2ID,1,1)
        return diagID, rodID, lodID

    # larger matrices have diagonal blocks equal to smdiagID
    # and off-diagonal blocks built from the smaller ones and three
    # zero matrices
    smdiagID, srodID, slodID = recursive_build_central_diff_matrix(level-1)
    diagID = mypy.get_pID_from_four_sub_pIDs(
                   smdiagID,srodID,slodID,smdiagID,level,level)
    zeroID = mypy.get_zero_pID(level-1,level-1)
    rodID = mypy.get_pID_from_four_sub_pIDs(
                   zeroID,zeroID,srodID,zeroID,level,level)
    lodID = mypy.get_pID_from_four_sub_pIDs(
                   zeroID,slodID,zeroID,zeroID,level,level)
    return diagID, rodID, lodID

## \brief builds a matrix which approximates the second derivative
def build_central_diff_matrix(level, deltaX):
    """
    This routine returns a specific tridiagonal matrix which is the matrix
    instantiation of the central difference formula for the second
    derivative, given that deltaX is the (fixed) spacing between points
    in the discretization. Multiplication of a vector by this matrix produces
    a result that is approximately the second derivative of that vactor.
    """
    scalarTypeStr = mypy.cvar.scalarTypeStr
    # deal with special case of 1x1 "tridiagonal" matrix
    if (level==0): 
        element = -2.0/deltaX/deltaX
        valStr = mypy.value_to_string(element,scalarTypeStr)
        return mypy.get_valID_from_valString(valStr)

    scaleFactor = 1.0/deltaX/deltaX
    scaleStr = mypy.value_to_string(scaleFactor,scalarTypeStr)
    scaleID = mypy.get_valID_from_valString(scaleStr)

    # deal with special case of 2x2 "tridiagonal" matrix
    if (level==1):
        diagID, rodID, lodID = recursive_build_central_diff_matrix(1)
        return mypy.scalar_mult(scaleID,diagID)

    # use recursion to build up the central difference matrix of
    # dimension 2^level x 2^level
    diagID, rodID, lodID = recursive_build_central_diff_matrix(level-1)
    unscaledDiffID = mypy.get_pID_from_four_sub_pIDs(
                diagID,rodID,lodID,diagID,level,level)

    # include the inverse delta^2 term in the returned matrix
    return mypy.scalar_mult(scaleID,unscaledDiffID)


if __name__ == '__main__':

    #*##################################################################*#
    # Set the level (matrices are 2**level by 2**level                   #
    # and the verbosity (0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG, 4=ALL)    #
    #*##################################################################*#

    if 4 == len(sys.argv):
        level = int(sys.argv[1]) 
        verbose = int(sys.argv[2])
        LARC_verbose = int(sys.argv[3])
    else:
        print("\nThis program requires three commandline integer inputs:")
        print("   level:  matrices will be 2**level by 2**level")
        print("   verbose:  the verbosity level for this program")
        print("   LARC_verbose:  the verbosity level for the LARC package")
        sys.exit()

    scalarTypeStr = mypy.cvar.scalarTypeStr
    if scalarTypeStr  in ('Integer','MPInteger'):
        print("\nThis routine does not work with integer types!")
        sys.exit(1)

    # read the parameter file into a python dictionary
    with open('../../InitParams/power_method.init_params','r') as init_file:
        init_param = json.load(init_file)
        for p in init_param['small']:
            if (verbose > 1):
                print('MatrixExponent: %d' %(p['matrix_exponent']))
                print('OpExponent: %d' %(p['op_exponent']))
                print('MaxLevel: %d' %(p['max_level']))
                print('RegionBitParam: %d' %(p['regionbitparam']))
                print('ZeroRegionBitParam: %d' %(p['zeroregionbitparam']))
                print('ReportIntervalSecs: %d' %(p['report_interval_seconds']))
                print('MinMemRequiredGiB: %d' %(p['min_memGiB_required']))
                print('Verbose: %d' %(p['verbose']))
                print('')
            matrix_exponent = p['matrix_exponent']
            op_exponent = p['op_exponent']
            max_level= p['max_level']
            regionbitparam = p['regionbitparam']
            zeroregionbitparam = p['zeroregionbitparam']
            report_interval_seconds = p['report_interval_seconds']
            min_memGiB_required = p['min_memGiB_required']
            p_verbose = p['verbose']

    # warn if the commandline value for LARC_verbose differs from the parameter file value for p_verbose        
    if (LARC_verbose > 0):
        if (LARC_verbose != p_verbose):
            print("NOTE: This program uses commandline (LARC_verbose = %d) " %LARC_verbose)
            print("      rather than the parameter file (p_verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    # initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,LARC_verbose)

    # may as well set the discretization distance to 1
    # (though varying this is recommended when testing)
    delta = 1.0

    # generate test matrix for viewing on screen (if small) or 
    # in file (if large)
    tridID = build_central_diff_matrix(level,delta)

    if (level<6):
        mypy.print_naive(tridID)
    else:
        mypy.fprint_naive(tridID,"./tridID.out")
