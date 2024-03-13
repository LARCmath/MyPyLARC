#!/usr/bin/env python3

 #*##############################################################*#
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
 #*##############################################################*#

from __future__ import print_function

import os
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mypy
import math
from ctypes import *

import potential
import centdiff

## \file hamilton_harm.py
#  \brief This routine calculates the Hamiltonian matrix H for a harmonic oscillator, and a suitable H' for the power method
#
#   In order to use the power method to find the lowest-energy eigenvector for
#   this potential, the program first calculates H, then modifies it by
#   subtracting H from an identity matrix times a suitable constant e: 
#   H' = eI - H. The power method finds the eigenvalue of H' with largest
#   magnitude (E'_max) and its associated eigenvector. The smallest eigenvalue
#   of H, E_min = e-E'_max, is associated to this same eigenvector.

if __name__ == '__main__':
    '''
    This routine generates the input matrix for the power method. It first
    calculates the Hamiltonian matrix H for a given potential (currently
    hardcoded to the harmonica ocillator), a given minimum and maximum
    point on the real axis (because we are hardcoding the harmonic oscillator,
    we make the interval symmetric about zero), and a discretization distance
    determined by the number of points in the discretization and the interval 
    [xmin,xmax]. 

    The power method can be used to find the maximal eigvenvalue of a matrix,
    but in this case we want the minimal eigenvalue. Thus we modify our
    Hamiltonian matrix to invert the magnitude of the eigenvalues without
    changing the eigenvectors. To do this, we assume the maximum of the
    potential function on the interval is at xmin, and that therefore the
    largest element in the Hamiltonian matrix is at the (0,0) location. We
    use the value A of the (0,0) element as a crude estimate of the maximal
    eigenvalue of the matrix. We calculate and save H'=A*I-H; the maximal
    eigenvalue E'_max of H' is related to the minimal eigenvalue E_0 of H by
    E_0 = A-E'_max, and the associated eigenvector is the same in both cases.
    
    The InfoStore is used to record the value A as a string in the header of
    the larcMatrix json file, where it can be read by the power method
    routine. Once it finds E'_max, it can use A to calculate E_0.
    '''
    if 3 == len(sys.argv):
        level = int(sys.argv[1])
        xmax = float(sys.argv[2])
    else:
        level = 3
        xmax = 10.0

    # since the harmonic oscillator is symmetric about x=0, we just set
    # xmin to be the negatie of xmax
    xmin = -xmax

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

    # We are going to discretize the real line between the two points xmin
    # and xmax. There will be a total of 2**level points in the discretization,
    # evenly spaced; the distance between adjacent points is given by delta.
    npts = 2**level
    delta = (xmax-xmin)/(npts-1)

    # create a list of x values to pass into the potential function
    x = [xmin + i*delta for i in range(npts)]

    # the harmonic oscillator hamiltonian, in scaled units, is
    # given by -1/2 d^2/dx^2 + 1/2 x^2

    # create and store the potential matrix (1/2 x^2)
    min_energy = 0.0
    scale_factor = 0.5
    vx = potential.harmonic_potential(x,scale_factor,min_energy)
    vx_str = mypy.map_to_str(vx, scalarTypeStr)
    potVecID = mypy.row_major_list_to_store(vx_str, level, 0, 1)
    # potVecID is a column vector, but we need to put the potential energy
    # values on the diagonal of a 2^level x 2^level matrix
    potMatID = potential.put_col_vector_on_diagonal(potVecID,level)
    #print("potMatID:")
    #mypy.print_naive(potMatID)

    # save the potential to disk for future plotting
    filename = "Data/harmonicPot_L" + str(level) + "_xmax" + str(xmax)
    file1 = open(filename,'w')
    for vpt in vx:
        file1.write(str(vpt)+"\n")
    file1.close()

    # create and store the kinetic (second derivative) matrix
    ddx2ID = centdiff.build_central_diff_matrix(level,delta)
    scale_factor = -0.5
    scaleStr = mypy.value_to_string(scale_factor,scalarTypeStr)
    scaleID = mypy.get_valID_from_valString(scaleStr)
    kinMatID = mypy.scalar_mult(scaleID,ddx2ID)
    #print("kinMatID:")
    #mypy.print_naive(kinMatID)

    # add to get Hamiltonian matrix
    hamMatID = mypy.matrix_add(kinMatID,potMatID)
    #print("hamMatID:")
    #mypy.print_naive(hamMatID)
    mypy.fprint_larcMatrixFile(hamMatID,
         "Data/basicHarmonic_L" + str(level) + "_xmax" + str(xmax) + ".json")

    # modify as describe above so that smallest eigenvalue becomes largest
    extremeID = mypy.get_scalarID_from_pID_and_coords(hamMatID,0,0)
    print("extreme value is ")
    mypy.print_naive(extremeID)
    bigDiagID = mypy.scalar_mult(extremeID, mypy.get_identity_pID(level))
    invHamID = mypy.matrix_diff(bigDiagID,hamMatID)
    #print("invHamID:")
    #mypy.print_naive(invHamID)

    #write the value of extremeID into the infoStore so that the 
    #power method will have access to it
    mypy.info_set(mypy.OTHERINFO,invHamID,
	"subtract the power-method-calculated eigenvalue from this to" +
        " get the correct eigenvalue for the ground state of the Hamiltonian")
    mypy.info_set(mypy.OTHERMATRIX,invHamID,
         mypy.get_scalar_value_string(extremeID))

    #help the power method by giving it a good starting point: tell it to use
    #the unit vector with a "1" in a position near the minimum energy (at x=0)
    minpos = int((0-xmin)/delta)
    mypy.info_set(mypy.COMMENT,invHamID, str(minpos) +
            " is a good position for the nonzero in the starting vector")

    mypy.fprint_larcMatrixFile(invHamID,
         "Data/invertHarmonic_L" + str(level) + "_xmax" + str(xmax) + ".json")
