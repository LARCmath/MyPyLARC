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
import math
from ctypes import *

import potential
import centdiff

## \file hamilton_morse.py
#
#  \brief This routine calculates the Hamiltonian matrix H for a Morse
#   potential, and a suitable H' for the power method.
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
    hardcoded to the Morse ocillator), a given minimum and maximum
    point on the real axis, and a discretization distance determined by the
    number of points in the discretization and the interval [xmin,xmax]. 

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
        rmax = float(sys.argv[2])
    else:
        level = 3
        rmax = 10.0
    rmin = 0.0

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = max(10,level)
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    mypy.create_report_thread(1800)
    verbose = 1
    mypy.initialize_larc(mat_store_exp, op_store_exp,
                  max_level,regionbitparam,zeroregionbitparam,verbose)

    scalarTypeStr = mypy.cvar.scalarTypeStr

    # We are going to discretize the real line between the two points rmin
    # and rmax. There will be a total of 2**level points in the discretization,
    # evenly spaced; the distance between adjacent points is given by delta.
    npts = 2**level
    delta = (rmax-rmin)/(npts-1)

    # create a list of r values to pass into the potential function
    r = [rmin + i*delta for i in range(npts)]

    # the morse oscillator hamiltonian is
    # given by -\hbar^2/2m d^2/dr^2 + V_{morse}(r;D_e,r_e,\alpha)
    # where: D_e is the well depth of the potential
    #           (thus -D_e is the minimum energy of the potential)
    #        r_e is the distance at which this minimum is achieved
    #           (the equilibrium bond distance)
    #        \alpha controls the width of the well

    # create and store the potential matrix
    # one set of suitable parameters
    #            r_e = 1.28 Å, D = 4.17 eV and α = 1.85 Å^{-1}
    # (Anil K. Shukla, Jean H. Futrell, "Ion Collision Theory", in
    # Encyclopedia of Spectroscopy and Spectrometry, 1999, Academic Press,
    # pp. 954--963:
    # "the Morse potential with these parameters matches the experimentally
    # determined scattering of protons by argon to better than 1% accuracy over
    # the interval from 0 to 5 Å.")

    # note that the reference listed the α units as Å^3, but this is 
    # obviously incorrect since it appears in the exponential multiplied by
    # a distance, not an inverse distance cubed.

    # now for the fun part: getting everything into consistent units
    # we will use Angstroms (1Å = 10^{-10} meters) as our distance unit,
    # and electron-volts (1eV = 1.60218e-19 Joules) as our energy unit.

    r_min_energy = 1.28
    well_depth = 4.17
    alpha = 1.85
    vr = potential.morse_potential(r,alpha,well_depth,r_min_energy)
    vr_str = mypy.map_to_str(vr, scalarTypeStr)
    potVecID = mypy.row_major_list_to_store(vr_str, level, 0, 1)
    # potVecID is a column vector, but we need to put the potential energy
    # values on the diagonal of a 2^level x 2^level matrix
    potMatID = potential.put_col_vector_on_diagonal(potVecID,level)
    #print("potMatID:")
    #mypy.print_naive(potMatID)

    # save the potential to disk for future plotting
    filename = "Data/morsePot_L" + str(level) + "_rmax" + str(rmax)
    file1 = open(filename,'w')
    for vpt in vr:
        file1.write(str(vpt)+"\n")
    file1.close()


    # create and store the kinetic (second derivative) matrix
    # from Wikipedia
    # speed of light c = 2.99792458e18 Å/s 
    # hbar = 6.582119569e−16 eV⋅s
    # mass of a proton: 1.007276466621(53) amu = 9.3827208816(29)e8 eV/c^2
    # mass of Argon: 39.9623831238(24) amu
    # reduced mass of two particles is given by (1/m1 + 1/m2)^{-1}:
    # (1/1.007276466621 + 1/39.9623831238)^{-1)
    #      = ( 0.992776098 + 0.025023533 )^{-1}
    #      = ( 1.017799631 )^{-1} =  0.982511656 amu

    # => conversion factor, calculated from proton mass numbers:
    #    931494102.416110814 (eV/c^2)/amu
    # using conversion factor, reduced mass is 9.15203813e8 eV/c^2

    # units of -hbar^2/2m: (eV*s)^2/(eV/c^2) == eV*s^2/c^2
    # -hbar^2/2m = -(6.582119569e-16)^2/2.0/9.15203813e8 
    #            = -2.366920756e-40 eV⋅s^2/c^2
    #  * c^2     => -2.1272822872e-3 eV⋅Å^2

    ddr2ID = centdiff.build_central_diff_matrix(level,delta)
    scale_factor = -0.0021272822872
    scaleStr = mypy.value_to_string(scale_factor,scalarTypeStr)
    scaleID = mypy.get_valID_from_valString(scaleStr)
    kinMatID = mypy.scalar_mult(scaleID,ddr2ID)
    #print("kinMatID:")
    #mypy.print_naive(kinMatID)

    # add to get Hamiltonian matrix
    hamMatID = mypy.matrix_add(kinMatID,potMatID)
    #print("hamMatID:")
    #mypy.print_naive(hamMatID)
    mypy.fprint_larcMatrixFile(hamMatID,
         "Data/basicMorse_L" + str(level) + "_rmax" + str(rmax) + ".json")

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
    #the unit vector with a "1" in a position near the minimum energy
    minpos = int((r_min_energy-rmin)/delta)
    mypy.info_set(mypy.COMMENT,invHamID, str(minpos) +
            " is a good position for the nonzero in the starting vector")

    # write the matrix to disk
    mypy.fprint_larcMatrixFile(invHamID,
         "Data/invertMorse_L" + str(level) + "_rmax" + str(rmax) + ".json")

