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
 # Additional contributors are listed in "LARCContributors".      #
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


from __future__ import print_function, division

import os
import glob
import sys
import numpy as np
import datetime
import json

## \file numpy_run.py
#
#  \brief Use numpy on Sycamore quantum supremacy computation.
#
if __name__ == '__main__':

    verbose = 0
    debug = 0

    #*################################*#
    #*    Basic Parameter Setting     *#
    #*################################*#

    
    print("\n    ######################################################################")
    print("    ##  Create matrices for Sycamore quantum supremacy computation.     ##")
    print("    ######################################################################\n")
    

    #*#########################################*#
    #*    Output paths and filename suffix     *#
    #*#########################################*#
    # make directory to hold output, if it doesn't already exist
    # (if parent directories don't exist, will make them too)
    # if directory already exists, print error message and abort

    output_path = "numpyFiles/"
    
    # If you want to add a subdirectory with a time stamp you could use
    #   now=datetime.datetime.now()
    #   output_path = "matrixFiles/test1.lev.%d.%d.%d.%d.%d/" %(
    #                  max_level,now.year,now.month,now.day,now.hour)
    
    if not os.path.isdir(output_path):
       os.umask(7) # corresponds to "chmod 770"
       os.makedirs(output_path)
       if verbose:       
              print("The new data is in is "+output_path)
    else:
       while True:
          print("Output directory "+output_path+" already exists.")
          user_input=input("Do you want to overwrite this directory (y/n)?")
          if user_input in['y','n']:
             break
          else:
             print ("Not a valid input.")
       if user_input=='n':
          print ("Rename this directory and try again.")
          sys.exit()
       else:
          # delete old files so new reports are not appended to old ones
          files=glob.glob(output_path+'/sycamore*')
          for f in files:
             os.remove(f)          



    # See if there is a Sycamore configuration file.
    try:
        with open("sycamore_config.json", "r") as config_fp:
            config_data = json.load(config_fp)
    except:
        config_data = None

    if config_data == None:
        print("No Sycamore config file found.")
        sys.exit()

    print(config_data)
    print()
    print("Processing Sycamore configuration file.")
    print("  Number of rows    = {0}".format(config_data["num_rows"]))
    print("  Number of columns = {0}".format(config_data["num_cols"]))
    print("  Number of qubits  = {0}".format(config_data["num_qubits"]))
    print("  Number of cycles  = {0}".format(config_data["num_cycles"]))
    print("  Layout = {0}".format(config_data["layout"]))

    print()
    print("Building atomic Sycamore gates.")
    # One qubit sqrt(X) gate.
    one_x = np.array([[0.5+0.5j, 0.5-0.5j], [0.5-0.5j, 0.5+0.5j]], dtype=np.complex128)
    # One qubit sqrt(Y) gate.
    one_y = np.array([[0.5+0.5j, -0.5-0.5j], [0.5+0.5j, 0.5+0.5j]], dtype=np.complex128)
    # One qubit sqrt(W) gate.
    root2 = np.sqrt(np.complex128(2))
    one_w = np.array([[0.5+0.5j, -0.5j], [0.5, 0.5+0.5j]], dtype=np.complex128)
    one_w[0,1] *= root2
    one_w[1,0] *= root2
    print("Sqrt(X) gate is:")
    print(one_x)
    print("Sqrt(Y) gate is:")
    print(one_y)
    print("Sqrt(W) gate is:")
    print(one_w)

    # Two qubit gate as sum of tensor products.
    two_1a = np.array([[1, 0], [0, 0]], dtype=np.complex128)
    two_1b = np.array([[0, 1], [0, 0]], dtype=np.complex128)
    two_1c = np.array([[0, 0], [1, 0]], dtype=np.complex128)
    two_1d = np.array([[0, 0], [0, 1]], dtype=np.complex128)
    two_2a = np.array([[1, 0], [0, 0]], dtype=np.complex128)
    two_2b = np.array([[0, 0], [1j, 0]], dtype=np.complex128)
    two_2c = np.array([[0, 1j], [0, 0]], dtype=np.complex128)
    two_2d = np.array([[0, 0], [0, 0]], dtype=np.complex128)
    two_2d[1,1] = (np.sqrt(np.complex128(3)) + np.complex128(1j)) / np.complex128(2)
    print("Two qubit gate is:")
    print(np.kron(two_1a, two_2a) + np.kron(two_1b, two_2b) + np.kron(two_1c, two_2c) + np.kron(two_1d, two_2d))
    
    print()
    print("Building two qubit system matrices for each link type.")
    system_size = config_data["num_qubits"]
    matrix_size = 2**system_size
    two_qubit_matrices = {}
    for link_type, pair_list in config_data["linkinfo"].items():
        print("Working on link type {0}.".format(link_type))
        link_matrix = np.identity(matrix_size, dtype=np.complex128)
        for qubit1, qubit2 in pair_list:
            part_matrix_a = np.identity(1, dtype=np.complex128)
            part_matrix_b = np.identity(1, dtype=np.complex128)
            part_matrix_c = np.identity(1, dtype=np.complex128)
            part_matrix_d = np.identity(1, dtype=np.complex128)
            for i in range(system_size-1, -1, -1):
                if i == qubit1:
                    part_matrix_a = np.kron(two_1a, part_matrix_a)
                    part_matrix_b = np.kron(two_1b, part_matrix_b)
                    part_matrix_c = np.kron(two_1c, part_matrix_c)
                    part_matrix_d = np.kron(two_1d, part_matrix_d)
                elif i == qubit2:
                    part_matrix_a = np.kron(two_2a, part_matrix_a)
                    part_matrix_b = np.kron(two_2b, part_matrix_b)
                    part_matrix_c = np.kron(two_2c, part_matrix_c)
                    part_matrix_d = np.kron(two_2d, part_matrix_d)
                else:
                    part_matrix_a = np.kron(np.identity(2, dtype=np.complex128), part_matrix_a)
                    part_matrix_b = np.kron(np.identity(2, dtype=np.complex128), part_matrix_b)
                    part_matrix_c = np.kron(np.identity(2, dtype=np.complex128), part_matrix_c)
                    part_matrix_d = np.kron(np.identity(2, dtype=np.complex128), part_matrix_d)
            next_matrix = part_matrix_a + part_matrix_b + part_matrix_c + part_matrix_d
            link_matrix = np.matmul(link_matrix, next_matrix)
        two_qubit_matrices[link_type] = link_matrix
        out_file_name = output_path + "sycamore_2gate_matrix_{0}.csv".format(link_type)
        # print("Writing link type {0} matrix to file '{1}'.".format(link_type, out_file_name))
        # np.savetxt(out_file_name, link_matrix, delimiter=",")
    part_matrix_a = None
    part_matrix_b = None
    part_matrix_c = None
    part_matrix_d = None
    next_matrix = None
    link_matrix = None

    print()
    print("Building individual cycle matrices.")
    cycle_matrix_list = []
    for k, (gate_string, link_type) in enumerate(config_data["circuit"]):
        print("Working on cycle # {0}.".format(k))
        first_matrix = np.identity(1, dtype=np.complex128)
        for i in range(system_size-1, -1, -1):
            if gate_string[i] == "X":
                first_matrix = np.kron(one_x, first_matrix)
            elif gate_string[i] == "Y":
                first_matrix = np.kron(one_y, first_matrix)
            elif gate_string[i] == "W":
                first_matrix = np.kron(one_w, first_matrix)
            else:
                print("BAD SINGLE QUBIT GATE: k={0} i={1} gate_string=\"{2}\".".format(k, i, gate_string))
                sys.exit(1)


        this_cycle_matrix = np.matmul(first_matrix, two_qubit_matrices[link_type])
        cycle_matrix_list.append(this_cycle_matrix)
        out_file_name = output_path + "sycamore_cycle_{0}_matrix.csv".format(k)
        # print("Writing cycle {0} matrix to file '{1}'.".format(k, out_file_name))
        # np.savetxt(out_file_name, this_cycle_matrix, delimiter=",")
    first_matrix = None
    this_cycle_matrix = None

    print()
    print("Computing matrix for entire circuit.")
    circuit_matrix = np.identity(matrix_size, dtype=np.complex128)
    for k, cycle_matrix in enumerate(cycle_matrix_list):
        print("Multiplying cycle matrix {0} into circuit...".format(k))
        circuit_matrix = np.matmul(circuit_matrix, cycle_matrix)
    out_file_name = output_path + "sycamore_full_circuit_matrix.csv"
    print("Writing full circuit matrix to file '{0}'.".format(out_file_name))
    np.savetxt(out_file_name, circuit_matrix, delimiter=",")

    if system_size > 6:
        sys.exit()
    out_file_name = output_path + "sycamore_row_by_row_dump.txt"
    print("Dumping full circuit matrix row by row to file '{0}'.".format(out_file_name))
    with open(out_file_name, "w") as fp:
        for i in range(2**system_size):
            for j in range(2**system_size):
                val = circuit_matrix[i,j]
                print("{0},{1} -> {2}".format(i, j, val), file=fp)


