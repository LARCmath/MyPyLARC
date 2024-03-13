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
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *
import datetime
import json

## \file sycamore.py
#
#  \brief Demonstrates LARC on Sycamore quantum supremacy computation.
#
if __name__ == '__main__':

    verbose = 0
    debug = 0

    #*################################*#
    #*    Basic Parameter Setting     *#
    #*################################*#

    #*##################################################################*#
    # Figure out the scalarType                                          #
    # In the Makefile you can compile with different scalarType values   #
    # Define string for using in formating filenames                     #
    #*##################################################################*#
    scalarTypeStr = mypy.cvar.scalarTypeStr

    #*##################################################*#
    #*   Find out if machine has a large amount of      *#
    #*   memory available so we can make bigger tables  *#
    #*##################################################*#
    memory_available = mypy.memory_available_GiB()
    if (verbose > 0):
        print("\nThe memory available is %ld GiB" %memory_available)
        print("We will use this to select which computing_env to read from parameter file.")
        print("You could write code to select computing_env automatically.")

    if (memory_available > 200):
        if (verbose > 0):
            print("\nThis memory is more than 200 GiB\n")
        computing_env = 'large'
    else:    
        if (memory_available > 50):
            if (verbose > 0):
                print("\nThis memory is between 50 and 200 GiB\n")
            computing_env = 'medium'
        else:
            if (verbose > 0):
                print("\nThis memory is less than 50 GiB\n")
            computing_env = 'small'
            
    if (verbose > 0):
        print("This program believes the computing_environment is %s" %computing_env)

    #*###################################*#
    #*    Print baseline usage report    *#
    #*###################################*#
    if (verbose > 0):
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mypy.memory_and_time_report(0, "stdout")

    # read the parameter file into a python dictionary
    #with open('../InitParams/REPLACE_THIS.init_params','r') as init_file:
    with open('../InitParams/tutorial.init_params','r') as init_file:
        init_param = json.load(init_file)
        for p in init_param[computing_env]:
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
            max_level= 12 #p['max_level']
            regionbitparam = p['regionbitparam']
            zeroregionbitparam = p['zeroregionbitparam']
            report_interval_seconds = p['report_interval_seconds']
            min_memGiB_required = p['min_memGiB_required']
            p_verbose = p['verbose']

    # warn if the commandline value for verbose differs from the parameter file value for verbose
    if (verbose > 0):
        if (verbose != p_verbose):
            print("NOTE: This program uses commandline (verbose = %d) " %verbose)
            print("      rather than the parameter file (verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    # initialize LARC
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)
    mypy.create_report_thread(report_interval_seconds)
    if (verbose):
        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("StopHogging check to see if program is too large to occur once every %d seconds.\n" %report_interval_seconds)

    half_level = max_level//2     


    
    print("\n    ######################################################################")
    print("    ##  Create matrices for Sycamore quantum supremacy computation.     ##" )
    print("    ######################################################################\n")
    

    #*#########################################*#
    #*    Output paths and filename suffix     *#
    #*#########################################*#
    # make directory to hold output, if it doesn't already exist
    # (if parent directories don't exist, will make them too)
    # if directory already exists, print error message and abort

    output_path = "matrixFiles/"
    
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


                  
    if scalarTypeStr not in ('Complex', 'MPComplex', 'MPRatComplex', 'Clifford'):
        raise Exception('Recompile with complex scalarType, not %s.'
                        %scalarTypeStr)


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
    print("Building two qubit system matrices for each link type.")
    system_size = config_data["num_qubits"]
    two_qubit_matrices = {}
    for link_type, pair_list in config_data["linkinfo"].items():
        link_matrixID = mypy.get_identity_pID(system_size)
        for qubit1, qubit2 in pair_list:
            next_matrixID = mypy.build_sycamore_2gate(qubit1, qubit2, system_size)
            link_matrixID = mypy.matrix_mult(link_matrixID, next_matrixID)
        two_qubit_matrices[link_type] = link_matrixID
        out_file_name = output_path + "sycamore_2gate_matrix_{0}.json".format(link_type)
        # print("Writing link type {0} matrix to file '{1}'.".format(link_type, out_file_name))
        mypy.fprint_larcMatrixFile(link_matrixID, out_file_name)

    print()
    print("Building individual cycle matrices.")
    cycle_matrix_list = []
    for k, (gate_string, link_type) in enumerate(config_data["circuit"]):
        gate_sequence = ["sqrt" + gate_name for gate_name in gate_string]
        first_matrixID = mypy.build_sycamore_gate_sequence(gate_sequence, system_size)
        this_cycle_matrixID = mypy.matrix_mult(first_matrixID, two_qubit_matrices[link_type])
        cycle_matrix_list.append(this_cycle_matrixID)
        out_file_name = output_path + "sycamore_cycle_{0}_matrix.json".format(k)
        # print("Writing cycle {0} matrix to file '{1}'.".format(k, out_file_name))
        mypy.fprint_larcMatrixFile(this_cycle_matrixID, out_file_name)

    print()
    print("Computing matrix for entire circuit.")
    circuit_matrixID = mypy.get_identity_pID(system_size)
    for k, cycle_matrixID in enumerate(cycle_matrix_list):
        print("Multiplying cycle matrix {0} into circuit...".format(k))
        circuit_matrixID = mypy.matrix_mult(circuit_matrixID, cycle_matrixID)
    out_file_name = output_path + "sycamore_full_circuit_matrix.json"
    # print("Writing full circuit matrix to file '{0}'.".format(out_file_name))
    mypy.fprint_larcMatrixFile(circuit_matrixID, out_file_name)

    if system_size > 6:
        sys.exit()
    out_file_name = output_path + "sycamore_row_by_row_dump.txt"
    print("Dumping full circuit matrix row by row to file '{0}'.".format(out_file_name))
    with open(out_file_name, "w") as fp:
        for i in range(2**system_size):
            for j in range(2**system_size):
                val = mypy.get_readableString_scalar_from_pID_and_coords(circuit_matrixID, i, j)
                print("{0},{1} -> {2}".format(i, j, val), file=fp)


