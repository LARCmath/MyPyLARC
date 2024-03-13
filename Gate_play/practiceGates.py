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

## \file practiceGates.py
#
#  \brief Demonstrates creation of reversible circuits.
#
# This program creates some matrices using routines from MyPyLARC/src/gates.c.
# It is a fairly simple demonstration of how to create a reversible circuit
# using the universal gate library for reversible computing
# (NOT, CNOT, CCNOT) and calculating results from such a circuit.
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
            max_level= p['max_level']
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

    half_level = max_level//2     

    if (verbose):
        # create a report thread
        mypy.create_report_thread(report_interval_seconds)
        # print out all sorts of information
        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("StopHogging check to see if program is too large to occur once every %d seconds.\n" %report_interval_seconds)
        print("\n    ######################################################################")
        print("    ##  Create some matrices using routines from MyPyLARC/src/gates.c   ##")
        print("    ##  The following routine creates a reversible circuit.             ##")
        print("    ##  The circuit consists of control wires and target wires          ##")
        print("    ##  the gates connect the controls and targets with CNots or CCNots ##")
        print("    ##  All gates commute in this circuit.                              ##")
        print("    ##  To show various LARC capabilities we do the computation twice   ##")
        print("    ##  with different ordering of the gates saving the result as json  ##")
        print("    ##  compressed matrix files.  We clean the matrix store between     ##")
        print("    ##  the two calculations so that the json files are not identical.  ##")
        print("    ##  The user can diff or look at the two compressed files to see    ##")
        print("    ##  that the internal matrixIDs in the saved files are not the same ##")
        print("    ##  although the matrices themselves are the same and will have the ##")
        print("    ##  same matrix ID when read into the same active matrix store.     ##")
        print("    ##  We show this by reading in the original file and checking that  ##")
        print("    ##  the matrixIDs are the same.                                     ##")
        print("    ######################################################################\n")

        print("**************************************************")
        print("*  Running a test problem with max_level %d      *" %(max_level))
        print("**************************************************\n")

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
          files=glob.glob(output_path+'/circuit*')
          for f in files:
             os.remove(f)          


    # create base suffix for file names
    file_suffix = "_gates%d" %max_level
    if verbose:      
        print("The file_suffix is "+file_suffix)
        print("\n")

    # Get a large identity matrix
    Imax_mID = mypy.get_identity_pID(max_level)

    # pick offsets for CCnot and Control Not Gates
    CntrlNot_offset = 2  # offset for the Control Not gate
    CCnotA_offset = 3  # offset for CCnot gate A
    CCnotB_offset = 4  # offset for CCnot gate B
              
    circuit_matrix1_mID = Imax_mID
    for i in range(half_level):
        # CCnot gates, with controls on control wires circularly shifted
        # by CCnotB_offset and CCnotA_offset
        circuit_matrix1_mID = mypy.matrix_mult(circuit_matrix1_mID,
           mypy.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
        # Cnot gates with control on control wires circularly shifted by CntrlNot_offset
        circuit_matrix1_mID = mypy.matrix_mult(circuit_matrix1_mID,
           mypy.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
    base_name = "circuit_matrix1"+file_suffix+".json"
    out_name1 = output_path+base_name
    mypy.fprint_larcMatrixFile(circuit_matrix1_mID,out_name1)
    print("The circuit 1 has been stored as a compressed matrix in the file")
    print(out_name1)
    
    print("\nWe are about to clean out the matrix store. Before we do,")
    print("we verify that the matrix is present in the store.")
    if (mypy.matrix_is_invalid(circuit_matrix1_mID)):
        print("The matrix ID for circuit 1 is no longer in the matrix store.")
    else:
        print("The matrix ID for circuit 1 is in the matrix store.")

    # Clean out the matrix store of everything but the basic preloaded matrices    
    mypy.clean_matrix_storage()

    print("\nWe have cleaned the matrix store, and now need to verify the matrix is gone.")
    if (mypy.matrix_is_invalid(circuit_matrix1_mID)):
        print("The matrix ID for circuit 1 is no longer in the matrix store.")
    else:
        print("The matrix ID for circuit 1 IS IN THE MATRIX STORE.")
        print("continuing on despite this...")
          
    circuit_matrix2_mID = Imax_mID
    for i in range(half_level):
        # Cnot gates with control on control wires circularly shifted by CntrlNot_offset
        circuit_matrix2_mID = mypy.matrix_mult(circuit_matrix2_mID,
           mypy.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
        # CCnot gates, with controls on control wires circularly shifted
        # by CCnotB_offset and CCnotA_offset
        circuit_matrix2_mID = mypy.matrix_mult(circuit_matrix2_mID,
           mypy.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
    base_name = "circuit_matrix2"+file_suffix+".json"
    out_name2 = output_path+base_name
    mypy.fprint_larcMatrixFile(circuit_matrix2_mID,out_name2)
    print("The circuit 2 has been stored as a compressed matrix in the file")
    print(out_name2)

    print("\nThe user might like to try the commandline to check whether these") 
    print("two compressed matrices represent the same matrix:")
    print("python ../larc/src/python/larcMatrix_files_matrix_equal.py  %s %s" %(out_name1,out_name2))
          

    # Second way to verify that two circuits produce identical matrices
    # Read in the first compressed matrix file and get its assigned matrixID
    circuit_matrix3_mID = mypy.read_larcMatrixFile(out_name1)
    print("\nA second way to verify equivalence is to read in our first stored json file") 
    print("and check to see if the read in matrix is assigned the same matrix ID as the") 
    print("matrix we have just calculated.")
    print("The second circuit has not been cleaned and its matrix ID is %d" %circuit_matrix2_mID)
    print("The matrix ID for the read in matrix from the first circuit is %d" %circuit_matrix3_mID)
    if (circuit_matrix3_mID == circuit_matrix2_mID):
        print("We see the matrixIDs are equal because the matrices are identical")
        print("which shows that the two circuits have the same effect.")
    else:
        print("We are confused, check the code for mistakes.")
