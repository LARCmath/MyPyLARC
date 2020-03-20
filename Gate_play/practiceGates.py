#!/usr/bin/env python3

 ##################################################################
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
 ##################################################################


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

if __name__ == '__main__':

    verbose = 0
    debug = 0

    ##################################################################
    ## Exit if the scalarType is not compatible with this program.  ##
    ##################################################################
    if mypy.cvar.scalarTypeDef != "r" and mypy.cvar.scalarTypeDef != "q":
        print("This program doesn't run with all scalarTypes.\n")
        print("Remake with TYPE=REAL or TYPE=MPRATIONAL.\n")
        sys.exit()

    ####################################
    ##    Basic Parameter Setting     ##
    ####################################

    ####################################################
    ##   Find out if machine is desktop workstation   ##
    ##   or a CPU-cycle servers (cs1-cs6)             ##
    ####################################################
    machine = os.uname()[1]
    cs = 0        # on desktop workstation, with smaller memory
    computing_env = 'desktop'
    if (machine.find('cs') >= 0):
        cs = 1    # on CPU-cycle server cs1-cs6, with larger memory
        computing_env = 'server'
        if (verbose):
            print("This machine is a CPU-cycle server")
    else:
        if (verbose):
            print("This machine is a desktop work station")

    if (debug):
        computing_env = 'debug'

    #######################################
    ##    Print baseline usage report    ##
    #######################################
    mypy.rusage_report(0, "stdout")

    #### read a parameter file into a dictionary
    with open('../InitParams/gate.init_params','r') as init_file:
        init_param = json.load(init_file)
        for p in init_param[computing_env]:
            if (verbose):
                print('MatrixExponent: %d' %(p['matrix_exponent']))
                print('OpExponent: %d' %(p['op_exponent']))
                print('MaxLevel: %d' %(p['max_level']))
                print('RoundSigBits: %d' %(p['rnd_sig_bits']))
                print('TruncToZeroBits: %d' %(p['trunc_to_zero_bits']))
                print('ReportIntervalSecs: %d' %(p['report_interval_seconds']))
                print('Verbose: %d' %(p['verbose']))
                print('')
            matrix_exponent = p['matrix_exponent']
            op_exponent = p['op_exponent']
            max_level= p['max_level']
            rnd_sig_bits = p['rnd_sig_bits']
            trunc_to_zero_bits = p['trunc_to_zero_bits']
            report_interval_seconds = p['report_interval_seconds']
            larc_verbose = p['verbose']
        if (larc_verbose==0 and verbose==1):
            larc_verbose=1
            
    print_naive = 0
    print_nonzeros = 0
    if computing_env == 'debug':
        if (verbose):
            print("Using small parameters for debugging")
            if print_naive:
                print("  will print files of naive matrices")
            else: 
                print("  not printing files of naive matrices")
                if print_nonzeros:
                    print("  will print files of nonzero matrices\n")
                else: 
                    print("  not printing files of nonzero matrices\n")

    elif computing_env == 'desktop':
        if (verbose):
            print("Problem size is small enough to run on desktop")
            if print_naive:
                print("  will print files of naive matrices")
            else: 
                print("  not printing files of naive matrices")
                if print_nonzeros:
                    print("  will print files of nonzero matrices\n")
                else: 
                    print("  not printing files of nonzero matrices\n")

    ## LARGE STORES for cs1l,cs4l,cs9l
    elif computing_env == 'server':
        if (verbose):
            print("Problem size is NOT small enough to run on desktop")
            if print_naive:
                print("  WARNING: will try to print files of naive matrices!!!")
            else: 
                print("  not printing files of naive matrices")
                if print_nonzeros:
                    print("  WARNING: will print files of nonzero matrices!!!\n")
                else: 
                    print("  not printing files of nonzero matrices\n")
    else:
        print("Error: unexpected value for 'computing_env' - exiting")
        exit(1)

    print("larc_verbose is %d" %larc_verbose)
    mypy.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,larc_verbose)
    mypy.create_report_thread(report_interval_seconds)
    if (verbose):
        print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
        print("Seppuku check to see if program is too large to occur once every %d seconds.\n" %report_interval_seconds)

    half_level = max_level//2     


    
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
    

    if verbose:
        print("**************************************************")
        print("*  Running a test problem with max_level %d      *" %(max_level))
        print("**************************************************\n")

    #############################################
    ##    Output paths and filename suffix     ##
    #############################################
    ## make directory to hold output, if it doesn't already exist
    ## (if parent directories don't exist, will make them too)
    ## if directory already exists, print error message and abort

    output_path = "matrixFiles/"
    
    ## If you want to add a subdirectory with a time stamp you could use
    #   now=datetime.datetime.now()
    #   output_path = "matrixFiles/test1.lev.%d.%d.%d.%d.%d/" %(
    #                  max_level,now.year,now.month,now.day,now.hour)
    
    if not os.path.isdir(output_path):
       os.umask(7) ## corresponds to "chmod 770"
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
          ## delete old files so new reports are not appended to old ones
          files=glob.glob(output_path+'/*')
          for f in files:
             os.remove(f)          


    ## create base suffix for file names
    file_suffix = "_gates%d" %max_level
    if verbose:      
        print("The file_suffix is "+file_suffix)
        print("\n")

    ## Get a large identity matrix
    Imax_mID = mypy.get_identity_matrixID(max_level)

    ## pick offsets for CCnot and Control Not Gates
    CntrlNot_offset = 2  ## offset for the Control Not gate
    CCnotA_offset = 3  ## offset for CCnot gate A
    CCnotB_offset = 4  ## offset for CCnot gate B
              
    circuit_matrix1_mID = Imax_mID
    for i in range(half_level):
        ## CCnot gates, with controls on control wires circularly shifted
        ## by CCnotB_offset and CCnotA_offset
        circuit_matrix1_mID = mypy.matrix_mult_matrixID(circuit_matrix1_mID,
           mypy.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
        ## Cnot gates with control on control wires circularly shifted by CntrlNot_offset
        circuit_matrix1_mID = mypy.matrix_mult_matrixID(circuit_matrix1_mID,
           mypy.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
    base_name = "circuit_matrix1"+file_suffix+".json"
    out_name1 = output_path+base_name
    mypy.write_larcMatrix_file_by_matID(circuit_matrix1_mID,out_name1)
    print("The circuit 1 has been stored as a compressed matrix in the file")
    print(out_name1)
    
    print("\nWe are about to clean out the matrix store. Before we do,")
    print("we verify that the matrix is present in the store.")
    if (mypy.matrix_is_invalid_matrixID(circuit_matrix1_mID)):
        print("The matrix ID for circuit 1 is no longer in the matrix store.")
    else:
        print("The matrix ID for circuit 1 is in the matrix store.")

    # Clean out the matrix store of everything but the basic preloaded matrices    
    mypy.clean_matrix_store()

    print("\nWe have cleaned the matrix store and verify the matrix is gone.")
    if (mypy.matrix_is_invalid_matrixID(circuit_matrix1_mID)):
        print("The matrix ID for circuit 1 is no longer in the matrix store.")
    else:
        print("The matrix ID for circuit 1 is in the matrix store.")
          
    circuit_matrix2_mID = Imax_mID
    for i in range(half_level):
        ## Cnot gates with control on control wires circularly shifted by CntrlNot_offset
        circuit_matrix2_mID = mypy.matrix_mult_matrixID(circuit_matrix2_mID,
           mypy.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
        ## CCnot gates, with controls on control wires circularly shifted
        ## by CCnotB_offset and CCnotA_offset
        circuit_matrix2_mID = mypy.matrix_mult_matrixID(circuit_matrix2_mID,
           mypy.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
    base_name = "circuit_matrix2"+file_suffix+".json"
    out_name2 = output_path+base_name
    mypy.write_larcMatrix_file_by_matID(circuit_matrix2_mID,out_name2)
    print("The circuit 2 has been stored as a compressed matrix in the file")
    print(out_name2)

    print("\nThe user might like to try the commandline to check whether these") 
    print("two compressed matrices represent the same matrix:")
    print("../larc/src/python/json_files_matrix_equal.py  %s %s" %(out_name1,out_name2))
          

    # Second way to verify that two circuits produce identical matrices
    # Read in the first compressed matrix file and get its assigned matrixID
    circuit_matrix3_mID = mypy.read_larcMatrix_file_return_matID(out_name1)
    print("\nA second way to verify equivalence is to read in our first stored json file") 
    print("and check to see if the read in matrix is assigned the same matrix ID as the") 
    print("matrix we have just calculated.")
    print("The second circuit has not been cleaned and its matrix ID is %d" %circuit_matrix2_mID)
    print("The matrix ID for the read in matrix from the first circuit is %d" %circuit_matrix3_mID)
    if (circuit_matrix3_mID == circuit_matrix2_mID):
        print("We see the matrixIDs are equal because the matrices are identical")
        print("which implies the two circuits have the same effect.")
    else:
        print("We are confused, check the code for mistakes.")
