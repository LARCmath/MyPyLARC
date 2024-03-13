#!/usr/bin/env python3

 ##################################################################
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
 ##################################################################


from __future__ import print_function, division

import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import MyPyLARC as mpl
import numpy as np
from ctypes import *
import datetime
import json

if __name__ == '__main__':

    verbose = 1
    debug = 0

    ##################################################################
    ## Exit if the scalarType is not compatible with this program.  ##
    ##################################################################
    if mpl.cvar.scalarTypeDef != "r" and mpl.cvar.scalarTypeDef != "q":
        print("This program doesn't run with all scalarTypes.\n")
        print("Remake with TYPE=REAL or TYPE=MPRATIONAL.\n")
        sys.exit()

    ######################################################
    ##   Find out if machine has a large amount of      ##
    ##   memory available so we can make bigger tables  ##
    ######################################################
    memory_available = mpl.memory_available_GiB()
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

    #######################################
    ##    Print baseline usage report    ##
    #######################################
    if (verbose > 0):
        print("")
        print("In the following baseline usage report")
        print("RSS, resident set size, refers to size of the process held in RAM.")
        print("HASHSTATS: hash occupancy means, variances and crash resolution chain lengths")
        mpl.memory_and_time_report(0, "stdout")

    ## read the parameter file into a python dictionary
    #with open('../InitParams/REPLACE_THIS.init_params','r') as init_file:
    with open('../../InitParams/power_method.init_params','r') as init_file:
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

    ## warn if the commandline value for verbose differs from the parameter file value for verbose        
    if (verbose > 0):
        if (verbose != p_verbose):
            print("NOTE: This program uses commandline (verbose = %d) " %verbose)
            print("      rather than the parameter file (verbose = %d)." %p_verbose)
            print("      The verbose key is:  0=SILENT, 1=BASIC, 2=CHATTY, 3=DEBUG.")

    ## initialize LARC
    mpl.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)

    mpl.create_report_thread(report_interval_seconds)
    print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
    print("stopHogging check to see if program is too large to occur once every %d seconds.\n" %report_interval_seconds)

    half_level = max_level//2     

    if verbose:
        print("**************************************************")
        print("*  Running the Test problem with max_level %d    *" %(max_level))
        print("**************************************************\n")

    #############################################
    ##    Output paths and filename suffix     ##
    #############################################
    ## make directory to hold output, if it doesn't already exist
    ## (if parent directories don't exist, will make them too)
    ## if directory already exists, print error message and abort
    now=datetime.datetime.now()
    output_path = "/scr/eigen/test1.lev.%d.%d.%d.%d.%d/" %(max_level,now.year,now.month,now.day,now.hour)
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
    file_suffix = "_Eig%d" %max_level
    if verbose:      
        print("The file_suffix is "+file_suffix)
        print("\n")


    ###############################################
    ##    Create some matrices                   ##
    ##  using routines from MyPyLARC/src/gates.c ##
    ###############################################
    ## Get convenient name for the largest identity matrix
    Imax_mID = mpl.get_identity_pID(max_level)


    ## pick offsets for CCnot and Control Not Gates
    CntrlNot_offset = 2  ## offset for the Control Not gate
    CCnotA_offset = 3  ## offset for CCnot gate A
    CCnotB_offset = 4  ## offset for CCnot gate B
              

    first_try_mID = Imax_mID
    for i in range(half_level):
        ## CCnot gates, with controls on first half circularly shifted
        ## by CCnotB_offset and CCnotA_offset
        first_try_mID = mpl.matrix_mult(first_try_mID,
           mpl.build_ccnot_gate((i+CCnotB_offset)%half_level,
           (i+CCnotA_offset)%half_level,i+half_level))
        ## CX gates, with control on first word circularly shifted by CntrlNot_offset
        first_try_mID = mpl.matrix_mult(first_try_mID,
           mpl.build_cnot_gate((i+CntrlNot_offset)%half_level,i+half_level,0))
    base_name = "first_try"+file_suffix+".json"
    out_name = output_path+base_name
    mpl.write_larcMatrix_file_by_matID(first_try_mID,out_name)
          
              
