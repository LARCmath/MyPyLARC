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
 ##################################################################


from __future__ import print_function

import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as myp
from ctypes import *
import json


if __name__ == '__main__':


    ################################################################
    ##   Get the name of the parameter file that you are writing  ##
    ##   e.g. 'eigen.init_params'                                 ##
    ################################################################
    
    if 2 == len(sys.argv):
        param_file_name = sys.argv[1] 

    else:
        print("\nThis program requires the name of the parameter")
        print("file you want to create, as an argument.  e.g.")
        print("  python3 write_init_param_file.py my_app.init_params")
        sys.exit(1)

    ######################################################################
    ##  Create Python dictionary for sets of initialization parameters  ##
    ##  for LARC Initialization of Matrix Store and Operation Stores    ##
    ######################################################################

    data = { }
    data['small'] = []
    data['small'].append({
        'matrix_exponent': 22,
        'op_exponent': 19, 
        'max_level': 8, 
        'rnd_sig_bits': 45,
        'trunc_to_zero_bits': 45,
        'report_interval_seconds': 180,
        'min_memGiB_required': 5,
        'verbose': 1
    })
    data['medium'] = []
    data['medium'].append({
        'matrix_exponent': 25,
        'op_exponent': 24, 
        'max_level': 10, 
        'rnd_sig_bits': 52,
        'trunc_to_zero_bits': 52,
        'report_interval_seconds': 360,
        'min_memGiB_required': 50,
        'verbose': 1
    })
    data['large'] = []
    data['large'].append({
        'matrix_exponent': 31,
        'op_exponent': 30, 
        'max_level': 32, 
        'rnd_sig_bits': 60,
        'trunc_to_zero_bits': 60,
        'report_interval_seconds': 3600,
        'min_memGiB_required': 500,
        'verbose': 1
    })
    data['desktop'] = []
    data['desktop'].append({
        'matrix_exponent': 22,
        'op_exponent': 19, 
        'max_level': 8, 
        'rnd_sig_bits': 54,
        'trunc_to_zero_bits': 54,
        'report_interval_seconds': 360,
        'min_memGiB_required': 50,
        'verbose': 1
    })
    data['server'] = []
    data['server'].append({
        'matrix_exponent': 31,
        'op_exponent': 30, 
        'max_level': 32, 
        'rnd_sig_bits': -1,   # default value
        'trunc_to_zero_bits': -1,  #default value
        'report_interval_seconds': 3600,
        'min_memGiB_required': 500,
        'verbose': 1
    })
    data['debug'] = []
    data['debug'].append({
        'matrix_exponent': 12,
        'op_exponent': 10, 
        'max_level': 5, 
        'rnd_sig_bits': 45,   # default value
        'trunc_to_zero_bits': 45,   # default value
        'report_interval_seconds': 180,
        'min_memGiB_required': 5,
        'verbose': 1
    })
    
    with open(param_file_name,'w') as outfile:
        json.dump(data,outfile,indent=" ")
        

    userInput = input("Would you like to see some of the parameters (y/n)?")
    if (userInput == 'n'):
        sys.exit(1)

        
    #### try to read a parameter file into a dictionary
    with open(param_file_name,'r') as init_params:
        init_param = json.load(init_params)
        print('The initialization parameters for running on a desktop are:')
        for param in init_param['desktop']:
            print('MatrixExponent: %d' %(param['matrix_exponent']))
            print('OpExponent: %d' %(param['op_exponent']))
            print('MaxLevel: %d' %(param['max_level']))
            print('RoundSigBits: %d' %(param['rnd_sig_bits']))
            print('TruncToZeroBits: %d' %(param['trunc_to_zero_bits']))
            print('ReportIntervalSeconds: %d' %(param['report_interval_seconds']))
            print('Minimum Memory Available in GiB Required: %d' %(param['min_memGiB_required']))
            print('Verbose: %d' %(param['verbose']))
            print('')

        print('The initialization parameters for running on a medium memory machine are:')
        for param in init_param['medium']:
            print('MatrixExponent: %d' %(param['matrix_exponent']))
            print('OpExponent: %d' %(param['op_exponent']))
            print('MaxLevel: %d' %(param['max_level']))
            print('RoundSigBits: %d' %(param['rnd_sig_bits']))
            print('TruncToZeroBits: %d' %(param['trunc_to_zero_bits']))
            print('ReportIntervalSeconds: %d' %(param['report_interval_seconds']))
            print('Minimum Memory Available in GiB Required: %d' %(param['min_memGiB_required']))
            print('Verbose: %d' %(param['verbose']))
            print('')

        print('The initialization parameters for running on a large memory machine are:')
        for param in init_param['large']:
            print('MatrixExponent: %d' %(param['matrix_exponent']))
            print('OpExponent: %d' %(param['op_exponent']))
            print('MaxLevel: %d' %(param['max_level']))
            print('RoundSigBits: %d' %(param['rnd_sig_bits']))
            print('TruncToZeroBits: %d' %(param['trunc_to_zero_bits']))
            print('ReportIntervalSeconds: %d' %(param['report_interval_seconds']))
            print('Minimum Memory Available in GiB Required: %d' %(param['min_memGiB_required']))
            print('Verbose: %d' %(param['verbose']))
            print('')

