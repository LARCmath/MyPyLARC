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

import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src"))
import MyPyLARC as mypy
from ctypes import *
import json

## \file write_init_param_file.py
#
# \brief This program creates a default parameter file for a project.
#
# The user may then edit the file to change parameters if desired.
#

if __name__ == '__main__':


    #*############################################################*#
    #*   Get the name of the parameter file that you are writing  *#
    #*   e.g. 'eigen.init_params'                                 *#
    #*############################################################*#
    
    if 2 == len(sys.argv):
        json_file_name = sys.argv[1] 

    else:
        print("\nThis program requires the name of the json")
        print("file that you will read into the Python dictionary.  e.g.")
        print("  python3 read_json_to_Python_dictionary.py <path to file>")
        print("\nExample:")
        print("  python3 read_json_to_Python_dictionary.py JsonFiles/In/int_lev3.json")
        sys.exit(1)



    #*##################################################################*#
    #*  Create Python dictionary for mock json file                     *#
    #*##################################################################*#
    #*
    #* In Python:
    #*   []:  list object or an index subscript, my_List[x].
    #*   {}:  dictionary object.
    #*          reference values in a dictionary the same way we
    #*          reference list, with subscripts. my_dict["value"]
    #*   a_list = ['one', 'two', 3, 4]
    #*   a_dict = { one: 1, two: 2 }
    #*
    #*   append is a list action  my_list.append
    #*
    #*  For example in the InitParam writer, we start by declaring
    #*  a dictionary called data, then give an item in the dictionary
    #*  which will be a list called 'small'. We append an item to the
    #*  list which is a nested dictionary and is given no name:
    #*    data = { }
    #*    data['small'] = []
    #*    data['small'].append({
    #*        'matrix_exponent': 22,
    #*        'op_exponent': 19, 
    #*        'max_level': 8, 
    #*
    #*  When this is written into a json file it is written as:
    #*     {
    #*  "small": [
    #*    {
    #*   "matrix_exponent": 22,
    #*   "op_exponent": 19,
    #*    ...
    #*   "verbose": 1
    #*  }
    #* ],
    #*
    #* Whereas the notation in a json file is to use
    #* The curly braces to define an object or list of objects
    #*     each object has a pair consisting of a
    #*           thekey:  followed by an object, list, or value
    #*     and each object in a list of objects is separated by a comma
    #*     e.g. {'amount':5, 'name':'sally'}
    #* The square braces are used to define a sequence of
    #*    objects, values, or lists
    #*     e.g. [2, 3, 'john', 7]
    #* Here is a mixed example:
    #* [{'value':6}, {'favorites':
    #* [1,3,3, {'first letter':'a', 'last digit':9}] }, 77]
    #*
    #*  In python a dictionary a_dict {one: 1, name: "cindy"}
    #*  translates in json to an object list {"one":1, "name":"cindy"}
    #*
    #*  In python a list a_list ['one','two',3,4]
    #*  translates in json to a list ["one","two",3,4]
    #*
    
    # LARCdata = {
    #     "matrixID_max":187,
    #     'matid':186,
    #     "table": {
    #         121:[0, 0, '1'],
    #         0:[0, 0, "0"],
    #         176:[1, 1, 121, 0, 0, 0],
    #         177:[1, 1, 0, 0, 121, 0],
    #         178:[1, 1, 0, 121, 0, 0],
    #         179:[1, 1, 0, 0, 0, 121],
    #         180:[2, 2, 176, 177, 178, 179],
    #         21:[1, 1, 0, 0, 0, 0],
    #         32:[2, 2, 21, 21, 21, 21],
    #         186:[3, 3, 180, 32, 32, 180],
    #         "end":0 }
    # }

    # print(LARCdata)
     

    # with open(json_file_name,'w') as outfile:
    #     json.dump(LARCdata,outfile,indent=" ")
        

    # sys.exit()

    
    #*** try to read a parameter file into a dictionary
    with open(json_file_name,'r') as in_file:
        LARC_matrix = json.load(in_file)

    print("Read in json file for compressed matrix %s" %json_file_name)
    print(LARC_matrix)
    print()

    print("The dictionary LARC_matrix['table'] is:")
    print(LARC_matrix['table'])
    print()

    # print("The value associated with key LARC_matrix['table']['178'] is:")
    # print(LARC_matrix['table']['179'])
    # print()

    print("The key:value pairs in the LARC_matrix['table'] are:")
    print(LARC_matrix['table'].items())
    print()

    items = LARC_matrix['table'].items()

    print("What does the list of items look like?")
    print(items)
    print()

    print("What if we list all the key value pairs?")
    for key, value in LARC_matrix['table'].items():
        print("  %s: %s" %(key,value))


    #  path_log_file =   path_log_file = path_log_dir + "/file2"
    # with open(path_log_file,'w') as outfile:
    #     print("This should go to a file",file=outfile)

    short_table = LARC_matrix['table']
    del short_table['end']

    temp_str = "abcd'efghi"
    print("\nTemp string is %s" %temp_str)
    temp_str = temp_str.replace('a','j')
    print("Temp string is %s" %temp_str)
    temp_str = temp_str.replace("'",'"')
    print("Temp string is %s\n" %temp_str)


    print("This is what a search for info in the dictionary found:")
    info_str = str(LARC_matrix.get('info','FAIL'))
    if (info_str == 'FAIL'):
        print("    *info was not found* ")
    else:
        print(info_str)
    
    
    # print("\nshort table is")
    # print(short_table)
    # print("\nLength of the short table is %d\n" %len(short_table))
    path_outfile = "JsonFiles/Out/copy4.json"
    with open(path_outfile,'w') as outfile:
        print("{",file=outfile)
        next_str = '  "matrixID_max":'+str(LARC_matrix['matrixID_max'])+','
        print(next_str,file=outfile)
        next_str = '  "matid":'+str(LARC_matrix['matid'])+','
        print(next_str,file=outfile)
        info_str = str(LARC_matrix.get('info','MISSING'))
        if (info_str != 'MISSING'):
            next_str = '  "info":{'
            for key, value in LARC_matrix['info'].items():
                print(next_str,file=outfile)
                next_str = '   "'+key+'":"'+value+'",'
            next_str = '      "end":"" },'
            print(next_str,file=outfile)
  #* This is where we put in the larc size using len(short_table)
        next_str = '  "table":{'
        print(next_str,file=outfile)
        for key, value in short_table.items():
            next_str = '    "'+key+'":'+str(value)+','
            next_str = next_str.replace("'",'"')
            print(next_str,file=outfile)
        next_str = '      "end":0 }'
        print(next_str,file=outfile)
        print("}",file=outfile)
        

  #  look in LAR/src/info_store.c and find out where we set the
  #  permitted types of entry.
  #  check with Steve and Mark to see if they like my name ideas.
            
    sys.exit()
    


    
    #*##################################################################*#
    #*  Create Python dictionary for sets of initialization parameters  *#
    #*  for LARC Initialization of Matrix Store and Operation Stores    *#
    #*##################################################################*#

    sys.exit()

    data = { }
    data['small'] = []
    data['small'].append({
        'matrix_exponent': 22,
        'op_exponent': 19, 
        'max_level': 8, 
        'regionbitparam': -1,
        'zeroregionbitparam': -1,
        'report_interval_seconds': 180,
        'min_memGiB_required': 5,
        'verbose': 1
    })
    data['medium'] = []
    data['medium'].append({
        'matrix_exponent': 25,
        'op_exponent': 24, 
        'max_level': 10, 
        'regionbitparam': -1,
        'zeroregionbitparam': -1,
        'report_interval_seconds': 360,
        'min_memGiB_required': 50,
        'verbose': 1
    })
    data['debug'] = []
    data['debug'].append({
        'matrix_exponent': 12,
        'op_exponent': 10, 
        'max_level': 5, 
        'regionbitparam': 45,   # default value
        'zeroregionbitparam': 45,   # default value
        'report_interval_seconds': 180,
        'min_memGiB_required': 5,
        'verbose': 1
    })
    
    with open(param_file_name,'w') as outfile:
        json.dump(data,outfile,indent=" ")
        

    userInput = input("Would you like to see some of the parameters (y/n)?")
    if (userInput == 'n'):
        sys.exit(1)

        
    #*** try to read a parameter file into a dictionary
    with open(param_file_name,'r') as init_params:
        init_param = json.load(init_params)
        print('The initialization parameters for running on a desktop are:')
        for param in init_param['desktop']:
            print('MatrixExponent: %d' %(param['matrix_exponent']))
            print('OpExponent: %d' %(param['op_exponent']))
            print('MaxLevel: %d' %(param['max_level']))
            print('regionbitparam: %d' %(param['regionbitparam']))
            print('zeroregionbitparam: %d' %(param['zeroregionbitparam']))
            print('ReportIntervalSeconds: %d' %(param['report_interval_seconds']))
            print('Minimum Memory Available in GiB Required: %d' %(param['min_memGiB_required']))
            print('Verbose: %d' %(param['verbose']))
            print('')

        print('The initialization parameters for running on a medium memory machine are:')
        for param in init_param['medium']:
            print('MatrixExponent: %d' %(param['matrix_exponent']))
            print('OpExponent: %d' %(param['op_exponent']))
            print('MaxLevel: %d' %(param['max_level']))
            print('regionbitparam: %d' %(param['regionbitparam']))
            print('zeroregionbitparam: %d' %(param['zeroregionbitparam']))
            print('ReportIntervalSeconds: %d' %(param['report_interval_seconds']))
            print('Minimum Memory Available in GiB Required: %d' %(param['min_memGiB_required']))
            print('Verbose: %d' %(param['verbose']))
            print('')

        print('The initialization parameters for running on a large memory machine are:')
        for param in init_param['large']:
            print('MatrixExponent: %d' %(param['matrix_exponent']))
            print('OpExponent: %d' %(param['op_exponent']))
            print('MaxLevel: %d' %(param['max_level']))
            print('regionbitparam: %d' %(param['regionbitparam']))
            print('zeroregionbitparam: %d' %(param['zeroregionbitparam']))
            print('ReportIntervalSeconds: %d' %(param['report_interval_seconds']))
            print('Minimum Memory Available in GiB Required: %d' %(param['min_memGiB_required']))
            print('Verbose: %d' %(param['verbose']))
            print('')

