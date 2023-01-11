#!/usr/bin/env python3

 #*################################################################
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
 #*################################################################

from __future__ import print_function

import numpy as np
import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
# sys.path.append(os.path.join(os.path.dirname(__file__),"../larc/src"))
# import pylarc
from ctypes import *
import json

##
# \file tut1_larc_overview.py
#
# \brief This program illustrates LARC initialization
# and some basic LARC operations.
#
# The specific operations shown are
# based on the user's response to various prompts.
# With each operation, this program provides an explanation
# of what is happening under the covers of the LARC package.

if __name__ == '__main__':

    verbose = 1


    #*##################################################
    #*   Print welcome message for LARC and MyPyLARC  ##
    #*##################################################

    # look in How to Write Python
    # give a menu
    #     0 to see read two page overview of LARC and MyPyLARC
    #     1 getting started in LARC, initialization, io and creating matrices
    #     2 recursive operations in LARC, the matrix store, the operations store
    #     3 writing your own recursive function for LARC

    loop = True

    while (loop == True):
        print("")            
        print("")            
        print("**********************************************************************")
        print("*           WELCOME to the MyPyLARC and the LARC Tutorial #3         *")
        print("**********************************************************************")
        print("*                                                                    *")
        print("*  Hashing Tutorial Menu:                                            *")
        print("*                                                                    *")
        print("*    0. exit tutorial,                                               *")
        print("*                                                                    *")
        print("*    1. Introduction - Hashing in LARC                               *")
        print("*              explain_hashing_uses_in_LARC()                        *")
        print("*                                                                    *")
        print("*    2. Multiplicative Fibonacci hashing - the basic building block  *")
        print("*                explain_multiplicative_Fibonacci_hash()             *")
        print("*                                                                    *")
        print("*    3. Hashing lists and accessing the                              *")
        print("*       MatrixStore and OperationStore                               *")
        print("*               explain_hashing_lists()                              *")
        print("*                                                                    *")
        print("*    4. Hashing the different scalarTypes                            *")
        print("*               explain_hashing_scalarTypes()                        *")
        print("*                                                                    *")
        print("*    5. Locality sensitive hashing techniques, snapping              *")
        print("*       scalars, regional representative retrieval,                  *")
        print("*       SPR mode with probabilistic success vs.                      *")
        print("*       MAR mode with a 'guarantee of success'                       *")
        print("*               explain_hashing_for_snapping_scalars()               *")
        print("*                                                                    *")
        print("*    6. How hash filters are used to speed searches                  *")
        print("*               explain_hash_filters()                               *")
        print("*                                                                    *")
        print("*    7. Getting statistics on your hash chains,                      *")
        print("*       theoretical expectations if hash was random.                 *")
        print("*               explain_hash_statistics()                            *")
        print("*                                                                    *")
        print("*    8. Hash and hash table demo                                     *")
        print("*                                                                    *")
        print("*    9. RESOURCES: papers, slides, documentation, tutorial routines  *")
        print("*                                                                    *")
        print("**********************************************************************")
        print("")
        
        user_input=input("Make a Tutorial selection by entering the item number: ")

        print("")
        print("")

        # if (user_input in['0','1','2','3','4'] == False):
        #     print("")
        #     print("")
        #     print("")
        #     print("        !!!!Please enter a number from the Tutorial Menu!!!!")

        if user_input=='0':
            print("\n\tExiting Hashing Tutorial.\n\tPlease visit again!\n")
            sys.exit()

        if user_input=='1':
            print("**********************************************************************")
            mypy. explain_hashing_uses_in_LARC()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")
            
        if user_input=='2':
            print("**********************************************************************")
            mypy.explain_multiplicative_Fibonacci_hash()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")

        if user_input=='3':
            print("**********************************************************************")
            mypy.explain_hashing_lists()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")

        if user_input=='4':
            print("**********************************************************************")
            mypy.explain_hashing_scalarTypes()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")

        if user_input=='5':
            print("**********************************************************************")
            mypy.explain_hashing_for_snapping_scalars()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")

        if user_input=='6':
            print("**********************************************************************")
            mypy.explain_hash_filters()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")

        if user_input=='7':
            print("**********************************************************************")
            mypy.explain_hash_statistics()
            print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")

        if user_input=='9':
            # print("**********************************************************************")
            print("")
            print("")
            mypy.list_explanatory_resources()
            print("")
            print("")
            # print("**********************************************************************")
            Userinput= input("                    Press <Enter> to return to Menu\n")
            
        # if user_input in['3']:
        # if user_input in['0','1','2','3','4']:
        # if user_input in['0','1','2','3']:
        #  break
        #        Userinput= input("Press <Enter> to continue\n")
        
        if user_input=='8':
            print("")
            print("")
            print("           Hash and Hash Table Demo    ")
            print("")
            print("------------------------------------------------")
            print("In this input/output demonstration we initialize the LARC")
            print("matrix store and operations store with some standard parameters,")
            print("read in a few matrices, then give an example hash chain output.")
            print("")
            print("For the purposes of this demo, we are initializing LARC with")
            print("tiny hash tables, so that more than one value may be found in")
            print("a hash chain.")
            print("------------------------------------------------")
            
            verbose = 0
            mat_store_exp = 5
            op_store_exp = 5
            max_level = 10
            regionbitparam = -1
            zeroregionbitparam = -1
            mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)
            if verbose:
                mypy.create_report_thread(report_interval_seconds)
                
            print("")
            print("LARC has been initialized with a 32-long MatrixStore hash table,")
            print("and a 32-long OperationStore hash table.")
            print("")
            scalarTypeStr = mypy.cvar.scalarTypeStr
            print("LARC is currently compiled with scalarType ")
            print(scalarTypeStr)
            print("")
            print("------------------------------------------------")
            print("")
            Userinput= input("Press <Enter> to continue\n")
            print("")
    
            # Calculate number of matrices created, then print part of matrix store
            num_matrices_made = mypy.num_matrices_created()
            #print("\n%d matrices have been created" %num_matrices_made)
            end = num_matrices_made - 1
            filename = "Data/Out/preload.%s.store" %scalarTypeStr
            iS_string = "After preload with parameters: " + str(mat_store_exp) + ", ";
            iS_string = iS_string + str(op_store_exp) + ", " + str(max_level) + "."
            mypy.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),iS_string)
            #print("Here is the matrix store after preloading is complete.")
            #print("")
            #os.system("cat "+filename)
            #print("")
            #Userinput= input("Press <Enter> to continue\n")
            #print("**********************************************************************")
            
            print("We now load some additional matrices into the store.")
            
            filename = "Data/In/sample.1.2.%s.json" %scalarTypeStr
            samp_pID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
            
            #*###########################
            # create a matrix in python #
            #*###########################
            if scalarTypeStr in ('Integer', 'MPInteger'):
                a_str = [1, 3, 5, 6,
                         8, 6, 3, 1,
                         -9, 11, 13, 15,
                     16, 13, 12, 10]
            elif scalarTypeStr == 'Boolean':
                a_str = [1, 0, 0, 0,
                         0, 0, 0, 1,
                         0, 1, 1, 1,
                         1, 1, 1, 0]
            elif scalarTypeStr in ('Complex', 'MPComplex', 'MPRatComplex'):
                a_str = [1+2j, 3+4j, 5+6j, 7+8j,
                         8+7j, 6+5j, 3+4j, 1+2j,
                         9+10j, 11+12j, 13+14j, 15+16j,
                         16+15j, 14+13j, 12+11j, 10+9j]
            elif scalarTypeStr in ('Real', 'MPReal', 'MPRational', 'Clifford'):
                a_str = [1, 3, 5, 6,
                         8, 6, 3, 1,
                         -9, 11, 13, 15,
                         16, 13, 12, 10]
            elif scalarTypeStr in ('Upper', 'Lower'):
                a_str = [0.1, 0.3, 0.5, 0.6,
                         0.8, 0.6, 0.3, 0.1,
                         0, 1, 0.013, 0.15,
                         0.16, 0.13, 0.12, 0.10]
            else:
                raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)
            
            # turn the string into an array
            #a_arr = list(map(str,a_str))
            a_arr = mypy.map_to_str(a_str,scalarTypeStr)
            
            # parameters for entering the python array into the store
            level = 2
            dim_whole = 2**level
            
            pID = mypy.row_major_list_to_store(a_arr, level, level, dim_whole)
            
            # make a parent matrix from four copies of the packedID matrix
            panel = [pID]*4   # alternatively panel=[pID,pID,pID,pID]
            pID_parent = mypy.get_pID_from_four_sub_pIDs(pID,pID,pID,pID,3,3)
            
            filename = "Data/In/nand.%s.json" %scalarTypeStr
            nand_pID = mypy.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
            
            filename_rmm = "Data/In/sample.1.1.%s.rmm" %scalarTypeStr
            
            sample_pID = mypy.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm)) 
            
            # make CNOT
            CNOT_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0]))
            CNOT_pID = mypy.row_major_list_to_store(CNOT_arr,level,level,dim_whole)
            
            # Calculate number of matrices created, then print part of matrix store
            num_matrices_made = mypy.num_matrices_created()
            print("\n%d matrices have been created" %num_matrices_made)
            start = end + 1
            end = num_matrices_made - 1
            #filename = "Data/Out/cnot.%s.store" %scalarTypeStr
            #mypy.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")
            #print("Here are the new matrices added to the store.")
            #print("(Note how some matrices may have already been in the store.")
            #print("")
            #os.system("cat "+filename)
            #print("")
            #Userinput= input("Press <Enter> to continue\n")
            #print("**********************************************************************")
            
            
            # build Zero matrices
            Z2_arr = list(map(str,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
            Z2_pID = mypy.row_major_list_to_store(Z2_arr,level,level,dim_whole)
            
            # build Identity matrices
            I2_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]))
            I2_pID = mypy.row_major_list_to_store(I2_arr,level,level,dim_whole)
            
            # build a doubly-controlled NOT (base unit for reversible computing)
            TOFFOLI_pID= mypy.get_pID_from_four_sub_pIDs(I2_pID,Z2_pID,Z2_pID,CNOT_pID,3,3)
            filename = "Data/Out/toffoli.%s.naive" %scalarTypeStr
            
            # use CCNOT to build an NAND
            #print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
            NOT_packedID = mypy.cvar.packedID_NOT;
            I1_packedID = mypy.get_identity_pID(1);
            not3_packedID = mypy.kronecker_product(I1_packedID,mypy.kronecker_product(I1_packedID, NOT_packedID));
            nand_from_Toff_packedID = mypy.matrix_mult(not3_packedID,TOFFOLI_pID);
            filename = "Data/Out/nandfromtoff.%s.naive" %scalarTypeStr
            mypy.fprint_naive(nand_from_Toff_packedID,os.path.join(os.path.dirname(__file__),filename))
            
            samp2_pID = mypy.scalar_mult(mypy.cvar.packedID_scalarM1,samp_pID)
            
            samp3_pID = mypy.matrix_add(samp_pID,samp2_pID)
            
            samp3_pID = mypy.adjoint(samp_pID)
            
            samp4_pID = mypy.matrix_mult(samp_pID,samp3_pID)
            samp4_pID = mypy.kronecker_product(samp_pID,samp_pID)
            samp4_pID = mypy.join(samp_pID,samp_pID)
            samp4_pID = mypy.stack(samp_pID,samp_pID)
            
            num_matrices_made = mypy.num_matrices_created()
            print("\n%d matrices have been created" %num_matrices_made)
            print("with hash tables of size",1<<mat_store_exp)
            #print("Previous range printed ended with matrixID %d\n" %end)
            #if (end == num_matrices_made-1) :
            #    print("Nothing new since last matrix store print\n")
            #else :
            #    start = end + 1
            #    end = num_matrices_made - 1
            #    filename = "Data/Out/nand.%s.store" %scalarTypeStr
            #    mypy.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__),filename),"Loaded NAND")
            #print("Here is the matrix store after preloading is complete.")
            #print("")
            #os.system("cat "+filename)
            print("")
            print("------------------------------------------------")
            print("")
            Userinput= input("Press <Enter> to continue\n")
            print("")
            
            # get the hashID and print the hash chain corresponding to a matrix
            print("We now find the hash chain containing one of the added matrices.")
            print("The matrix has packedID",nand_pID,".")
            print("Since this is a nonscalar matrix, the hash chain will be taken")
            print("from the nonscalar hash table in the matrixStore, and the")
            print("function hash_pID will perform a recursive hash on the four")
            print("submatrix packedIDs that define this matrix.")
            hashID = mypy.hash_pID(nand_pID)
            print("the hash_pID function returns the value",hashID,".")
            print("")
            print("------------------------------------------------")
            print("")
            Userinput= input("Press <Enter> to continue\n")
            print("")
            comment = "hash chain for nand"
            #filename = "Data/Out/hashChain.beforeMatrixRemove"
            #out_path =  os.path.join(os.path.dirname(__file__),filename)
            mypy.print_matrix_hash_chain_info(hashID, comment)
            print("")
            print("We see that the value",nand_pID,"appears in this hash chain.")
            print("")
            print("")
            print("------------------------------------------------")
            print("")
            Userinput= input("                    Press <Enter> to return to Menu\n")
            #print("**********************************************************************")
            # Userinput= input("Press <Enter> to continue\n")
            # sys.exit(0)
            

