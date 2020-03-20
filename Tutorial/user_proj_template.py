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

 #Your_Project_template
 # Some example application "playgrounds" for eigen, fft, toeplitz
 # can be found under mypylarc


 #Imports of utility python code
from __future__ import print_function

import os
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../src"))
import MyPyLARC as mypy
import numpy as np
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print("This code provides a template for using basic LARC routines \n")

    print("At this point in the code one specifies basic LARC parameters ")
    print("It would be good to have these as command line parameters ")
    print("and have them echo-ed out at the beginning of each run ")
    print("It also would be good to start with toy values during debugging ")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    rnd_sig_bits = 40   # (default value -1 =>62 bits)
    trunc_to_zero_bits = 20  # (default value -1 =>62 bits)
    verbose = 1

    mypy.create_report_thread(1800)
    mypy.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)

    print("mat_store_exp = ",mat_store_exp)
    print("op_store_exp = ",op_store_exp)
    print("max_level = ",max_level)
    print("rnd_sig_bits = ",rnd_sig_bits,"; ( -1 implies a default 62 bits)")
    print("trunc_to_zero_bits = ",trunc_to_zero_bits,"; ( -1 for  62 bits)")
    print("verbose = ",verbose)

    # scalarType is set at compile time (e.g. make TYPE=INTEGER) 
    # this command returns the scalar type in case you need to know it
    scalarTypeStr = mypy.cvar.scalarTypeStr;

    #Create, Convert and Load auxiliary scalars
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C/Python subroutines:...
	     # ...

	     #SPECIFY IN C OR PYTHON or numpy
	     #const_pi = 3.141592653      #Real
	     #const_root2 = sqrt(2.)      #Real
	     #constant_eipi/4 = [1/const_root2,1/const_root2]  #Complex
	     #CONVERT AND LOAD INTO LARC
	     #...
	     #Or RECURSIVE GENERATION

    #Create, Convert and Load auxiliary row vectors
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C/Python subroutines:...
	     # ...
	     #SPECIFY IN C OR PYTHON or numpy
	     # List1 = [1., 2., 3., 4.,]
	     # V1Row = [1., 2., 3., 4.,]
	     #CONVERT AND LOAD INTO LARC
	     #...
	     # or RECURSIVE GENERATION

    #Create, Convert and Load auxiliary column vectors
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C/Python subroutines:...
	     # ...
	     # List1 = [1., 2., 3., 4.,]
	     # V1Col = [[1.], [2.], [3.], [4.],]
	     #...
	     # or RECURSIVE GENERATION

    #Create, Convert and load auxiliary matrices
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     #...
	     # List1 = [1., 2., 3., 4.,]
	     # M = [ [1. , 2.] , [3. , 4.] ]
	     # ...
	     #CONVERT AND LOAD INTO LARC
	     #...
	     # or RECURSIVE GENERATION


    #Perform LARC Linear Operations (or LARC/PYTHON hybrid operations???)
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     # ...


    #Evaluate Properties of LARC objects (norms, is_symmetric, is_toeplitz, ...
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     # ...

    #Output Results
    	     # See list of LARC subroutines:...
    	     # See list of numpy subroutines:...
    	     # See list of C / Python subroutines:...
	     # ...

