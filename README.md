##                   README.md
###################################################################
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

Welcome to the MyPyLARC Python and C Package which is build
on top of the LARC library 
   LARC = Linear Algebra via Recursive Compression.
LARC is useful for carrying out linear algebra on very
large power of 2 dimensioned matrices and vectors which
have some internal repetition or structure.  LARC recursively
compresses matrices and carries out matrix operations
while in the compressed format.

The core code of LARC is written in C and there is a
a SWIG-generated Python interface for ease of use.
MyPyLARC is primarily written in Python but also has some
C code in MyPyLARC/src.

MyPyLARC illustrates how a person could build their own
package on top of the LARC library. 

For more details see:
* slide deck MyPyLARC/larc/doc/aboutLARC.pdf
* paper MyPyLARC/doc/LARCandMyPyLARC.pdf
* tutorials in MyPyLARC/Tutorial 
* examples in MyPyLARC/FFT_play   Eigen_play  Gate_play ...
* and project templates you could copy and modify
  Tutorial/user_proj_template.py and an example of how
  this is used in Tutorial/user_proj_max-A-AT.py
* explanations in MyPyLARC/Eigen_play/how.we.did.this or in
  MyPyLARC/Tutorial/newuser_instructions
* doxygen documentation for larc in MyPyLARC/larc/html

Many of the above are far from finished;
MyPyLARC is a work in progress.

If you want to contact the developers at CCS you can send
email to larc@super.org.

