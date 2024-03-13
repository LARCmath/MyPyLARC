##                   README.md
  Copyright (C) 2014-2024, Institute for Defense Analyses        
  4850 Mark Center Drive, Alexandria, VA; 703-845-2500           
  This material may be reproduced by or for the US Government    
  pursuant to the copyright license under the clauses at DFARS   
  252.227-7013 and 252.227-7014.                                 
                                                                 
  LARC : Linear Algebra via Recursive Compression                
  Authors:                                                       
    - Steve Cuccaro (IDA-CCS)                                    
    - John Daly (LPS)                                            
    - John Gilbert (UCSB, IDA adjunct)                           
    - Mark Pleszkoch (IDA-CCS)                                     
    - Jenny Zito (IDA-CCS)                                       
                                                                 
  Additional contributors are listed in "LARCcontributors".      
                                                                 
  Questions: larc@super.org                                      
                                                                 
  All rights reserved.                                           
                                                                 
  Redistribution and use in source and binary forms, with or     
  without modification, are permitted provided that the          
  following conditions are met:                                  
    - Redistribution of source code must retain the above        
      copyright notice, this list of conditions and the          
      following disclaimer.                                      
    - Redistribution in binary form must reproduce the above     
      copyright notice, this list of conditions and the          
      following disclaimer in the documentation and/or other     
      materials provided with the distribution.                  
    - Neither the name of the copyright holder nor the names of  
      its contributors may be used to endorse or promote         
      products derived from this software without specific prior 
      written permission.                                        
                                                                 
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         
  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    
  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   
  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             
                                                                 

### INTRO:

Welcome to the MyPyLARC Python and C Package which is built
on top of the LARC library 

    LARC = Linear Algebra via Recursive Compression.

LARC is useful for carrying out linear algebra on very
large power-of-2-dimensioned matrices and vectors which
have some internal repetition or structure.  LARC recursively
compresses matrices and carries out matrix operations
while in the compressed format.
MyPyLARC has tutorials and sample applications to help
you learn about LARC and it illustrates how a developer
could build their own package on top of the LARC library.

LARC version 2.1 has not been optimized for any particular
application.  MyPyLARC version 2.1 is still an evolving
package.


### LANGUAGES USED:

MyPyLARC is primarily
written in Python but also has some C code in MyPyLARC/src.
There are also several Jupyter Notebook demonstrations.
The core code of LARC is written in C and there is a
a SWIG-generated Python3 interface for ease of use.
The python module pylarc.py imports SWIG-generated python
wrappers for the C code and contains additional
routines that are written directly in Python.


### DOCUMENTATION AND EXAMPLES:

For more details see:
\* slide deck MyPyLARC/larc/doc/aboutLARC.pdf
\* paper MyPyLARC/doc/LARCandMyPyLARC.pdf
\* poster MyPyLARC/doc/LARCposterSIAMAnnualMeeting2020.pdf
\* tutorials in MyPyLARC/Tutorial 
\* examples in MyPyLARC/FFT\_play   Eigen\_play  Gate\_play ...
\* and project templates you could copy and modify
  Tutorial/user\_proj\_template.py and an example of how
  this is used in Tutorial/user\_proj\_max-A-AT.py
\* explanations in MyPyLARC/Eigen\_play/how.we.did.this or in
  MyPyLARC/Tutorial/newuser\_instructions
\* doxygen documentation for larc in MyPyLARC/larc/html
\* sample Jupyter Notebook code running LARC in
  MyPyLARC/Eigen\_play/physics\_work/EigenRun.ipynb.
  
To view the doxygen documentation,

        gio open MyPyLARC/html/index.html

### DEVELOPERS:

If you want to contact the developers at CCS you can send
email to larc@super.org.  Much of the code is far from finished,
so we would appreciate your comments and questions.


### OPERATING SYSTEMS:

We have tested LARC and MyPyLARC on Linux-Centos7 and
Ubuntu-18.04 systems and include our best guess below
about how you would make LARC and MyPyLARC work on
other systems.


### NECESSARY SOFTWARE:

As of the August 2022 version of MyPyLARC and LARC,
you will need to have software
packages in your environment for the following:
\*   Python 3.6 or later
\*   GNU GMP 5.0 or later
\*   GNU MPFR 4.0.0 or later
\*   GNU MPC 1.1.0 or later
\*   SWIG 3.0.0 or later (to create Python interface for C routines)
\*   Doxygen-1.8.13 exact version (for documentation)
\*	* 1.8.x works but with warnings that can be ignored

We also recommend
\*   GCC 10 or later

Python3 and the multiprecision libraries can be obtained from the
Anaconda package, version 2020.02 or later, or may be downloaded and
installed individually.  
These should be available on line for upload at several
places including: https://anaconda.com, https://www.swig.org,
https://www.doxygen.nl, and https://gcc.gnu.org.
In the instructions below we will tell you how to update
your MyPyLARC/local/Makefile.conf and your
MyPyLARC/larc/local/Makefile.conf to reflect the locations
of this software in your environment.

### GETTING STARTED:

If you are reading this, you probably have already cloned the git repository
GitHub/LARCmath/MyPyLARC. If so, you may skip step I below.

#### I. Clone MyPyLARC
\*    cd [directory location you want to put repository]
\*    git clone https://github.com/LARCmath/MyPyLARC  MyPyLARC
-      NOTE: You could designate another name than MyPyLARC:
\*    git clone https://github.com/LARCmath/MyPyLARC  MyProjectName
-    We will use the name MyPyLARC in the discussion below.

There is also separate repository for GitHub/LARCmath/LARC
which MyPyLARC will update automatically and keep in a
subdirectory called larc. Before the first time you "make"
it is probably easier to clone this yourself.

#### II. Clone LARC
\*    cd MyPyLARC
\*    git clone https://github.com/LARCmath/LARC larc
-     NOTE: on this last step we are specifying that the
            subdirectory containing the LARC repository is
	    named "larc".  This is necessary since each time
	    you "make" MyPyLARC will check to see if you
            have the latest version of LARC and if not update
	    the repository in the larc subdirectory.

#### III. Modify the configuration files to reflect your local environment.

Before compiling, it is necessary to create two local configuration files
in MyPyLARC/local/Makefile.conf and MyPyLARC/larc/local/Makefile.conf. These
files define the paths to the libraries required by LARC and MyPYLARC. We
have provided some example files in the MyPyLARC/local subdirectory, but it
is likely that you will need to modify these examples for your system, 
depending on the install locations for your libraries. 

MyPyLARC/larc/doc/build\_notes.pdf contains some basic information on the
libraries that are required, and also their locations in a "standard" Linux
installation.

\* cp larc/local/Makefile.conf.ubuntu1804 larc/local/Makefile.conf
\* cp local/Makefile.conf.ubuntu1804 local/Makefile.conf
\* [edit these files as needed]


#### IV. Compiling

The following Makefile commands tell LARC to define its
scalarType to be one of the pre-defined C datatypes:
\*    make TYPE=COMPLEX (long double complex)
\*    make TYPE=REAL (long double)
\*    make TYPE=INTEGER (int64\_t)

as well as five multiprecision Scalartypes
\*    make TYPE=MPINTEGER    (mpz\_t)
\*    make TYPE=MPREAL       (256-bit mpfr\_t)
\*    make TYPE=MPCOMPLEX    (256-256-bit mpc\_t)
\*    make TYPE=MPRATIONAL   (mpq\_t)
\*    make TYPE=MPRATCOMPLEX (larc\_mpratcomplex\_t)

or one of the specialty types:
\*    make TYPE=BOOLEAN      (int64\_t)
\*    make TYPE=CLIFFORD     (clifford\_t)
\*    make TYPE=UPPER        (larc\_exponent\_scalar\_t)
\*    make TYPE=LOWER        (larc\_exponent\_scalar\_t)

the MPRATCOMPLEX type is constructed from two copies of mpq\_t
for the real and imaginary parts of a complex rational number.
If no TYPE option is specified, the default is REAL (C long double).

Due to our setup, the larc directory contains a current copy of the LARC
software. If this software is later updated, the Makefile will "git pull" the
newer version of LARC automatically. If this behavior is not desired, simply
comment out this line of the MyPyLARC makefile by prefixing a "#" character:

        if [ -d larc ]; then cd larc && git pull https://$(GITLARC) master; else git clone https://$(GITLARC) ./larc; fi

#### V. Tutorial

You should now be able to run all the programs in the MyPyLARC tutorial and 
example package. We recommend that you first look in the Tutorial subdirectory
for various tutorials which explain in basic terms how to use LARC. The
Tutorial/newuser\_instructions file contains suggestions for which tutorials
you should try first. There are also several applications of LARC which are
found in subdirectories of the MyPyLARC main directory.


#### VI. Index of Contents of MyPyLARC directory

bin                exampleMPL, a binary executable

Count\_triangles   graph adjacency matrix python programs and data subdirectory

doc                LARCandMyPyLARC paper and SIAM poster

Eigen\_play        Routines for solving for eigenvalues including
                   power\_method.py MyPyLARC implementation of power method;
                   README  documentation of proposed hybrid-Lanczos method;
                   Data subdirectory;  Experiments subdirectory;
		   The physics\_work subdirectory (including Jupyter notebook
		   code) includes a project for finding minimum energy
		   eigenstates.
		   
FFT\_play          Routines for experiments with large Fast Fourier Transforms

Gate\_play         Creation reversible circuits using the universal gate library

Hash\_play         Details of locality-sensitive hash and ideas on
                   possible implementation of locality-preserving hash
		   
html               Files supporting doxygen documentation.

InitParams         Parameter files for initialization used by code in other
                   project directories such as Tutorial, FFT_play, and code
		   for user to make their own initialization file.
		   
larc               Sub-clone of the LARC repository with workhorse LARC routines

lib                File to ensure the creation of a lib directory

local              Directory holding configuration files for compiling

Makefile           File governing the compilation of MyPyLARC

obj                Directory containing a few .o files

README.md          This file of documentation

Roots\_play        Directory holding a test program using the roots of unity.

src                Source directory for .c and python SWIG wrappers

Sycamore\_play     Directory contain experiments for simulating a portion
                   of the Sycamore quantum circuit.

tests              Some small nearly empty test directories.

Tutorial           Directory of Tutorials and Examples (see README.md in that
                   directory for list of routines)


#### VII. List of C code math operations from larc/src/matmath.h
          (For details, also see the doxygen page on MyPyLARC/larc/matmath.h)

void get\_array\_of\_scalars\_in\_larcMatrixFile(scalarType \*\*scalars\_ptr,
                                   int64\_t \*numScalars, char \*path);
				   
int64\_t matrix\_add(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_diff(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_mult(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_mult\_clean(int64\_t A\_pID, int64\_t B\_pID,
                                   mat\_level\_t cleanThresh);
				   
int64\_t scalar\_mult(int64\_t A\_pID, int64\_t B\_pID);

int64\_t scalar\_divide(int64\_t A\_pID, int64\_t B\_pID);

int64\_t kronecker\_product(int64\_t A\_pID, int64\_t B\_pID);

int64\_t join(int64\_t A\_pID, int64\_t B\_pID);

int64\_t stack(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_entrySquared(int64\_t m\_pID, char \*scale\_factor);

int64\_t iHadamard\_times\_matrix(int64\_t A\_pID);

int64\_t matrix\_basischange\_A\_by\_B(int64\_t B\_pID, int64\_t A\_pID);

int64\_t matrix\_saxpy(int64\_t A\_pID, int64\_t B\_pID, 
                      int64\_t scalar\_a\_pID, int64\_t scalar\_b\_pID);
			      
int64\_t vector\_dot\_product(int64\_t A\_pID, int64\_t B\_pID, int verbose);

int64\_t matrix\_times\_iHadamard(int64\_t A\_pID);

char \*tracenorm(int64\_t m\_pID, char \*scale\_factor);

char \*traceID(int64\_t m\_pID);

char \*get\_scalar\_value\_string(int64\_t m\_pID);

char \*matrix\_count\_entries(int64\_t mat\_pID, char \*scalar\_str);

char \*get\_list\_\of\_scalars\_in\_larcMatrixFile(char \*path);

int64\_t random\_bool\_matrix\_from\_count(mat\_level\_t row\_level,
	mat\_level\_t col\_level, char\* numOnes);
					
int64\_t matrix\_element\_with\_maxNorm(int64\_t mat\_pID);

int64\_t normID(int64\_t mat\_pID, int whichNorm);

int64\_t create\_const\_matrix(int64\_t constID, mat\_level\_t row\_level,
	mat\_level\_t col\_level);

int64\_t apply\_function\_to\_matrix\_values(int64\_t m\_pID,
        void (*func)(scalarType*, const scalarType), op\_type\_t op\_memz);

