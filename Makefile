##                   Makefile
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
# standard directory paths
BINDIR = bin
OBJDIR = obj
SRCDIR = src
LIBDIR = lib
TESTDIR = tests

# Get the local paths 
include local/Makefile.conf
#    GMPIDIR (gmp-6 include),
#    GMPLDIR (gmp-6 lib),
#    MPIDIR (Anaconda3 include),
#    MPLDIR (Anaconda3 lib)
#    GITLARC (path to for cloning git repository for LARC)

## NOT INCLUDED YET for DOXYGEN
# DOCDIR_HTML = html

# directory paths to LARC Linear Algebra via Recursive Compression package
LARCLIBDIR = larc/lib
ABSLARCLIBDIR = `pwd`/larc/lib
LARCDIR = larc/src

# compiler options
CC = gcc
OPTS = -O2 -g -pg -Wall -fPIC
CFLAGS = $(OPTS) -I$(SRCDIR) -I$(MPIDIR) -I$(GMPIDIR) -I$(LARCDIR) -std=gnu99

## Libraries
LIBS = -L$(MPLDIR) -L$(GMPLDIR) -L$(LARCLIBDIR) -lcurses -lm -lrt -lmpc -lmpfr -lgmp -lpthread -llarc -ltinfo
## SWIG NOTE 1:
# SWIG generates a python package from the C code called swiglarc
# OURLIBS specifies use of the non-shared library liblarc.a, rather than
# the default shared library liblarc.so. This is used to make sure the python
# routines from C code in both LARC and MyPyLARC are accessible to python
# routines.
# Both should be included in the _MyPyLARC.so library.
# The reason it's important to use specify the .a file from LARC instead
# of the .so file is so that LD_LIBRARY_PATH does not need to be set
# to include the directory where liblarc.so resides.
# If we had used the shared .so file then we would have also had to
# include a path to the directory that it is in.
OURLIBS = -L$(MPLDIR) -L$(GMPLDIR) -L$(LARCLIBDIR) -lm -lrt -lmpc -lmpfr -lgmp -lpthread -l:liblarc.a -lcurses -ltinfo

## SWIG NOTE 2:
# The file MyPyLARC_wrap.o is created by SWIG, and therefore may or may not
# be present in the SRCDIR. We use SRCSNOWRAP to make sure we're not expecting
# the *wrap.o file to be in the OBJDIR or in OBJS.
SRCS = $(wildcard $(SRCDIR)/*.c)
SRCSNOWRAP = $(filter-out %_wrap.c, $(SRCS))
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCSNOWRAP))
INCS = $(wildcard $(SRCDIR)/*.h $(LARCDIR)/*.h)

# Python
PYINC = -I$(PYINCDIR)
PYINCDIR := $(shell python3 -c 'import sys; print(sys.prefix+"/include/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+sys.abiflags)')

## Testing, Unit Tests, and Examples
# as a simple test, take the code larc/obj/exampleLARC.o
# and create bin
BIN_LIST = larcExample
BINS = $(patsubst %,$(BINDIR)/%,$(BIN_LIST))

## Library to implement both larc and MyPyLARC code
# this library will be used by python routines
TARGET = lib/libMyPyLARC.a
SO_TARGET = $(patsubst %.a, %.so, $(TARGET))

## makefile magic
.PHONY: all build clean larc #tests install

## 
all: larc $(SRCS) $(BINS) $(SO_TARGET) $(SRCDIR)/_mplSWIG.so # tests $(TARGET)

## TYPE was calculated in the larc/Makefile
export TYPE

## git clone a new copy of larc if needed
larc: 
	if [ -d larc ]; then cd larc && git pull ssh://$(GITLARC) master; else git clone ssh://$(GITLARC) ./larc; fi
	cd larc && $(MAKE)

## Create the MyPyLARC libraries
$(TARGET): $(OBJS)
	@test -d $(LIBDIR)/ || mkdir -p $(LIBDIR)
	ar rcs $@ $(OBJS)
	ranlib $@

$(SO_TARGET): $(TARGET) $(OBJS)
	$(CC) -shared -o $@ $(OBJS)

## silently create the following directories if need be
build: 
	@mkdir -p obj
	@mkdir -p bin
	@mkdir -p lib

## make clean commands
clean: 
	rm -rf $(OBJS) $(TESTS) $(BINS) 
	rm -f $(TARGET) $(SO_TARGET)
	rm -f $(SRCDIR)/*.pyc $(SRCDIR)/_mplSWIG.so $(SRCDIR)/*_wrap.o $(SRCDIR)/*_wrap.c $(SRCDIR)/mplSWIG.py
	rm -f tests/tests.log
	find . -name "*.gc*" -exec rm {} \;
	rm -rf `find . -name "*.dSYM" -print`
	rm -rf src/__pycache__
	if [ -d larc ]; then cd larc && make clean; fi

## If we add the DOXYGEN stuff above we need to add them to the clean
#	rm -rf $(DOCDIR_HTML)


# we keep all the SWIG dependent stuff in SRCDIR to avoid problems later
# OURLIBS links liblarc.a directly into _mplSWIG.so, so we don't need to set
# LD_LIBRARY_PATH to find the shared liblarc.so
# since we are directly linking in the .a file, it must be in the dependencies
$(SRCDIR)/_mplSWIG.so: $(SRCDIR)/mplSWIG_wrap.o $(OBJS) $(LARCLIBDIR)/liblarc.a
	$(CC) $(CFLAGS) -shared $^ -o $@ $(OURLIBS)

$(SRCDIR)/mplSWIG_wrap.o: $(SRCDIR)/mplSWIG_wrap.c
	$(CC) $(CFLAGS) $(PYINC) -c $< -o $@


# using SRCSNOWRAP avoids a dependency loop (mplSWIG_wrap.c shouldn't depend
# on itself)
$(SRCDIR)/mplSWIG_wrap.c: $(SRCDIR)/mplSWIG.i $(SRCSNOWRAP) $(INCS)
	swig -python $<

# The binary executables do not depend on the python wrapper.
# Use of the ABSLARCLIBDIR path allows us to call the executables from anywhere
# (not just the directory the Makefile resides in)
$(BINS): $(OBJS) $(BINDIR)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -Wl,-rpath=$(ABSLARCLIBDIR)

$(BINDIR): 
	@test -d $(BINDIR)/ || mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCS)
	@test -d $(OBJDIR)/ || mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c -o $@ $< 


