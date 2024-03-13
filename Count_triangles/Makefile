#*                   Makefile
#*#################################################################
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

# standard directory paths
BINDIR = bin
OBJDIR = obj
SRCDIR = src


# directory paths to LARC Linear Algebra via Recursive Compression package
LARCLIBDIR = ../larc/lib
LARCDIR = ../larc/src

# Get the local paths
ifeq ("$(wildcard ../local/Makefile.conf)","")
  $(info  ERROR: You need to create a "local/Makefile.conf" file for both MyPyLARC (./local/Makefile.conf) and LARC (./larc/local/Makefile.conf).)
  $(info  NOTE: This needs to be done in stages. First, create local/Makefile.conf.  Second, run "make", which will fail but will pull LARC from github.  Third, create larc/local/Makefile.conf.  Fourth, run "make clean" and then "make" again.)
  $(info  NOTE: Sample files for various environments can be found in the "local" directory.)
  $(error Terminating makefile)
endif
include ../local/Makefile.conf

# Which binaries to build
BIN_LIST = triangle_counter
BINS = $(patsubst %,$(BINDIR)/%,$(BIN_LIST))

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))
INCS = $(wildcard $(SRCDIR)/*.h $(LARCDIR)/*.h)

# compiler options
CC = gcc
OPTS = -O2 -g -pg -Wall -fPIC
CFLAGS = $(OPTS) -I$(SRCDIR) -I$(MPIDIR) -I$(GMPIDIR) -I$(LARCDIR) -std=gnu99

LIBS = -L$(MPLDIR) -L$(GMPLDIR) -L$(LARCLIBDIR) -lcurses -lm -lrt -lmpc -lmpfr -lgmp -lpthread -llarc -ltinfo

# Rules

#* makefile magic
.PHONY: all clean

all: $(BINS)

$(OBJDIR): 
	@test -d $(OBJDIR)/ || mkdir -p $(OBJDIR)

$(BINDIR): 
	@test -d $(BINDIR)/ || mkdir -p $(BINDIR)


$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCS)
	@test -d $(OBJDIR)/ || mkdir -p $(OBJDIR)
	$(CC) -c -o $@ $(CFLAGS) $<

$(BINS): $(OBJS) $(BINDIR)
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LIBS) -Wl,-rpath=$(MPLDIR):$(GMPLDIR):$(LARCLIBDIR)


#* make clean commands
clean: 
	rm -rf $(OBJS) $(BINS) 

