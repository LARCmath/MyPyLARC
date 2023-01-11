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
 #*##############################################################*#


from __future__ import print_function, division

import random
import json

## \file sycamore_generate.py
#
#  \brief Creates a Sycamore configuration file which specifies
#  the information needed to perform a single Sycamore run.
#


class Links:
    def __init__(self, linkchars, verbose):
        self.verbose = verbose
        self.linkinfo = {}
        for c in linkchars:
            self.linkinfo[c] = []

    def add_link(self, c, q1, q2):
        self.linkinfo[c].append((q1, q2))
        if self.verbose > 0:
            print("  Link {0:2d} - {1:2d} has type {2}.".format(q1, q2, c))


if __name__ == '__main__':

    verbose = 1
    debug = 0

    #*################################*#
    #*    Basic Parameter Setting     *#
    #*################################*#

    num_rows = 3
    num_cols = 2
    layout = 0
    num_cycles = 13
    output_filename = "sycamore_config.json"

    #*################################*#
    #*    JSON Data Initialization    *#
    #*################################*#

    jdata = {}
    jdata["num_rows"] = num_rows
    jdata["num_cols"] = num_cols
    num_qubits = num_rows * num_cols
    jdata["num_qubits"] = num_qubits
    jdata["layout"] = "SUPREMACY REGIME" if layout==0 else "CLASSICALLY VERIFIABLE"
    jdata["num_cycles"] = num_cycles

    #*################################*#
    #*    Qubit Layout Computation    *#
    #*################################*#

    if verbose > 0:
        print("The following is the {0} by {1} Sycamore layout:".format(num_rows, num_cols))
        print()
        for i in range(num_rows):
           if (i % 2) == 1:
               print("    ", end="")
           for j in range(num_cols):
               qubit_index = (i * num_cols) + j
               print("      {0:2d}".format(qubit_index), end="")
           print()
           print()

    # Determine the connections for each link type.
    if verbose > 0:
        print("Determining connections for each link type:")

    # Do the supremacy configuration ABCD.
    if layout == 0:
        qlinks = Links("ABCD", verbose)
        qpattern = "ABCDCDAB"
        for i in range(num_rows - 1):
            for j in range(num_cols):
                qubit_index = (i * num_cols) + j
                if (i % 2) == 0:
                    if j > 0:
                        qlinks.add_link("C", qubit_index, qubit_index + num_cols - 1)
                    qlinks.add_link("A", qubit_index, qubit_index + num_cols)
                else:
                    qlinks.add_link("D", qubit_index, qubit_index + num_cols)
                    if j < (num_cols - 1):
                        qlinks.add_link("B", qubit_index, qubit_index + num_cols + 1)

    # Do the classical configuration EFGH.
    if layout == 1:
        qlinks = Links("EFGH", verbose)
        qpattern = "EFGHEFGH"
        for i in range(num_rows - 1):
            for j in range(num_cols):
                qubit_index = (i * num_cols) + j
                if (i % 2) == 0:
                    if j > 0:
                        qlinks.add_link("H" if ((i+2*j)%4)==0 else "G",
                                        qubit_index, qubit_index + num_cols - 1)
                    qlinks.add_link("E" if ((i+2*j)%4)==0 else "F",
                                    qubit_index, qubit_index + num_cols)
                else:
                    qlinks.add_link("H" if ((i+2*j)%4)==1 else "G",
                                    qubit_index, qubit_index + num_cols)
                    if j < (num_cols - 1):
                        qlinks.add_link("F" if ((i+2*j)%4)==1 else "E",
                                        qubit_index, qubit_index + num_cols + 1)

    jdata["linkinfo"] = qlinks.linkinfo

    #*################################*#
    #*    Circuit Specification       *#
    #*################################*#

    jdata["circuit"] = []
    last_1qubit_cycle = [random.randint(0,2) for _1 in range(num_qubits)]

    for cycle_index in range(num_cycles):
        next_1qubit_cycle = [(i + random.randint(1,2)) % 3 for i in last_1qubit_cycle]
        jcycle = [None, None]
        jcycle[0] = "".join("XYW"[i] for i in next_1qubit_cycle)
        jcycle[1] = qpattern[cycle_index % 8]
        jdata["circuit"].append(jcycle)
        last_1qubit_cycle = next_1qubit_cycle

    #*################################*#
    #*    Output JSON File            *#
    #*################################*#

    with open(output_filename, "w") as outfp:
        json.dump(jdata, outfp, indent=4)


