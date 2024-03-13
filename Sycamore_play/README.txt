/******************************************************************
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC (Linear Algebra via Recursive Compression)                *
 * Authors:                                                       *
 *   - Steve Cuccaro (IDA-CCS)                                    *
 *   - John Daly (LPS)                                            *
 *   - John Gilbert (UCSB, IDA adjunct)                           *
 *   - Mark Pleszkoch (IDA-CCS)                                   *
 *   - Jenny Zito (IDA-CCS)                                       *
 *                                                                *
 * Additional contributors are listed in "LARCcontributors".      *
 *                                                                *
 * Questions: larc@super.org                                      *
 *                                                                *
 * All rights reserved.                                           *
 *                                                                *
 * Redistribution and use in source and binary forms, with or     *
 * without modification, are permitted provided that the          *
 * following conditions are met:                                  *
 *   - Redistribution of source code must retain the above        *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer.                                      *
 *   - Redistribution in binary form must reproduce the above     *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer in the documentation and/or other     *
 *     materials provided with the distribution.                  *
 *   - Neither the name of the copyright holder nor the names of  *
 *     its contributors may be used to endorse or promote         *
 *     products derived from this software without specific prior *
 *     written permission.                                        *
 *                                                                *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       *
 * DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   *
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       *
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   *
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, *
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
 *                                                                *
 *****************************************************************/

Google carried out an experiment using the Sycamore quantum hardware.
The experiment created a quantum state on 53 qubits with the intention
of demonstrating quantum supremacy (a problem that can be done faster
on a quantum computer than on a classical computer).  Originally they
claimed that they could do their computation in 200 seconds compared
to an estimated 10,000 years for a "state-of-the art classical
supercomputer".  However, researchers at IBM Watson responded almost
immediately with a description of how this computation could be
simulated classically in a few days using large amounts of secondary
storage on the Summit supercomputer at Oak Ridge National Labs.

The Sycamore circuit is built on a grid of 53 qubits.  Each
qubit can interact with neighbors in each of the cardinal directions.
There are only four different gates that make up the circuit:
three 1-qubit gates sqrt_X, sqrt_Y, and sqrt_W, where W = (X+Y)/sqrt(2);
and a 2-qubit gate to provide entanglement (see Google paper for
details).

Each cycle in the circuit starts by applying a random 1-qubit gate to
every qubit in the grid pattern. Then there is a step in which a set
of 2-qubit gates are applied simultaneously to pairs of adjacent
qubits in the grid.  There are four different sets A, B, C, and D of
links in the grid which give the locations of the two qubit gates.
The lines of links that are parallel to each other in one direction
contain only links labelled A and B.  And those lines of links
that are in the perpendicular direction contain only links labelled
C and D.   The A's and B's alternate along a given line, and are
in a checkerboard alignment with the lines below and above them.
Similarly the C's and D's.   See the figure below.

We show the set up of a portion of the array of qubits in the Sycamore
circuit tipped at 45 degrees from the figures in the Google reference
to allow us to show the lines of links as vertical and horizontal
in our drawing below.


      /                 qubit   A    qubit   B    qubit  A
     /
top edge                 C             D           C
  /
 /        qubit   A    qubit   B    qubit  A     qubit   B

             C            D           C           D

 qubit A   qubit B      qubit   A    qubit   B    qubit  A

              D           C           D              C
...
           qubit   A    qubit   B    qubit  A     qubit B


In the first step, each qubit is paired
with its neighbor in set A. In the second step, each qubit is paired
with its neighbor in set B. The eight step pattern of cardinal
directions is: A B C D C D A B.  This pattern was selected by Google
because it makes the circuit hard to simulate classically.

To explore LARC's abilities on a sample quantum simulation, this
directory contains a program demonstrating how one might code a small
portion of the Sycamore circuit using LARC.  

References:
* Quantum supremacy using a programmable superconducting processor
  Frank Arute, Kunal Arya, […] John M. Martinis, (of Google)
  Nature 574, 505-510 (Oct 2019).
* Leveraging Secondary Storage to Simulate Deep 54-qubit Sycamore
  Circuits
  Edwin Pednault, John A. Gunnels, Giacomo Nannicini, Lior Horesh,
  Robert Wisnieff, (of IBM Watson)
  Quantum Physics  > arXiv:1910.09534
  Oct 2019.


CONTENTS OF SYCAMORE_PLAY DIRECTORY
===================================

README.txt - This file.

sycamore_generate.py - Python 3 program to generate a Sycamore circuit.
                       Currently, this uses a small 3 by 4 qubit layout.
                       Single qubit gates are selected at random, and
                       then written out to "sycamore_config.json".

sycamore_config.json - JSON file containing a Sycamore circuit. This file
                       contains the information on which single qubit
                       operations to apply during each cycle, so that
                       the exact same Sycamore circuit computation can
                       be replicated.

sycamore_run.py - Python 3 program to build the operator matrices corresponding
                  to the Sycamore circuit in "sycamore_config.json"



BASIC INSTRUCTIONS
=============
How to run a Sycamore circuit simulation.
Inside the top level MyPyLARC directory, compile in complex type
   make TYPE=COMPLEX

Then move into the Sycamore_play subdirectory.
You can use the existing configuration file if you like:
         sycamore_config.json 

The configuration file sycamore_config.json contains
the details of the Sycamore circuit you would like to simulate, e.g.
    - number and arrangement of qubits,
    - randomly generated X,Y,W 1-qubit gates for each cycle
    - which simulation pattern to use for the 2-qubit S gate
       (either [ABCD] "quantum supremacy" gate pattern
      or the [EFGH] pattern, and
    - number of simulation cycles to execute.

If you want a different circuit configuration than what is
currently in sycamore_config.json, then make changes in in
the program sycamore_generate.py  and then run
         python sycamore_generate.py.
	 
To run the simulation (which reads sycamore_config.json) type
       python sycamore_run.py

Output from the simulation during various steps are 
saved in files in a subdirectory matrixFiles.

The final density matrix (of size 2^n by 2^n, where n
is the number of qubits) will be output in two formats.
A file in LARC recursive format
     matrixFiles/sycamore_full_circuit_matrix.json
A file expressed in standard (row major) format 
      matrixFiles/sycamore_row_by_row_dump.txt
      
