 /*****************************************************************
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


In this directory we are using LARC to count the triangles in
graphs and to test our ablility to read sparse matrix market
files we grab from the internet.

Matrix Market Format
====================
Matrix Market format for sparse matrices looks like:

%%MatrixMarket matrix coordinate integer general
% optional comments go here
% and may continue on subsequent lines
4 4 10
1 1 100
1 2 200
1 3 300
1 4 400
2 1 500
3 1 600
4 1 700
2 2 800
3 3 900
4 4 1000

The first line is a header that must have the form
	%%MatrixMarket <object> <format> <field> <symmetry>
LARC currently requires that <object> be "matrix", <format> be "coordinate"
(which together imply a sparse matrix), and <symmetry> be "general".
<field> can have the values "real", "double", "complex", "integer",
or "pattern" (the latter gives only the locations of nonzeros and not their
values, and is not currently supported by LARC).

The header line may be followed by any number (include zero) of comment lines
starting with "%". The first non-comment line consists of three integers
	numberOfRows, numberOfCols, numberOfNonzeroEntries.
which give the size of the matrix and the number of lines which follow.

These following lines have one matrix entry per line with the
entries (row i, col j, value). For complex type, the value field
consists of two real numbers; for example the line
	3 2 14 .3 
would be used for the complex number 14 + 0.3*J.

The data we are using for matrix market files was downloaded
in a different format TSV (with one pair per line - defining an edge)
from the MIT graph collection at:
  https://graphchallenge.mit.edu/data-sets

We then add two lines to the top of the file (our matrix
market reader accepts an edge per line with the end vertices
separated by a space or a tab.

Eventually we want to be able to read their larger files
of 50Meg or more while leaving them in compressed format.

---------------------------------------------------------------
How does LARC count triangles
=============================
LARC counts the number of triangles in a graph as follows.
Start with A = adjacency matrix which as a 1 in entry i,j
if there is an edge between vertex i and j.

The nth power of the adjacency matrix has as its i,j entry
the number of paths with n edges starting at vertex i and
ending at vertex j.

We are interested in the number of triangles.
For each triangle a,b,c, there are 6 paths
a,b,c,a; a,c,b,a; b,c,a,b; etc.
Each of these paths start and end at the same vertex, so
they are counted in an element along the diagonal of the
third power of the adjacency matrix A.

Thus LARC can take 1/6 of the Trace(A^3) to count the
number of triangles.


---------------------------------------------------------------
References:

1. "Design, Generation, and Validation of Extreme Scale Power-Law Graphs"
    by Kepner et al. 	https://arxiv.org/pdf/1803.01281.pdf

2. Paper Kronecker Product Graphs for modeling networks
   http://jmlr.csail.mit.edu/papers/volume11/leskovec10a/leskovec10a.pdf


