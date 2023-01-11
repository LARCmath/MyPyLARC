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


