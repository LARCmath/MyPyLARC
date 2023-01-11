#!/usr/bin/env python

import numpy as np
import sys

# The following code does not use LARC, it uses numpy to generate
# some a random matrix with know eigenvalues and vectors so that
# we can run other tests using LARC.


# the following code is from StackOverflow, and returns a random orthonormal
# matrix of dimension (dim,dim)
def rvs(dim=3):
        # random_state = np.random
        H = np.eye(dim)                 # identity of size dim,dim
        D = np.ones((dim,))             # all-ones vector of size dim (?)
        print("D is:")
        print(D)
        print("")
        for n in range(1, dim):
                # x = random_state.normal(size=(dim-n+1,))
                x = np.random.normal(size=(dim-n+1,))
                print("x is:")
                print(x)
                print("")
                D[n-1] = np.sign(x[0])
                print("D is:")
                print(D)
                print("")
                # np.sqrt(x*x).sum() is sqrt of sum of squares of elements of x
                x[0] -= D[n-1]*np.sqrt((x*x).sum())
                # this causes x[0] to change sign
                print("If neg add norm x to first term, else subtract:")
                print(x)
                print("")
                # Householder transformation
                # take identity subtract 2 times normalized outer product of x
                Hx = (np.eye(dim-n+1) - 2.*np.outer(x, x)/(x*x).sum())
                print("Hx is:")
                print(Hx)
                print("")
                mat = np.eye(dim)
                # replace a subblock in bottom right of Identity with Hx
                mat[n-1:, n-1:,] = Hx
                print("mat is:")
                print(mat)
                print("")
                # matrix multiply H times mat
                H = np.dot(H, mat)
                print("replace H with H times mat:")
                print(H)
                print("")
        # fix the last sign such that the determinant is 1
        # dim%2 is dim mod 2 and D.prod() is the product of all elements in D
        # D[-1] is just reverse indexing and is the last term of D
        D[-1] = (-1)**(1-(dim%2))*D.prod()
        print("After loop D is:")
        print(D)
        print("")
        # Equivalent to np.dot(np.diag(D), H) but faster, apparently
        # H.T is the transpose of H,
        # v*B is a matrix with row i of B replaced by v_i times that row.
        H = (D*H.T).T
        input("hit return to continue")
        return H


# the following code generates a random real-valued square matrix of dimension
# input by the user (default=8) with known real eigenvalues and eigenvectors.
#
# The matrix E has right eigenvectors that are the rows of the orthonormal
# matrix M, and M[i] has eigenvalue v[i].
if __name__ == '__main__':

        if 2 == len(sys.argv):
                dim = int(sys.argv[1])
        else:
                print("no input dimension given: using default of 8")
                dim = 8;

        print("generating random orthonormal matrix M of dimension %dx%d"
               %(dim,dim))
        M = rvs(dim)
        print(M)
        print("\nthe last row of this matrix m_last is:")
        m_last = M[dim-1]
        print(m_last)
        print("\nconfirming orthonormality of M via M @ M.T: ")
        # @ sign is matrix multiply
        # note: not M*(M.T) - that is element-wise multiplication
        shouldBeye = M @ M.T
        print(shouldBeye)
        print("\nconfirming orthonormality of M via M.T @ M: ")
        alsoShouldBeye = M.T @ M
        print(alsoShouldBeye)
        print("\ngenerating random row vector v of dimension 1x%d" %dim)
        v = np.random.normal(size=dim)
        v = np.sort(v)
        print(v)
        print("\ngenerating diagonal matrix D with these values")
        # replaces each element of diagonal with appropriate element of v
        D = v*np.eye(dim)
        print(D)
        print("\ngenerating matrix E with known eigenvalues, eigenvectors")
        E = (M.T) @ D @ M
        print(E)
        print("\nconfirming m_last is eigenvector with largest eigenvalue")
        w = m_last @ E
        print("\tw = m_last @ E, eigenvalue should be %g" %v[dim-1])
        val = w / m_last        # element-wise division
        for n in range(dim):
                print("n=%d: w[n]/m_last[n] = %g" %(n,val[n]))
