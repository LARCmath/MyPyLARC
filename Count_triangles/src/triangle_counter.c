//                    triangle_counter.c
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC : Linear Algebra via Recursive Compression                *
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


// Standard Libaries
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

// Our header files structures and functions
#include "global.h"
#include "larc.h"
#include "hash.h"
#include "io.h"
#include "json.h"
#include "matmath.h"
#include "matrix_store.h"
#include "op_store.h"
#include "organize.h"
#include "scalars.h"
#include "exampleLARC.h"

/*!
 * \file exampleLARC.c
 * \brief Test program that runs through several LARC capabilities
 */

// Sources of Triangle Data:
//      /g01/larc/TriangleData/MIT_GraphChallenge
//      /g01/larc/TriangleData/Stanford_Data   (SNAP)

// Test file for matrix market format:
// ../larc/tests/dat/in/matrixMarketExchange1 was copied to Data subdirectory




// The following function is contained within a CPP block that makes it
// invisible to SWIG. This is important because its arguments are of type
// scalarType, which may be a C structure which should not be passed to
// Python code.
#ifndef SWIG
/*!
 * \ingroup larc
 * \brief A function that returns the norm squared of a scalar argument
 * \param y A pointer to the return value
 * \param x The input scalarType value
 */
static void square_scalar(scalarType *y, const scalarType x)
{
        sca_conj(&scratchVars.calc_conj,x);
        sca_mult(y,scratchVars.calc_conj,x);
}
#endif

/*!
 * \ingroup LARC
 * \brief A function which uses apply_function_to_matrix_values to square
 * each element of a matrix
 *
 * The function takes as input a matrixID for a matrix A, and generates a new
 * matrix B whose elements B_{ij} = |A_{ij}|^2. The new matrix is added to the
 * matrixStore and its matrixID returned.
 *
 * \param input_ID The matrixID of the matrix whose elements are to be squared
 * \param op_memz One of the operation types FUNC_A, FUNC_B, ... left available for memoizing user-defined functions
 * \result The matrixID of the new matrix
 */
int64_t square_matrix_elements(int64_t input_ID, op_type_t op_memz)
{
        return apply_function_to_matrix_values(input_ID,
                square_scalar, op_memz);
}

int main (int argc, char *argv[])
{
  initialize_larc(25, 26, 10, -1, -1, 1);

  // read in a matrix and get its matrixID. The matrix is size 8x8
  char* path_to_matrix = "Data/matrixMarketExchange1";
  int64_t A_ID = read_matrixMarketExchange_file(path_to_matrix);

  printf("\nWe compare two different methods for performing a nonlinear\n");
  printf("matrix operation: taking each element of a matrix and squaring\n");
  printf("it. The result of the two operations should be identical and\n");
  printf("therefore have the same matrixID once stored in the matrixStore.\n");

  // We call the LARC function matrix_entrySquared with scale factor
  // set to 1. This should have the same result as the local functions above.
  int64_t A2_ID = matrix_entrySquared(A_ID,"1");

  // We call the function square_scalar_matID, which takes two arguments, the
  // input matrix and a type of operation, the latter of which determines how
  // the calculation is memoized. We must use one of the operation types set
  // aside for user-defined functions, because using any operation type for
  // more than one operation can result in wrong answers.
  int64_t A2x_ID = square_matrix_elements(A_ID,FUNC_A); 

  // We compare the two matrixIDs, which should be identical
  if (A2_ID == A2x_ID) { printf("SUCCESS!\n"); }
  else
  {
     printf("ERROR in %s: the two ways of calculating the same ",__func__);
     printf("matrix return different results!\n");
  }
}
