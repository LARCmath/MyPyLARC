//gate.c 
/******************************************************************
 *                                                                *
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
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



// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Our header files structures and functions
// #include "../larc/src/global.h"
// #include "../larc/src/matrix_store.h"
#include "gate.h"

/*!
 * \file gate.c
 * \brief These functions provide basic reversible computing capabilities on n-bit circuits.
 *
 * The CNOT and CCNOT functions return the matrixIDs for matrices
 * appropriately sized to the circuit which provide the function
 * specified by the control(s) and target. The NOT function also produces
 * a matrix of the correct size, but instead takes as input an array of
 * target values and returns the matrixID for the matrix which applies
 * NOT to all given targets simultaneously.
 */

int64_t build_not_gates(int* targets, int num_targets)
{

  uint64_t max_level = max_level_allowed_matrixStore();

  // create bit representation of targets (to avoid need to sort target list)
  uint64_t bitflags = 0;
  for (int i = 0; i < num_targets; ++i)
  { bitflags ^= ((uint64_t)1 << targets[i]); }
  
  // LARC has global variables for the matrixIDs of the 2x2 identity (packedID_I1)
  // and NOT matrix (packedID_NOT). We also need the scalar1 for the initial state
  // of out_pID.

  int64_t out_pID = packedID_scalar1;

  // recursively construct 2^n by 2^n matrix starting with scalar 1 value
  for (int bit_index = max_level-1; bit_index >= 0; --bit_index)
  {
    int64_t pID_loop = ((bitflags>>bit_index)&1) ? packedID_NOT : packedID_I1;
    out_pID = kronecker_product(pID_loop, out_pID);
  }
  return out_pID;
}

int64_t build_cnot_gate(int control, int target, int is_reverse_logic)
{

  if (target == control)
  {
    printf("ERROR in %s: control and target cannot have the same value.\n",
        __func__);
    return -1;
  }

  int64_t max_level = max_level_allowed_matrixStore();
  int64_t tmp0_pID, tmp1_pID;

  // LARC has global variables for the matrixIDs of the 2x2 identity (packedID_I1)
  // and NOT matrix (packedID_NOT). We also need packedID_scalar0 and packedID_scalar1
  // for referencing the scalar values 0 and 1.

  // When t0_ptr is I1 and t1_ptr is NOT, construct a controlled NOT gate.
  // When t1_ptr is I1 and t0_ptr is NOT, construct a "reverse-logic" controlled
  //  NOT gate that triggers when the control qubit is 0 instead of 1.
  int64_t t0_pID = (is_reverse_logic) ? packedID_NOT : packedID_I1;
  int64_t t1_pID = (is_reverse_logic) ? packedID_I1 : packedID_NOT;

  int64_t panel[4];
  panel[0] = packedID_scalar1;
  panel[1] = panel[2] = panel[3] = packedID_scalar0;
  int64_t packedID_P00 = get_pID_from_four_sub_pIDs(panel[0], panel[1],
          panel[2], panel[3], 1, 1);
  panel[3] = panel[0];
  panel[0] = panel[1];
  int64_t packedID_P11 = get_pID_from_four_sub_pIDs(panel[0], panel[1],
          panel[2], panel[3], 1, 1);

  // recursively construct 2^n x 2^n matrix starting with scalar 1 value
  tmp0_pID = tmp1_pID = packedID_scalar1;

  for (int bit_index = max_level-1; bit_index >= 0; --bit_index)
  {
    if (bit_index == control)
    {
      tmp0_pID = kronecker_product(packedID_P00, tmp0_pID);
      tmp1_pID = kronecker_product(packedID_P11, tmp1_pID);
    }
    else if (bit_index == target)
    {
      tmp0_pID = kronecker_product(t0_pID, tmp0_pID);
      tmp1_pID = kronecker_product(t1_pID, tmp1_pID);
    }
    else
    {
      tmp0_pID = kronecker_product(packedID_I1, tmp0_pID);
      tmp1_pID = kronecker_product(packedID_I1, tmp1_pID);
    }
  }

  // CNOT is sum of the two temporary matrices
  return matrix_add(tmp0_pID, tmp1_pID);
}

int64_t build_ccnot_gate(int control1, int control2, int target)
{
  mat_level_t max_level = max_level_allowed_matrixStore();
  //printf("The max_level_allowed_matrixStore function inside the build_ccnot_gate returns %d\n",max_level);

  if (target == control1)
  {
    printf("ERROR in %s: control1 and target cannot have the same value.\n",
        __func__);
    return -1;
  }

  if (target == control2)
  {
    printf("ERROR in %s: control2 and target cannot have the same value.\n",
        __func__);
    return -1;
  }

  if (control1 == control2)
  {
    printf("ERROR in %s: control1 and control2 cannot have the same value.\n",
        __func__);
    return -1;
  }

  // LARC has global variables for the matrixIDs of the 2x2 identity (packedID_I1)
  // and NOT matrix (packedID_NOT). We also need packedID_scalar0 and packedID_scalar1
  // for referencing the scalar values 0 and 1.

  int64_t  tmp01_pID, tmp10_pID, tmp11_pID;
  int64_t  tmp1x_sum_pID;

  int64_t panel[4];
  panel[0] = packedID_scalar1;
  panel[1] = panel[2] = panel[3] = packedID_scalar0;
  int64_t packedID_P00 = get_pID_from_four_sub_pIDs(panel[0], panel[1],
             panel[2], panel[3], 1, 1);
  panel[3] = panel[0];
  panel[0] = panel[1];
  int64_t packedID_P11 = get_pID_from_four_sub_pIDs(panel[0], panel[1],
             panel[2], panel[3], 1, 1);

  // recursively construct Ccnot gate starting with scalar 1 value
  tmp01_pID = tmp10_pID = tmp11_pID = packedID_scalar1;
  for (int bit_index = max_level-1; bit_index >= 0; --bit_index)
  {
    if (bit_index == control1)
    {
      tmp01_pID = kronecker_product(packedID_P00, tmp01_pID);
      tmp10_pID = kronecker_product(packedID_P11, tmp10_pID);
      tmp11_pID = kronecker_product(packedID_P11, tmp11_pID);
    }
    else if (bit_index == control2)
    {
      tmp01_pID = kronecker_product(packedID_I1,  tmp01_pID);
      tmp10_pID = kronecker_product(packedID_P00, tmp10_pID);
      tmp11_pID = kronecker_product(packedID_P11, tmp11_pID);
    }
    else if (bit_index == target)
    {
      tmp01_pID = kronecker_product(packedID_I1,  tmp01_pID);
      tmp10_pID = kronecker_product(packedID_I1,  tmp10_pID);
      tmp11_pID = kronecker_product(packedID_NOT, tmp11_pID);
    }
    else
    {
      tmp01_pID = kronecker_product(packedID_I1, tmp01_pID);
      tmp10_pID = kronecker_product(packedID_I1, tmp10_pID);
      tmp11_pID = kronecker_product(packedID_I1, tmp11_pID);
    }
  }

  // Ccnot is the sum of the three temporary matrices
  tmp1x_sum_pID = matrix_add(tmp10_pID, tmp11_pID);
  return matrix_add(tmp01_pID, tmp1x_sum_pID);

}

