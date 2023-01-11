//sycamore.c
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



// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Our header files structures and functions
#include "larc.h"
#include "matrix_store.h"
#include "matmath.h"
#include "global.h"
#include "fft.h"
#include "sycamore.h"

/*!
 * \brief These functions provide the quantum gates used in the Google vs. IBM
 * quantum supremacy debate.
 */

#ifdef IS_COMPLEX

int64_t build_sycamore_gate_sequence(char **gate_list, int system_size)
{
    // Get scalar values contained in the gates. In most cases, we must
    // initialize these values and store them in the matrixStore. There are
    // a few values that are preloaded, and for those values we can just look
    // up the matrixID for the global variable of that value.

#if 0
    scalarType scalar0_5i0_5;
    sca_init(&scalar0_5i0_5);
    sca_set_2ldoubles(&scalar0_5i0_5, 0.5L, 0.5L);
    tempPTR = get_valMatPTR_from_val(scalar0_5i0_5);
    int64_t valID_0_5i0_5 = tempPTR->packedID;
    sca_clear(&scalar0_5i0_5);

    scalarType scalar0_5iM0_5;
    sca_init(&scalar0_5iM0_5);
    sca_set_2ldoubles(&scalar0_5iM0_5, 0.5L, -0.5L);
    tempPTR = get_valMatPTR_from_val(scalar0_5iM0_5);
    int64_t valID_0_5iM0_5 = tempPTR->packedID;
    sca_clear(&scalar0_5iM0_5);

    scalarType scalarM0_5iM0_5;
    sca_init(&scalarM0_5iM0_5);
    sca_set_2ldoubles(&scalarM0_5iM0_5, -0.5L, -0.5L);
    tempPTR = get_valMatPTR_from_val(scalarM0_5iM0_5);
    int64_t valID_M0_5iM0_5 = tempPTR->packedID;
    sca_clear(&scalarM0_5iM0_5);
#endif

    int64_t valID_0_5i0_5 = get_valID_from_valString("0.5+I*0.5");
    int64_t valID_0_5iM0_5 = get_valID_from_valString("0.5-I*0.5");
    int64_t valID_M0_5iM0_5 = get_valID_from_valString("-0.5-I*0.5");

    int64_t valID_1 = packedID_scalar1;
//  1.0/sqrt(2) is not preloaded, so we must provide it. Some of our Clifford
//  scalarTypes have sqrt(2) as an algebraic extension to the rationals, so we
//  can use sca_set_enum to get a good sqrt(2). The call does not add the value
//  to the matrixStore, so we do that with get_valMatPTR_from_val.
    sca_set_enum(&(scratchVars.quick_use),SCALAR_ENUM_INV_SQRT2);
    mats_ptr_t valPTR_inv_sqrt_2 = get_scalarPTR_for_scalarVal(scratchVars.quick_use);
    int64_t valID_inv_sqrt_2 = valPTR_inv_sqrt_2->packedID;
    mats_ptr_t valPTR_0iM1 = get_scalarPTR_for_scalarVal(scalar0iM1);
    int64_t valID_0iM1 = valPTR_0iM1->packedID;
//  kronecker product of scalars is equivalent to a scalar multiply
    int64_t valID_0iM1_div_sqrt_2 = kronecker_product(valID_0iM1,
         valID_inv_sqrt_2);
    
    // get pointers to addresses for required 2x2 matrices
    int64_t packedID_I = get_identity_pID(1);

    int64_t panels_sqrtX[4] = { valID_0_5i0_5,  valID_0_5iM0_5,
                                  valID_0_5iM0_5, valID_0_5i0_5 };
    int64_t packedID_sqrtX = get_pID_from_four_sub_pIDs(panels_sqrtX[0],
       panels_sqrtX[1], panels_sqrtX[2], panels_sqrtX[3], 1, 1);

    int64_t panels_sqrtY[4] = { valID_0_5i0_5,  valID_M0_5iM0_5,
                                  valID_0_5i0_5, valID_0_5i0_5 };
    int64_t packedID_sqrtY = get_pID_from_four_sub_pIDs(panels_sqrtY[0],
       panels_sqrtY[1], panels_sqrtY[2], panels_sqrtY[3], 1, 1);

    int64_t panels_sqrtW[4] = { valID_0_5i0_5,    valID_0iM1_div_sqrt_2,
                                  valID_inv_sqrt_2, valID_0_5i0_5 };
    int64_t packedID_sqrtW = get_pID_from_four_sub_pIDs(panels_sqrtW[0],
       panels_sqrtW[1], panels_sqrtW[2], panels_sqrtW[3], 1, 1);

    // iteratively construct 2^n by 2^n matrix starting with scalar 1 value
    int64_t out_pID = valID_1;

    int64_t packedID_loop;
    for (int index = system_size-1; index >= 0; --index)
    {
        if (0 == strncmp(gate_list[index], "I", strlen("I")+1)) {
            packedID_loop = packedID_I;
        } else if (0 == strncmp(gate_list[index], "sqrtX", strlen("sqrtX")+1)) {
            packedID_loop = packedID_sqrtX;
        } else if (0 == strncmp(gate_list[index], "sqrtY", strlen("sqrtY")+1)) {
            packedID_loop = packedID_sqrtY;
        } else if (0 == strncmp(gate_list[index], "sqrtW", strlen("sqrtW")+1)) {
            packedID_loop = packedID_sqrtW;
        } else {
            printf("ERROR in %s: unrecognized gate (%s) at index %d.\n",
                __func__, gate_list[index], index);
            return -1;
        }
        out_pID = kronecker_product(packedID_loop, out_pID);
    }
    return out_pID;
}

int64_t build_sycamore_2gate(int target1, int target2, int system_size)
{

  if ((target1 < 0) || (target1 >= system_size))
  {
    printf("ERROR in %s: target1 (= %d) out of valid range 0 to %d.\n",
        __func__, target1, system_size - 1);
    return -1;
  }

  if ((target2 < 0) || (target2 >= system_size))
  {
    printf("ERROR in %s: target2 (= %d) out of valid range 0 to %d.\n",
        __func__, target2, system_size - 1);
    return -1;
  }

  if (target1 == target2)
  {
    printf("ERROR in %s: target1 (= %d) and target2 (= %d) cannot have the same value.\n",
        __func__, target1, target2);
    return -1;
  }

  // Because the supremacy gate is symmetric between the qubits,
  // we can sort the targets.
  int min_target, max_target;
  if (target1 < target2) {
      min_target = target1;
      max_target = target2;
  } else {
      min_target = target2;
      max_target = target1;
  }

  // First, we build the submatrices at level system_size - max_target - 1.
  int stage1 = system_size - max_target - 1;
  int64_t mat_stage1_zero_pID = get_zero_pID(stage1, stage1);
  int64_t mat_stage1_id_pID = get_identity_pID(stage1);
  int64_t val_1j_pID = get_scalarPTR_for_scalarVal(scalar0i1)->packedID;
  int64_t mat_stage1_1j_pID = scalar_mult(val_1j_pID, mat_stage1_id_pID);
  scalarType omega_scalar;
  sca_init(&omega_scalar);
  sca_set_enum(&omega_scalar, SCALAR_ENUM_SQRT3);
  sca_add(&omega_scalar, omega_scalar, scalar0i1);
  sca_mult(&omega_scalar, omega_scalar, scalar0_5);
  int64_t val_omega_pID = get_scalarPTR_for_scalarVal(omega_scalar)->packedID;
  sca_clear(&omega_scalar);
  int64_t mat_stage1_omega_pID = scalar_mult(val_omega_pID, mat_stage1_id_pID);

  // Second, we build the submatrices at level system_size - max_target.
  int stage2 = stage1 + 1;
  int64_t panels_stage2_00[4] = { mat_stage1_id_pID,   mat_stage1_zero_pID,
                                    mat_stage1_zero_pID, mat_stage1_zero_pID };
  int64_t mat_stage2_00_pID = get_pID_from_four_sub_pIDs(
     panels_stage2_00[0], panels_stage2_00[1], panels_stage2_00[2],
     panels_stage2_00[3], stage2, stage2);

  int64_t panels_stage2_01[4] = { mat_stage1_zero_pID, mat_stage1_zero_pID,
                                    mat_stage1_1j_pID,   mat_stage1_zero_pID };
  int64_t mat_stage2_01_pID = get_pID_from_four_sub_pIDs(
     panels_stage2_01[0], panels_stage2_01[1], panels_stage2_01[2],
     panels_stage2_01[3], stage2, stage2);

  int64_t panels_stage2_10[4] = { mat_stage1_zero_pID, mat_stage1_1j_pID,
                                    mat_stage1_zero_pID, mat_stage1_zero_pID };
  int64_t mat_stage2_10_pID = get_pID_from_four_sub_pIDs(
     panels_stage2_10[0], panels_stage2_10[1], panels_stage2_10[2],
     panels_stage2_10[3], stage2, stage2);

  int64_t panels_stage2_11[4] = { mat_stage1_zero_pID, mat_stage1_zero_pID,
                                    mat_stage1_zero_pID, mat_stage1_omega_pID };
  int64_t mat_stage2_11_pID = get_pID_from_four_sub_pIDs(
     panels_stage2_11[0], panels_stage2_11[1], panels_stage2_11[2],
     panels_stage2_11[3], stage2, stage2);

  // Third, we build the submatrices at level system_size - min_target - 1.
  int stage3 = system_size - min_target - 1;
  int64_t mat_stage3_expander_pID = get_identity_pID(stage3 - stage2);
  int64_t mat_stage3_00_pID = kronecker_product(mat_stage3_expander_pID, mat_stage2_00_pID);
  int64_t mat_stage3_01_pID = kronecker_product(mat_stage3_expander_pID, mat_stage2_01_pID);
  int64_t mat_stage3_10_pID = kronecker_product(mat_stage3_expander_pID, mat_stage2_10_pID);
  int64_t mat_stage3_11_pID = kronecker_product(mat_stage3_expander_pID, mat_stage2_11_pID);

  // Fourth, we build the submatrices at level system_size - min_target.
  int stage4 = stage3 + 1;
  int64_t panels_stage4[4] = { mat_stage3_00_pID, mat_stage3_01_pID,
                                 mat_stage3_10_pID, mat_stage3_11_pID };
  int64_t mat_stage4_pID = get_pID_from_four_sub_pIDs( panels_stage4[0],
      panels_stage4[1], panels_stage4[2], panels_stage4[3], stage4, stage4);

  // Lastly, we build the matrix at level system_size.
  int64_t mat_stage5_expander_pID = get_identity_pID(system_size - stage4);
  int64_t mat_supremacy_pID = kronecker_product(mat_stage5_expander_pID, mat_stage4_pID);

  return mat_supremacy_pID;
}

#endif

