//gate.c 

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

/* The functions in gate.c provide basic reversible computing capabilities on *
 * n-bit circuits. The CNOT and CCNOT functions return the matrixIDs          *
 * for matrices appropriately sized to the circuit which provide the function *
 * specified by the control(s) and target. The NOT function instead takes as  *
 * input an array of target values and returns the matrixID for the           *
 * matrix which applies NOT to all given targets simultaneously.              */

// Create a matrix to apply multiple NOT gates (other gates are identity)
int64_t build_not_gates(int* targets, int num_targets)
{

  int64_t max_level = maximum_level();

  // create bit representation of targets (to avoid need to sort target list)
  int64_t bitflags = 0;
  for (int i = 0; i < num_targets; ++i) { bitflags ^= (1 << targets[i]); }
  
  // get pointers to addresses for required 2x2 matrices
  mat_ptr_t matptr_I1 = get_identity_matrix_ptr(1);
  mat_ptr_t matptr_NOT = get_matPTR_from_matID(matID_NOT, "NOT matrix", __func__, 0);

  // recursively construct 2^n by 2^n matrix starting with scalar 1 value
  mat_ptr_t out_ptr = get_matPTR_from_matID(matID_scalar1, "scalar 1", __func__, 0);
  for (int bit_index = max_level-1; bit_index >= 0; --bit_index)
  {
    mat_ptr_t matptr_loop = ((bitflags>>bit_index)&1) ? matptr_NOT : matptr_I1;
    out_ptr = kronecker_product(matptr_loop, out_ptr);
  }
  return out_ptr->matrixID;
}

// When is_reverse_logic is 0, result is the normal CNOT; when it is 1, we
// instead get a controlled NOT gate that triggers when the control qubit is
// 0 instead of 1.
int64_t build_cnot_gate(int control, int target, int is_reverse_logic)
{

  if (target == control)
  {
    printf("ERROR in %s: control and target cannot have the same value.\n",
        __func__);
    return -1;
  }

  int64_t max_level = maximum_level();
  mat_ptr_t tmp0_ptr, tmp1_ptr;

  // get pointers to addresses for required 2x2 matrices
  mat_ptr_t matptr_I1 = get_identity_matrix_ptr(1);
  mat_ptr_t matptr_NOT = get_matPTR_from_matID(matID_NOT, "NOT matrix", __func__, 0);
  // When t0_ptr is I1 and t1_ptr is NOT, construct a controlled NOT gate.
  // When t1_ptr is I1 and t0_ptr is NOT, construct a "reverse-logic" controlled
  //  NOT gate that triggers when the control qubit is 0 instead of 1.
  mat_ptr_t t0_ptr = (is_reverse_logic) ? matptr_NOT : matptr_I1;
  mat_ptr_t t1_ptr = (is_reverse_logic) ? matptr_I1 : matptr_NOT;

  mat_ptr_t panel[4];
  panel[0] = get_matPTR_from_matID(matID_scalar1, "scalar 1", __func__, 0);
  panel[1] = panel[2] = panel[3] = get_matPTR_from_matID(matID_scalar0, "scalar 0", __func__, 0);
  mat_ptr_t matptr_P00 = get_matPTR_from_array_of_four_subMatPTRs(panel,1,1);
  panel[3] = panel[0];
  panel[0] = panel[1];
  mat_ptr_t matptr_P11 = get_matPTR_from_array_of_four_subMatPTRs(panel,1,1);

  // recursively construct 2^n x 2^n matrix starting with scalar 1 value
  tmp0_ptr = tmp1_ptr = get_matPTR_from_matID(matID_scalar1, "scalar 1", __func__, 0);
  for (int bit_index = max_level-1; bit_index >= 0; --bit_index)
  {
    if (bit_index == control)
    {
      tmp0_ptr = kronecker_product(matptr_P00, tmp0_ptr);
      tmp1_ptr = kronecker_product(matptr_P11, tmp1_ptr);
    }
    else if (bit_index == target)
    {
      tmp0_ptr = kronecker_product(t0_ptr, tmp0_ptr);
      tmp1_ptr = kronecker_product(t1_ptr, tmp1_ptr);
    }
    else
    {
      tmp0_ptr = kronecker_product(matptr_I1, tmp0_ptr);
      tmp1_ptr = kronecker_product(matptr_I1, tmp1_ptr);
    }
  }

  // CNOT is sum of the two temporary matrices
  mat_ptr_t out_ptr = matrix_add(tmp0_ptr, tmp1_ptr);
  return out_ptr->matrixID;

}

int64_t build_ccnot_gate(int control1, int control2, int target)
{
  mat_level_t max_level = maximum_level();
  // printf("The maximum_level function inside the build_ccnot_gate returns %d\n",max_level);

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

  mat_ptr_t  tmp01_ptr, tmp10_ptr, tmp11_ptr;
  mat_ptr_t  tmp1x_sum_ptr;

  // get pointers to addresses for required 2x2 matrices
  mat_ptr_t matptr_I1 = get_identity_matrix_ptr(1);
  mat_ptr_t matptr_NOT = get_matPTR_from_matID(matID_NOT, "NOT matrix", __func__, 0);

  mat_ptr_t panel[4];
  panel[0] = get_matPTR_from_matID(matID_scalar1, "scalar 1", __func__, 0);
  panel[1] = panel[2] = panel[3] = get_matPTR_from_matID(matID_scalar0, "scalar 0", __func__, 0);
  mat_ptr_t matptr_P00 = get_matPTR_from_array_of_four_subMatPTRs(panel,1,1);
  panel[3] = panel[0];
  panel[0] = panel[1];
  mat_ptr_t matptr_P11 = get_matPTR_from_array_of_four_subMatPTRs(panel,1,1);

  // recursively construct Ccnot gate starting with scalar 1 value
  tmp01_ptr = tmp10_ptr = tmp11_ptr = get_matPTR_from_matID(matID_scalar1, "scalar 1", __func__, 0);
  for (int bit_index = max_level-1; bit_index >= 0; --bit_index)
  {
    if (bit_index == control1)
    {
      tmp01_ptr = kronecker_product(matptr_P00, tmp01_ptr);
      tmp10_ptr = kronecker_product(matptr_P11, tmp10_ptr);
      tmp11_ptr = kronecker_product(matptr_P11, tmp11_ptr);
    }
    else if (bit_index == control2)
    {
      tmp01_ptr = kronecker_product(matptr_I1,  tmp01_ptr);
      tmp10_ptr = kronecker_product(matptr_P00, tmp10_ptr);
      tmp11_ptr = kronecker_product(matptr_P11, tmp11_ptr);
    }
    else if (bit_index == target)
    {
      tmp01_ptr = kronecker_product(matptr_I1,  tmp01_ptr);
      tmp10_ptr = kronecker_product(matptr_I1,  tmp10_ptr);
      tmp11_ptr = kronecker_product(matptr_NOT, tmp11_ptr);
    }
    else
    {
      tmp01_ptr = kronecker_product(matptr_I1, tmp01_ptr);
      tmp10_ptr = kronecker_product(matptr_I1, tmp10_ptr);
      tmp11_ptr = kronecker_product(matptr_I1, tmp11_ptr);
    }
  }

  // Ccnot is the sum of the three temporary matrices
  tmp1x_sum_ptr = matrix_add(tmp10_ptr, tmp11_ptr);
  mat_ptr_t out_ptr = matrix_add(tmp01_ptr, tmp1x_sum_ptr);
  return out_ptr->matrixID;

}

