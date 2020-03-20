//gate.h
#ifndef MPL_GATE_H
#define MPL_GATE_H

#include <inttypes.h>
#include "larc.h"
#include "global.h"
#include "matmath.h"

/* The functions in gate.c provide basic reversible computing capabilities on *
 * n-bit circuits. The CNOT and Ccnot functions return the matrixIDs        *
 * for matrices appropriately sized to the circuit which provide the function *
 * specified by the control(s) and target. The NOT function instead takes as  *
 * input an array of target values and returns the matrixID for the           *
 * matrix which applies NOT to all given targets simultaneously.              */

// Create a matrix to apply multiple NOT gates (other gates are identity)
int64_t build_not_gates(int *targets, int num_targets);

// When is_reverse_logic is 0, result is the usual CNOT gate; when it is 1, we
// instead get a controlled NOT gate that triggers when the control qubit is
// 0 instead of 1.
int64_t build_cnot_gate(int control, int target, int is_reverse_logic);

// The usual Ccnot gate
int64_t build_ccnot_gate(int control1, int control2, int target);

#endif
