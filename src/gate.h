//gate.h
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

/*!
 * \brief Create a matrix which applies multiple NOT gates to a circuit
 *
 * This routine takes a list of target wires and applies a NOT gate to each
 * target. Since all NOT gates commute, doing them individually would be
 * inefficient.
 *
 * \param targets An array of target wire indices
 * \param num_targets The size of the input array
 * \result The matrixID of the requested gate matrix
 */
int64_t build_not_gates(int *targets, int num_targets);

/*!
 * \brief Create a matrix which applies a controlled-NOT gate to a circuit
 *
 * This routine produces a matrix which applies a single controlled-NOT gate
 * to a circuit. Only single CNOTs are created as not all CNOT gates commute.
 * When is_reverse_logic is 0, result is the normal CNOT, where the target
 * wire is NOTted only when the control wire has value 1. If is_reverse_logic
 * is 1, we instead get a controlled NOT gate that triggers when the control
 * qubit is 0 instead of 1.
 *
 * \param control The index of the wire which controls the NOT gate
 * \param target The index of the wire to which the NOT gate may be applied
 * \param is_reverse_logic When 1, NOT triggers on a 0 instead of a 1
 * \result The matrixID for the gate matrix
 *
 */
int64_t build_cnot_gate(int control, int target, int is_reverse_logic);

// The usual Ccnot gate
/*!
 * \brief Create a matrix which applies a doubly-controlled-NOT gate to a circuit
 *
 * This routine produces a matrix which applies a single doubly-controlled-NOT
 * gate to a circuit. Only single CCNOTs are created as not all CCNOT gates
 * commute. The target wire is NOTted only if both control wires have value 1.
 *
 * \param control1 The index of one wire which controls the NOT gate
 * \param control2 The index of the other wire which controls the NOT gate
 * \param target The index of the wire to which the NOT gate may be applied
 * \result The matrixID for the gate matrix
 *
 */
int64_t build_ccnot_gate(int control1, int control2, int target);


#endif
