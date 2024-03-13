//                    sycamore_testrun.c
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

int main()
{

#ifndef IS_COMPLEX
    printf("recompile MyPyLARC with a complex type and run again...\n");
    exit(0);
#endif

    // paramaters for large machine copied from Python default.init_params
    int matrix_exponent = 31;
    int op_exponent = 30;
    int regionbitparam = -1;
    int zeroregionbitparam = -1;
    int larc_verbose = 3;
    int system_size = 12;
//    int report_period = 3600;
    initialize_larc(matrix_exponent,op_exponent,12,
        regionbitparam,zeroregionbitparam,larc_verbose);
//    create_report_thread(report_period);

    printf("\n    ######################################################################");
    printf("    ##  Create matrices for Sycamore quantum supremacy computation.     ##");
    printf("    ######################################################################\n");
    
    int no_cycles_to_use = 4;
    // will read in cycle matrices generated by previous Python code run
    uint64_t circuitMatrixID = get_identity_pID(system_size);
    uint64_t cycleMatrixID;
    char cycleMatrix[50];
    for (int i=0;i<no_cycles_to_use;++i)
    {
        printf("Multiplying cycle matrix %d into circuit...",i);
        sprintf(cycleMatrix,"matrixFiles/sycamore_cycle_%d_matrix.json",i);
        cycleMatrixID = read_larcMatrixFile(cycleMatrix);
        circuitMatrixID = matrix_mult(circuitMatrixID, cycleMatrixID);
//        empty_op_store();
    }
    sprintf(cycleMatrix,"matrixFiles/cycle%d_product_matrix.json",no_cycles_to_use);
    fprint_larcMatrixFile(circuitMatrixID,cycleMatrix);
}

