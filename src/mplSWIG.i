%module mplSWIG

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

/*
  The following functions all return char * (string) type from malloc'd
  memory and are declared in headers in this .i file. This is how you
  tell SWIG it's ok for the wrapper functions to free the memory.
*/
%newobject info_get;                                  /* info_store.h */
%newobject create_log_dir;                            /* larc.h */
%newobject sca_get_str;                               /* larc.h */
%newobject tracenorm;                                 /* matmath.h */
%newobject trace;                                     /* matmath.h */
%newobject get_scalar_value_string;                   /* matmath.h */
%newobject matrix_count_entries;                      /* matmath.h */
%newobject get_list_of_scalars_in_larcMatrixFile;     /* matmath.h */
%newobject get_readableString_scalar_from_pID_and_coords;/* matrix_store.h */

%{
#include "../larc/src/type.h"
#include "../larc/src/larc.h"
#include "../larc/src/fft.h"
#include "../larc/src/global.h"
#include "../larc/src/hash.h"
#include "../larc/src/info_store.h"
#include "../larc/src/io.h"
#include "../larc/src/matmath.h"
#include "../larc/src/matrix_store.h"
#include "../larc/src/op_store.h"
#include "../larc/src/organize.h"
#include "../larc/src/scalars.h"
#include "version.h"
#include "gate.h"
#include "sycamore.h"
#include <complex.h>
#include <gmp.h>
#include <pthread.h>
%}

#define SWIGWORDSIZE64
// #enddef

%include "stdint.i"
%include "complex.i"
%include "cstring.i"
%include "carrays.i"
%include "cpointer.i"

/* This tells SWIG to output int64_t ** as list of integer lists. 
   In this setup, -1's mark the end of each list and a NULL marks the end
   of the final list */
/* As of now, locate_entries_larcMatrixFile is the only routine that makes use of this.*/
%typemap(out) int64_t ** {
    int len;
    int i, j;
    len = 0;
    while (NULL != $1[len]) len++;
    $result = PyList_New(len);
    for (i = 0; i < len; i++){
        j = 0;
        PyObject *l = PyList_New(0);
        while (0 <= $1[i][j]){
            PyList_Append(l, PyLong_FromLong($1[i][j]));
            j ++;
        }
        PyList_SetItem($result, i, l);
    }
}

/* This tells SWIG to input list of strings as char **. */
%typemap(in) char ** {
    /* Check if is a list */
    if (PyList_Check($input)){
        int size = PyList_Size($input);
        int i = 0;
        $1 = (char **) malloc((size+1)*sizeof(char *));
        for (i = 0; i < size; i++){
            PyObject *o = PyList_GetItem($input, i);
            if (PyUnicode_Check(o))
                $1[i] = PyBytes_AsString(PyUnicode_AsEncodedString(o, "utf-8", "strict"));
            else {
                PyErr_SetString(PyExc_TypeError, "list must contain strings");
                free($1);
                return NULL;
            }
        }
        $1[i] = 0;
    }
    else {
        PyErr_SetString(PyExc_TypeError, "argument is not a list");
        return NULL;
    }
}

/* This cleans up the char ** array we malloc'd before the function call. */
%typemap(freearg) char ** {
    free((char *) $1);
}


// # NOTE: global.h uses the USE_INTEGER/REAL/COMPLEX that is #defined in
// # type.h, so type.h must be before it in the list

%include "../larc/src/type.h"
%include "../larc/src/larc.h"
%include "../larc/src/fft.h"
%include "../larc/src/global.h"
%include "../larc/src/hash.h"
%include "../larc/src/info_store.h"
%include "../larc/src/io.h"
%include "../larc/src/matmath.h"
%include "../larc/src/matrix_store.h"
%include "../larc/src/op_store.h"
%include "../larc/src/organize.h"
%include "../larc/src/scalars.h"
%include "version.h"
%include "gate.h"
%include "sycamore.h"

%array_class(complex, complexArray);
%array_class(long int, int64Array); // works because SWIGWORDSIZE64 defined
%array_class(int, intArray);
%array_class(double, doubleArray);
