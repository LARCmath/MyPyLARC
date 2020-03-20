%module mplSWIG

%{
#include "../larc/src/type.h"
#include "../larc/src/fft.h"
#include "../larc/src/experimental.h"
#include "../larc/src/global.h"
#include "../larc/src/io.h"
#include "../larc/src/matmath.h"
#include "../larc/src/larc.h"
#include "../larc/src/info_store.h"
#include "../larc/src/matrix_store.h"
#include "../larc/src/op_store.h"
#include "../larc/src/organize.h"
#include "../larc/src/scalars.h"
#include "gate.h"
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
%include "../larc/src/fft.h"
%include "../larc/src/experimental.h"
%include "../larc/src/global.h"
%include "../larc/src/io.h"
%include "../larc/src/matmath.h"
%include "../larc/src/larc.h"
%include "../larc/src/info_store.h"
%include "../larc/src/matrix_store.h"
%include "../larc/src/op_store.h"
%include "../larc/src/organize.h"
%include "../larc/src/scalars.h"
%include "gate.h"

%array_class(complex, complexArray);
%array_class(long int, int64Array); // works because SWIGWORDSIZE64 defined
%array_class(int, intArray);
%array_class(double, doubleArray);
