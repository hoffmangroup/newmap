#include <stdint.h>
#include <stdbool.h>

#include "AwFmIndex.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h> /* includes stdio, string, errno, and stdlib */

static PyObject* py_count_kmers(PyObject* self, PyObject* args) {
    const char *indexFileName;
    PyObject *inputKmerList;
    uint8_t numThreads;
    Py_ssize_t numKmers;

    if (!PyArg_ParseTuple(args, "sOb",  /* string object predicate */
                      &indexFileName,
                      &inputKmerList,
                      &numThreads))
        return NULL;

    if (!PyList_Check(inputKmerList)) {
        PyErr_SetString(PyExc_TypeError, "Second argument must be a list of kmer byte strings");
        return NULL;
    }

    numKmers = PyList_Size(inputKmerList);

    // Create the FM index to query from file
    struct AwFmIndex *index;
    enum AwFmReturnCode returnCode = 
        awFmReadIndexFromFile(&index, indexFileName, true);

    if(awFmReturnCodeIsFailure(returnCode)){
        PyErr_SetString(PyExc_OSError,
            "Could not load reference index from file\n");
        return NULL;
    }

    // NB: This returns nothing
    struct AwFmKmerSearchList *searchList = 
        awFmCreateKmerSearchList(numKmers);

    searchList->count = numKmers;
    for(size_t i = 0; i < numKmers; i++) {
        PyObject *kmerBytes = PyList_GetItem(inputKmerList, i);
        // Ensure all elements are byte strings to avoid seg faults
        if (!PyBytes_Check(kmerBytes)) {
            PyErr_SetString(PyExc_TypeError,
                    "All elements of the kmer list must be byte strings");
            return NULL;
        }

        searchList->kmerSearchData[i].kmerLength = PyBytes_Size(kmerBytes);
        // NB: The kmer strings are never deallocated via freeing the kmer
        // search list, so we can use as-is
        searchList->kmerSearchData[i].kmerString = PyBytes_AsString(kmerBytes);
    }

    awFmParallelSearchCount(index, searchList, numThreads);

    PyObject* returnList = PyList_New(numKmers);
    for(size_t kmerIndex = 0; kmerIndex < searchList->count; kmerIndex++) {
        const uint32_t numHitPositions = searchList->kmerSearchData[kmerIndex].count;
        
        PyList_SetItem(returnList, kmerIndex,
                       PyLong_FromUnsignedLong(numHitPositions));
    }
    
    // Free FM index data structures
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);

    return returnList;
}

static PyMethodDef countKmersMethods[] = {
    {"count_kmers", py_count_kmers, METH_VARARGS,
     "Returns a count of the number of times each kmer appears in the "
     "reference index."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef countKmersModule = {
    PyModuleDef_HEAD_INIT,
    "_c_newmap_count_kmers", /* name of module */
    NULL,               /* module documentation, may be NULL */
    -1,                 /* size of per-interpreter state of the module,
          or -1 if the module keeps state in global variables. */
    countKmersMethods
};

PyMODINIT_FUNC
PyInit__c_newmap_count_kmers(void) /* name is important for ref on import */
{
    return PyModule_Create(&countKmersModule);
}
