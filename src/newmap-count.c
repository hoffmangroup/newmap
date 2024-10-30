#include <stdint.h>
#include <stdbool.h>

#include "AwFmIndex.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h> /* includes stdio, string, errno, and stdlib */

void createIndex(const char *indexFileName, struct AwFmIndex **index) {
    enum AwFmReturnCode returnCode =
        awFmReadIndexFromFile(index, indexFileName, true);

    if(awFmReturnCodeIsFailure(returnCode)){
        PyErr_SetString(PyExc_OSError,
            "Could not load reference index from file\n");
    }
}

void createCountList(PyObject *list, struct AwFmKmerSearchList *searchList) {
    for(size_t kmerIndex = 0; kmerIndex < searchList->count; kmerIndex++) {
        const uint32_t numHitPositions = searchList->kmerSearchData[kmerIndex].count;

        PyList_SetItem(list, kmerIndex,
                       PyLong_FromUnsignedLong(numHitPositions));
    }
}

static PyObject* py_count_kmers(PyObject* self, PyObject* args) {
    const char *indexFileName;
    PyObject *inputKmerList;
    uint8_t numThreads;
    Py_ssize_t numKmers;

    if (!PyArg_ParseTuple(args, "sOb",  /* string object unsigned char */
                      &indexFileName,
                      &inputKmerList,
                      &numThreads))
        return NULL;

    if (!PyList_Check(inputKmerList)) {
        PyErr_SetString(PyExc_TypeError, "Second argument must be a list of kmer byte strings or memoryviews");
        return NULL;
    }

    numKmers = PyList_Size(inputKmerList);

    struct AwFmIndex *index;
    createIndex(indexFileName, &index);

    // NB: This returns nothing
    struct AwFmKmerSearchList *searchList =
        awFmCreateKmerSearchList(numKmers);

    searchList->count = numKmers;
    for(size_t i = 0; i < numKmers; i++) {
        // Parse a list of bytes
        PyObject *kmerBytes = PyList_GetItem(inputKmerList, i);

        if(PyBytes_AsStringAndSize(kmerBytes,
                                   &searchList->kmerSearchData[i].kmerString,
                                   &searchList->kmerSearchData[i].kmerLength)
                                   == 0) {  // NB: Returns 0 on success
            // Ensure the byte string have non-zero length
            if(searchList->kmerSearchData[i].kmerLength == 0){
                PyErr_SetString(
                    PyExc_ValueError,
                    "All elements of the kmer list must have non-zero length");
                return NULL;
            }
        }
        else {
            PyErr_SetString(
                PyExc_TypeError,
                "All elements of the kmer list must be byte strings");
            return NULL;
        }
    }

    awFmParallelSearchCount(index, searchList, numThreads);

    PyObject* returnList = PyList_New(numKmers);
    createCountList(returnList, searchList);

    // Free FM index data structures
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);

    return returnList;
}

static PyObject* py_count_kmers_from_sequence(PyObject* self, PyObject* args) {
    const char *indexFileName;
    const char *inputByteSequence;
    PyObject *inputIndexList;
    PyObject *inputLengthList;
    uint8_t numThreads;
    Py_ssize_t numKmers;

    if (!PyArg_ParseTuple(args, "syOOb",  /* string, read-only bytes-like, object, object, unsigned char */
                      &indexFileName,
                      &inputByteSequence,
                      &inputIndexList,
                      &inputLengthList,
                      &numThreads))
        return NULL;

    if (!PyList_Check(inputIndexList)) {
        PyErr_SetString(PyExc_TypeError, "Third argument must be a list of integers");
        return NULL;
    }

    if (!PyList_Check(inputLengthList)) {
        PyErr_SetString(PyExc_TypeError, "Fourth argument must be a list of integers");
        return NULL;
    }

    // Get the number of kmers to search counts for
    numKmers = PyList_Size(inputIndexList);

    // Verify the lists are the same size
    if(numKmers != PyList_Size(inputLengthList)){
        PyErr_SetString(PyExc_ValueError, "Both lists of indices and lengths must be the same length");
        return NULL;
    }

    // Create our index
    struct AwFmIndex *index;
    createIndex(indexFileName, &index);

    // NB: This returns nothing
    // Create the search list
    struct AwFmKmerSearchList *searchList =
        awFmCreateKmerSearchList(numKmers);

    // Fill the search list with the kmers
    searchList->count = numKmers;
    for(size_t i = 0; i < numKmers; i++) {
        // Get the index in the sequence from the index list
        // And get length of the corresponding kmer length from the length list
        PyObject *indexObj = PyList_GetItem(inputIndexList, i);
        PyObject *lengthObj = PyList_GetItem(inputLengthList, i);

        if (!PyLong_Check(indexObj) ||
            !PyLong_Check(lengthObj)) {
          PyErr_SetString(PyExc_TypeError, "All elements of the the index and "
                                           "length lists must be integers");
          return NULL;
        }

        size_t kmer_index = PyLong_AsSize_t(indexObj);
        size_t kmer_length = PyLong_AsSize_t(lengthObj);

        // TODO: Assert or check that the index + length does not go out of
        // bounds
        searchList->kmerSearchData[i].kmerLength = kmer_length;
        searchList->kmerSearchData[i].kmerString =
            &inputByteSequence[kmer_index];
    }

    awFmParallelSearchCount(index, searchList, numThreads);

    PyObject* returnList = PyList_New(numKmers);
    createCountList(returnList, searchList);

    // Free FM index data structures
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);

    return returnList;
}

static PyMethodDef countKmersMethods[] = {
    {"count_kmers", py_count_kmers, METH_VARARGS,
     "Returns a count of the number of times each kmer appears in the "
     "reference index."},
    {"count_kmers_from_sequence", py_count_kmers_from_sequence, METH_VARARGS,
     "Returns a count of the number of times kmer position and length appears "
     "from the provided byte sequence in the reference index."},
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
