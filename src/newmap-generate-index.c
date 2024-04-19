#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>  /* command line parsing */
#include <wchar.h>

#include "AwFmIndex.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h> /* includes stdio, string, errno, and stdlib */


int generate_fm_index(char *fastaInputFileName,
                      char *indexFileName,
                      uint8_t suffixArrayCompressionRatio,
                      uint8_t kmerLengthInSeedTable);
static PyObject* py_generate_fm_index(PyObject* self, PyObject* args) {
    char* fastaInputFileName;
    char* indexFileName;
    uint8_t suffixArrayCompressionRatio;
    uint8_t kmerLengthInSeedTable;

    if (!PyArg_ParseTuple(args, "ssBB",  /* string string u_byte u_byte */
                        &fastaInputFileName,
                        &indexFileName,
                        &suffixArrayCompressionRatio,
                        &kmerLengthInSeedTable
                        ))
        return NULL;

    enum AwFmReturnCode returnCode = generate_fm_index(
        fastaInputFileName,
        indexFileName,
        suffixArrayCompressionRatio,
        kmerLengthInSeedTable
    );

    if(returnCode == AwFmFileOpenFail){
        PyErr_SetString(PyExc_FileNotFoundError,
            "Could not open fasta file to create index\n");
        return NULL;
    }

    /* Return None in Python */
    Py_RETURN_NONE;
}

static PyMethodDef generateIndexMethods[] = {
    {"generate_fm_index", py_generate_fm_index, METH_VARARGS,
     "Generates a FM index from a fasta file."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef generateIndexModule = {
    PyModuleDef_HEAD_INIT,
    "_c_newmap_generate_index", /* name of module */
    NULL,               /* module documentation, may be NULL */
    -1,                 /* size of per-interpreter state of the module,
          or -1 if the module keeps state in global variables. */
    generateIndexMethods
};

PyMODINIT_FUNC
PyInit__c_newmap_generate_index(void) /* name is important for ref on import */
{
  return PyModule_Create(&generateIndexModule);
}


enum AwFmReturnCode generate_fm_index(char *fastaInputFileName,
                      char *indexFileName,
                      uint8_t suffixArrayCompressionRatio,
                      uint8_t kmerLengthInSeedTable) {

    struct AwFmIndex *index;

    struct AwFmIndexConfiguration config = {
            .suffixArrayCompressionRatio = suffixArrayCompressionRatio,
			.kmerLengthInSeedTable = kmerLengthInSeedTable,
			.alphabetType = AwFmAlphabetDna,
			.keepSuffixArrayInMemory = false,
			.storeOriginalSequence = false,
    };

    enum AwFmReturnCode returnCode = awFmCreateIndexFromFasta(
            &index,
            &config,
            fastaInputFileName,
            indexFileName);

    return returnCode;
}
