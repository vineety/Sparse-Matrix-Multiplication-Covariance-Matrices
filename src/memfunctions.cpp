#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "matrix_def.h"

void destroy(struct sparsemat *matrix)
{
    /*This routine deallocates dynamically
    allocated memory to a sparsemat in a
    CSR format*/
    if (matrix != 0)
        {
            if (matrix->rowPtr != 0)
                {
                    free(matrix->rowPtr);
                    matrix->rowPtr = NULL;
                }

            if (matrix->colInd != 0)
                {
                    free(matrix->colInd);
                    matrix->colInd = NULL;
                }

            if (matrix->values != 0)
                {
                   free(matrix->values);
                   matrix->values = NULL;
                }
            matrix->nzmax = 0;
            matrix->rows = 0;
            matrix->cols = 0;

        }
}

void destroy(struct iarray *matrix)
{
    /*This routine deallocates
    dynamically allocated memory to
    an integer array (i4arry) */
    if (matrix != NULL)
        {
            if (matrix->array != NULL)
                {
                    free(matrix->array);
                    matrix->array = NULL;
                }
            matrix->rows = 0;
            matrix->cols = 0;
            matrix = NULL;
        }
}
void destroy(struct darray *matrix)
{
    /*This routine deallocates
    dynamically allocated memory to
    an integer array (iarray) */
    if (matrix != NULL)
        {
            if (matrix->array != NULL)
                {
                    free(matrix->array);
                    matrix->array = NULL;
                }
            matrix->rows = 0;
            matrix->cols = 0;
            matrix = NULL;
        }
}


void modifyalloc(struct sparsemat *matrix, int imemSize)
{
    /**************************************************************************
    Written by: Vineet Yadav
    Date: 12/15/2014

    This function encapsulates (wrapper) for the realloc function.
    Its main purpose is to increase or decrease the size of allocated memory
    by calloc or malloc. Here, it is used to increase or decrease the size of
    ColInd and value vectors of the output sparse matrix in CSR format as the number of
    non-zeros in the output matrix is unknown. The function uses two pointers for
    increasing the size of memory and in the end set them to NULL and return the
    original pointers pointing to the memory with increased size.
    Input:
    (1) matrix: matrixa in CSR format
    (2) imemSize: size of the readjusted memory that is by how many numbers should
    the memory of matrix should be increased.
    Thus if:
    imemsize = 200 then original allocation would be increased or decreased and would
    now have space for storing 200 double or int values.
    ****************************************************************************/


    // get out of the function if imemSize is <= 0 and
    // set rows, cols and nzmax of the matrix to zero
    if (imemSize <= 0)
        {
            free(matrix->rowPtr);
            matrix->rowPtr = NULL;
            free(matrix->colInd);
            matrix->colInd = NULL;
            free(matrix->values);
            matrix->values = NULL;
            matrix->rows = 0;
            matrix->cols = 0;
            matrix->nzmax = 0;
            matrix = NULL;
            return;
        }
    // Only implement this section if initially allocated pointer is
    // not allocated to the desired size
    if (imemSize != matrix->nzmax)
        {
            int *colA = NULL;
            //int *rowA = NULL;
            double *valA = NULL;

            //rowA = (int*)realloc(matrix->rowPtr, matrix->rows*sizeof(int));

            colA = (int*)realloc(matrix->colInd, imemSize * sizeof(int));
            valA = (double*)realloc(matrix->values, imemSize * sizeof(double));

            // Exchange pointers and set pointers used for reallocation to NULL
            matrix->colInd = colA;
            colA = NULL;
            matrix->values = valA;
            valA = NULL;
            //matrix->rowPtr = rowA;   rowA = NULL;
        }
}

struct sparsemat *create(int rows, int cols, int nzmax)
{
    struct sparsemat* matrix;
    matrix = (struct sparsemat*) calloc(1, sizeof(struct sparsemat));
    matrix->rows = rows;
    matrix->nzmax = nzmax;
    matrix->cols = cols;
    matrix->colInd = (int*)calloc(nzmax, sizeof(int));
    matrix->rowPtr = (int*)calloc(rows + 1, sizeof(int));
    matrix->values = (double*)calloc(nzmax, sizeof(double));
    return matrix;
}

struct darray *create(int rows, int cols)
{
    struct darray* matrix;
    matrix = (struct darray*) calloc(1, sizeof(struct darray));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->array = (double*)calloc(rows*cols, sizeof(double));
    return matrix;
}
