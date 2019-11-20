#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "matrix_def.h"
#include "functions.h"

static void sparseMatrixOne(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct sparsemat * const matrixc)

/**************************************************************************
Written by: Vineet Yadav
Date: 12/15/2014

This function computes C = A*B where A and B are in CSR form. C is returned in
a CSR form. The function handles the special case where A and/or B are (1,1) in
dimensions and non-zeros(nzmax) in A and B > 0. In this special case a scalar *
matrix multiplication is performed without resorting to FG Gustavson's Algorithm
for performing sparse matrix multiplication.

Input:

(1) matrixa in CSR sparse format
(2) matrixb in CSR format
(3) matrixc is a dynamically allocated CSR structure whose colInd, values and rowPtr
members are dynamically allocated by calloc function at runtime by this function.

Output:
(1) On output matrixc is updated and will contain result from the matrix multiplication
of A*B.
 ****************************************************************************/
{
    int i;
    // Case where A is a (1,1) matrix
    if (matrixa->rows == 1 && matrixa->cols == 1)
        {

            matrixc->nzmax = matrixb->nzmax;
            matrixc->rows = matrixb->rows;
            matrixc->cols = matrixb->cols;

            // Memory allocation
            matrixc->rowPtr = (int *) calloc(matrixc->rows + 1, sizeof (int));
            matrixc->colInd = (int *) calloc(matrixc->nzmax, sizeof (int));
            matrixc->values = (double *) calloc(matrixc->nzmax, sizeof (double));

            // compute row pointers and column index and values simultaneously
            #pragma omp parallel
            {
                #pragma omp for nowait
                for (i = 0; i < matrixc->rows + 1; ++i)
                    matrixc->rowPtr[i] = matrixb->rowPtr[i];
                #pragma omp for nowait
                for (i = 0; i < matrixc->nzmax; ++i)
                    {
                        matrixc->colInd[i] = matrixb->colInd[i];
                        matrixc->values[i] = matrixa->values[0] * matrixb->values[i];
                    }
            }

        }
    // Case where B is a (1,1) matrix
    else if (matrixb->rows == 1 && matrixb->cols == 1)
        {
            matrixc->nzmax = matrixa->nzmax;
            matrixc->rows = matrixa->rows;
            matrixc->cols = matrixa->cols;

            // Memory allocation

            matrixc->rowPtr = (int *) calloc(matrixc->rows + 1, sizeof (int));
            matrixc->colInd = (int *) calloc(matrixc->nzmax, sizeof (int));
            matrixc->values = (double *) calloc(matrixc->nzmax, sizeof (double));

            // compute row pointers and column index and values simultaneously
            #pragma omp parallel
            {
                #pragma omp for nowait
                for (i = 0; i < matrixc->rows + 1; ++i)
                    matrixc->rowPtr[i] = matrixa->rowPtr[i];
                #pragma omp for nowait
                for (i = 0; i < matrixc->nzmax; ++i)
                    {
                        matrixc->colInd[i] = matrixa->colInd[i];
                        matrixc->values[i] = matrixb->values[0] * matrixa->values[i];
                    }
            }
        }
}

static void sparseMatrix0(struct sparsemat * const matrixc, const int rows, const int cols)
{
    /**************************************************************************
    Written by: Vineet Yadav
    Date: 12/3/2014

    This function computes C = A*B where A and B are in CSR form. C is returned in
    a CSR form. The function handles the special case where non-zeros(nzmax) in A
    and/or B == 0. In this special case no matrix multiplication is performed.

    Input:

    (1) matrixa in CSR sparse format
    (2) matrixb in CSR format
    (3) matrixc is a dynamically allocated CSR structure whose rowPtr pointer is
    dynamically allocated by calloc function at runtime by this function.

    Output:
    (1) On output matrixc is updated and will contain result from the matrix multiplication
    of A*B.
     ****************************************************************************/


    matrixc->nzmax = 0;
    matrixc->rows = rows;
    matrixc->cols = cols;
    //Allocate memory for row pointers and initialize with zeros
    matrixc->rowPtr = (int *) calloc(matrixc->rows + 1, sizeof (int));
}

void sparse_nosym(const struct sparsemat * const  matrixA, const struct sparsemat
                  * const  matrixB, struct sparsemat *  const matrixC, int imemSize)
{
    /*
    This function computes sparse X sparse matrix multiplication for two
    conformable sparse matrices stored in compressed sparse row (CSR) format with
    zero based indices.

    Mathematically perform: C = A*B where a and b are given in sparse CSR format
    and C is obtained in CSR form. The method here is based on FG Gustavsons algorithm
    for performing sparse matrix multiplication given in paper:
    Two fast algorithms for sparse matrices: Multiplication and permuted transposition
    (published in ACM in 1978).
    This function modifies this algorithm to do parallel sparse X sparse
    matrix multiplication. It uses OpenMP for performing matrix multiplication.

    Steps:
    (1)  Serial: Use a load balancer to sub-divide work by rows for number of threads
    specified by the user.
    (2)  Parallel: Use OpenMp and or MPI to compute column indices and matrix entries for the
    output matrix.
    (3)  Serial: Collect results in a structure.

    Inputs (Unmodified and Modified by the Function): All inputs are required

    (1)   matrixA in CSR format.
    (2)   matrixB in CSR format.
    (3)   matrixC a dynamically allocated CSR structure. On return it will have results of
    (4)   imemSize size of initially allocated memory specified at a percentage of nnz in the output matrix.
	 For e.g. if the the output matrix is (4, 5) matrix then if the the user specify imemsize as 5 then
	 out of 20 entries in the output matrix 5% would are assumed to be non-zeros (i.e. ceil(0.05*20)).

     If a thread runs out of memory then each thread would increase the memory by doubling
    the previous allocation by adding space for more entries.
    */
    int i, j, scratch, counter, threads, subDivision;
    /*rowDistribute determines the subDivision of work i.e. the amount of work
    or rows of distance  matrix processed by each thread*/
    struct iarray rowDistribute = { 0 };
    struct sparsemat constant = { 0, 0, 0, NULL, NULL, NULL };
    /*The structure in which final results from the work done by each thread
    would be collected. Array of structure that would have results computed by each thread.*/
    struct sparsemat *dimensions = NULL;

    /* handles special case when non zeros in matrix A or matrix B is 0*/
    if (matrixA->nzmax == 0 || matrixB->nzmax == 0)
        {
            sparseMatrix0(matrixC, matrixA->rows, matrixB->cols);
            return;
        }
    /*handles special case when non zeros in matrix A or matrix B is (1,1)*/
    if ((matrixA->rows == 1 && matrixA->cols == 1) || (matrixB->rows == 1 && matrixB->cols == 1))
        {
            sparseMatrixOne(matrixA, matrixB, matrixC);
            return;
        }
    /*Threads are set to max number of threads*/
    threads = omp_get_max_threads();
    /*Subdivision is also set to omp_get_max_threads but can be higher or lower
    depending on desired load balancing*/
    subDivision = threads;
    /*Number of rows to be processed by each thread. Distribution of the work is
    determined on the basis of number of threads and variable subDivision*/
    rowDistribute = limits(matrixA->rows, subDivision);
    /*Update subDivision based on number of Processors*/
    subDivision = rowDistribute.rows;
    /*Assign dimensions of matrixC*/
    matrixC->rows = matrixA->rows;
    matrixC->cols = matrixB->cols;

    /*if memIncrese is set to be > than 10% than default sparsity
    is assumed (O(1) estimate). Default sparsity of 10% equally divided
    across rows.*/

    if (imemSize>10 || imemSize<0)
        {
            imemSize = ceil((matrixA->rows*matrixB->cols*0.1) / threads);
        }
    else
        {
            imemSize = ceil((matrixA->rows*matrixB->cols*((double)imemSize / (double) 100)) / threads);
        }
    if (imemSize < matrixB->cols)
        {
            imemSize = matrixB->cols;
        }
    /*Each thread would write result to a unique structure because of which
    dimensions is equal to subDivision*/
    dimensions = (struct sparsemat *) malloc(subDivision * sizeof(struct sparsemat));

    #pragma omp parallel firstprivate (imemSize)
    {
        #pragma omp for
        /*Total number of subDivisions. Each subDivision is uniquely
        given to each thread and the results of its computation is stored in the
        array of structures (dimension)*/
        for (i = 0; i < subDivision; ++i)
            {
                /*Initialize each structure of the array with 0 or NULL*/
                dimensions[i] = constant;
                /*Do matrix multiplication for subdivided rows*/
                sparsework_nosym(matrixA, matrixB, &dimensions[i],
                                 rowDistribute.array[i], rowDistribute.array[i + rowDistribute.rows], imemSize);
                /*compute non-zero entries in each structure*/
                matrixC->nzmax += dimensions[i].nzmax;
            }
        /*collect results only done by master thread*/
        #pragma omp single
        {
            destroy(&rowDistribute);
            matrixC->rowPtr = (int *)calloc(matrixC->rows + 1, sizeof(int));
            matrixC->colInd = (int *)calloc(matrixC->nzmax, sizeof(int));
            matrixC->values = (double *)calloc(matrixC->nzmax, sizeof(double));
            /*Counters to increase for storing results in the final Array*/
            scratch = 0;
            counter = 1;
            for (i = 0; i < subDivision; ++i)
                {
                    if (&dimensions[i] != NULL)
                        {
                            for (j = 0; j < dimensions[i].nzmax; ++j)
                                {
                                    matrixC->colInd[scratch] = dimensions[i].colInd[j];
                                    matrixC->values[scratch] = dimensions[i].values[j];
                                    scratch += 1;
                                }
                            /*Special condition for accounting for case when a row has 0 entries*/
                            if (dimensions[i].rows == 0)
                                {
                                    matrixC->rowPtr[counter] = matrixC->rowPtr[counter - 1];
                                    counter += 1;
                                }
                            /*General Case*/
                            for (j = 0; j < dimensions[i].rows; ++j)
                                {
                                    matrixC->rowPtr[counter] = matrixC->rowPtr[counter - 1] + dimensions[i].rowPtr[j];
                                    counter += 1;
                                }
                        }
                    /*delete allocated arrays after exiting from the parallel region*/
                    destroy(&dimensions[i]);
                }
        }
    }
    /*delete allocated array of structures*/
    free(dimensions);
    dimensions = NULL;
}
void sparse_sym(const struct sparsemat * const  matrixA, const struct sparsemat
                * const  matrixB, struct sparsemat *  const matrixC, int imemSize)
{
    /*
    This function computes sparse X sparse matrix multiplication for two
    conformable sparse matrices stored in compressed sparse row (CSR) format with
    zero based indices and return the upper triengular portion of the symmetric output
	 matrix. It is important that the output matrix that results from the multiplication of two
	 sparse matrices results in a symmetric output matrix.

    Mathematically perform: C = A*B where a and b are given in sparse CSR format
    and C is obtained in CSR form. The method here is based on FG Gustavsons algorithm
    for performing sparse matrix multiplication given in paper:
    Two fast algorithms for sparse matrices: Multiplication and permuted transposition
    (published in ACM in 1978).
    This function modifies this algorithm to do parallel sparse X sparse
    matrix multiplication. It uses OpenMP for performing matrix multiplication.

    Steps:
    (1)  Serial: Use a load balancer to sub-divide work by rows for number of threads
    specified by the user.
    (2)  Parallel: Use OpenMp and or MPI to compute column indices and matrix entries for the
    output matrix.
    (3)  Serial: Collect results in a structure.

    Inputs (Unmodified and Modified by the Function): All inputs are required

    (1)   matrixA in CSR format.
    (2)   matrixB in CSR format.
    (3)   matrixC a dynamically allocated CSR structure. On return it will have results of
    (4)   imemSize size of initially allocated memory specified at a percentage of nnz in the outmut matrix.
	 For e.g. if the the output matrix would be (4, 5) matrix then if the the user specify imemsize is 5 then
	 out of 20 entries in the output matrix 5% would be non-zeros (i.e. ceil(0.05*20))

    If a thread runs out of memory then each thread would increase the memory by doubling
    the previous allocation by adding space for more entries.
    */
    int i, j, scratch, counter, threads, subDivision;
    /*rowDistribute determines the subDivision of work i.e. the amount of work
    or rows of distance  matrix processed by each thread*/
    struct iarray rowDistribute = { 0 };
    struct sparsemat constant = { 0, 0, 0, NULL, NULL, NULL };
    /*The structure in which final results from the work done by each thread
    would be collected. Array of structure that would have results computed by each thread.*/
    struct sparsemat *dimensions = NULL;
    if (matrixA->nzmax == 0 || matrixB->nzmax == 0)
        {
            sparseMatrix0(matrixC, matrixA->rows, matrixB->cols);
            return;
        }
    /*Threads are set to max number of threads*/
    threads = omp_get_max_threads();
    /*Subdivision is also set to omp_get_max_threads but can be higher or lower
    depending on desired load balancing*/
    subDivision = threads;
    /*Number of rows to be processed by each thread. Distribution of the work is
    determined on the basis of number of threads and variable subDivision*/
    rowDistribute = limits(matrixA->rows, subDivision);
    /*Update subDivision based on number of Processors*/
    subDivision = rowDistribute.rows;
    /*Assign dimensions of matrixC*/
    matrixC->rows = matrixA->rows;
    matrixC->cols = matrixB->cols;

    /*if memIncrese is set to be > than 10% than default sparsity
    is assumed (O(1) estimate). Default sparsity of 10% equally divided
    across rows.*/
    if (imemSize>10 || imemSize<0)
        {
            imemSize = ceil((matrixA->rows*matrixB->cols*0.1) / threads);
        }
    else
        {
            imemSize = ceil((matrixA->rows*matrixB->cols*((double)imemSize / (double) 100)) / threads);
        }
    if (imemSize < matrixB->cols)
        {
            imemSize = matrixB->cols;
        }
    /*Each thread would write result to a unique structure because of which
    dimensions is equal to subDivision*/
    dimensions = (struct sparsemat *) malloc(subDivision * sizeof(struct sparsemat));

    #pragma omp parallel firstprivate (imemSize)
    {
        #pragma omp for
        /*Total number of subDivisions. Each subDivision is uniquely
        given to each thread and the results of its computation is stored in the
        array of structures (dimension)*/
        for (i = 0; i < subDivision; ++i)
            {
                /*Initialize each structure of the array with 0 or NULL*/
                dimensions[i] = constant;
                /*Do matrix multiplication for subdivided rows*/
                sparsework_sym(matrixA, matrixB, &dimensions[i],
                               rowDistribute.array[i], rowDistribute.array[i + rowDistribute.rows], imemSize);
                /*compute non-zero entries in each structure*/
                matrixC->nzmax += dimensions[i].nzmax;
            }
        /*collect results only done by master thread*/
        #pragma omp single
        {
            destroy(&rowDistribute);
            matrixC->rowPtr = (int *)calloc(matrixC->rows + 1, sizeof(int));
            matrixC->colInd = (int *)calloc(matrixC->nzmax, sizeof(int));
            matrixC->values = (double *)calloc(matrixC->nzmax, sizeof(double));
            /*Counters to increase for storing results in the final Array*/
            scratch = 0;
            counter = 1;
            for (i = 0; i < subDivision; ++i)
                {
                    if (&dimensions[i] != NULL)
                        {
                            for (j = 0; j < dimensions[i].nzmax; ++j)
                                {
                                    matrixC->colInd[scratch] = dimensions[i].colInd[j];
                                    matrixC->values[scratch] = dimensions[i].values[j];
                                    scratch += 1;
                                }
                            /*Special condition for accounting for case when a row has 0 entries*/
                            if (dimensions[i].rows == 0)
                                {
                                    matrixC->rowPtr[counter] = matrixC->rowPtr[counter - 1];
                                    counter += 1;
                                }
                            /*General Case*/
                            for (j = 0; j < dimensions[i].rows; ++j)
                                {
                                    matrixC->rowPtr[counter]
                                        = matrixC->rowPtr[counter - 1] + dimensions[i].rowPtr[j];
                                    counter += 1;
                                }
                        }
                    /*delete allocated arrays after exiting from the parallel region*/
                    destroy(&dimensions[i]);
                }
        }
    }
    /*delete allocated array of structures*/
    free(dimensions);
    dimensions = NULL;
}

void sparsemult_single(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct sparsemat * matrixc, int imemsize)
/*
This serial function computes sparse X sparse matrix multiplication for two
conformable sparse matrices stored in compressed sparse row (CSR) format with
zero based indices and return results in matrixc.


Mathematically perform: C = A*B where A and B are given in sparse CSR format
and C is obtained in CSR form. The method here is based on FG Gustavsons algorithm
for performing sparse matrix multiplication given in paper:
Two fast algorithms for sparse matrices: Multiplication and permuted transposition
(published in ACM in 1978).

Inputs (Unmodified and Modified by the Function): All inputs are required

(1)   matrixA in CSR format.
(2)   matrixB in CSR format.
(3)   matrixC a dynamically allocated CSR structure. On return it will have results of matrix multiplication
(4)   imemSize size of initially allocated memory specified at a percentage of nnz in the outmut matrix.
For e.g. if the the output matrix would be (4, 5) matrix then if the the user specify imemsize is 5 then
out of 20 entries in the output matrix 5% would be non-zeros (i.e. ceil(0.05*20)). We have allowed user to specify
this as a % however this can also be done internally without the need to specify it earlier. If more memory
is required then the routine doubles the size of memory*/


/*loop counters and scratch variables*/
{

	int i, j, k, nzrows = 0;
	int col_num_a;

	/* Output Matrix */

	int *workarray = NULL;
	double value = 0;

	if (imemsize>10 || imemsize<0)
	{
		imemsize = ceil(matrixa->rows*matrixb->cols*0.1);
	}
	else
	{
		imemsize = ceil((matrixa->rows*matrixb->cols*((double)imemsize / (double)100)));
	}
	if (imemsize < matrixb->cols)
	{
		imemsize = matrixb->cols;
	}

	const int increasefactor = 2;

	/* Workarray that contains indicator of non zero entries*/
	/* Initilize different fields of output sparse matrix without a function */

	matrixc->rows = matrixa->rows;
	matrixc->cols = matrixb->cols;
	matrixc->nzmax = 0;
	matrixc->rowPtr = (int*)calloc(matrixa->rows + 1, sizeof(int));

	matrixc->colInd = (int*)calloc(imemsize, sizeof(int));
	matrixc->values = (double*)calloc(imemsize, sizeof(double));

	workarray = (int*)calloc(matrixb->cols, sizeof(int));
	//memInitialize(workarray, matrixb->cols, -1);
	/**************Parallel Part *******************************************/
	for (i = 0; i<matrixa->rows; ++i)
	{
		for (j = matrixa->rowPtr[i]; j <matrixa->rowPtr[i + 1]; ++j)
		{
			value = matrixa->values[j];
			col_num_a = matrixa->colInd[j];
			for (k = matrixb->rowPtr[col_num_a]; k<matrixb->rowPtr[col_num_a + 1]; ++k)
			{
				/*This would give coloumn no of the entry in the B matrix*/

				if (workarray[matrixb->colInd[k]] == 0)
				{
					matrixc->colInd[nzrows] = matrixb->colInd[k];
					workarray[matrixb->colInd[k]] = nzrows;
					matrixc->values[nzrows] = value*matrixb->values[k];
					nzrows += 1;
				}
				else
				{
					matrixc->values[workarray[matrixb->colInd[k]]] += value * matrixb->values[k];
				}
			}
		}
		matrixc->rowPtr[i + 1] = nzrows;
		for (k = matrixc->rowPtr[i]; k<nzrows; ++k)
		{
			workarray[matrixc->colInd[k]] = 0;
		}
		/* Check if more memory is required*/
		if ((imemsize - nzrows)<matrixb->cols)
		{
			modifyalloc(matrixc, increasefactor*imemsize);
			imemsize = increasefactor * imemsize;
		}
	}
	matrixc->nzmax = matrixc->rowPtr[matrixc->rows];
	modifyalloc(matrixc, nzrows);
	free(workarray);
	workarray = NULL;
}
