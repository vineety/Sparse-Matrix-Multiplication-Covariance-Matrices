#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "matrix_def.h"
#include "functions.h"

static void sparseMatrix0(struct sparsemat* const matrixc, const int rows, const int cols)
{
  matrixc->nzmax = 0;
  matrixc->rows = rows;
  matrixc->cols = cols;
  matrixc->rowPtr = (int*)calloc(((size_t)matrixc->rows + 1), sizeof(int));
}

void sparse_sym(const struct sparsemat* const matrixA, const struct sparsemat* const matrixB, struct sparsemat* const matrixC, int imemSize)
{
  int i, j, scratch, counter, threads, subDivision;
  struct iarray rowDistribute = { 0 };
  struct sparsemat constant = { 0, 0, 0, NULL, NULL, NULL };
  struct sparsemat* dimensions = NULL;
  
  if (matrixA->nzmax == 0 || matrixB->nzmax == 0)
  {
    sparseMatrix0(matrixC, matrixA->rows, matrixB->cols);
    return;
  }
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
#else
  threads = 1;
#endif

  subDivision = threads;
  limits(matrixA->rows, subDivision, &rowDistribute);
  subDivision = rowDistribute.rows;
  
  matrixC->rows = matrixA->rows;
  matrixC->cols = matrixB->cols;
  
  if (imemSize > 10 || imemSize < 0)
  {
    imemSize = (int)ceil((double)(matrixA->rows * (size_t)matrixB->cols * 0.1) / threads);
  }
  else
  {
    imemSize = (int)ceil((double)(matrixA->rows * (size_t)matrixB->cols * ((double)imemSize / 100.0)) / threads);
  }
  if (imemSize < matrixB->cols)
  {
    imemSize = matrixB->cols;
  }
  
  dimensions = (struct sparsemat*)calloc(subDivision, sizeof(struct sparsemat));
  if (dimensions == NULL)
  {
    fprintf(stderr, "Memory allocation failed for dimensions\n");
    return;
  }
  
  // Parallel block for computation
#ifdef USE_OPENMP
#pragma omp parallel firstprivate(imemSize) private(i)
#endif
{
#ifdef USE_OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (i = 0; i < subDivision; ++i)
  {
    dimensions[i] = constant;
    
    
    sparsework_sym(matrixA, matrixB, &dimensions[i], rowDistribute.array[i], rowDistribute.array[i + rowDistribute.rows], imemSize);
    
    
    // Update nzmax atomically
#ifdef USE_OPENMP
#pragma omp atomic
#endif
    matrixC->nzmax += dimensions[i].nzmax;
  }
}

// Allocate memory for matrixC based on the total nzmax
matrixC->rowPtr = (int*)calloc(((size_t)matrixC->rows + 1), sizeof(int));
matrixC->colInd = (int*)calloc((size_t)matrixC->nzmax, sizeof(int));
matrixC->values = (double*)calloc((size_t)matrixC->nzmax, sizeof(double));

// Memory allocation failure check
if (matrixC->rowPtr == NULL || matrixC->colInd == NULL || matrixC->values == NULL)
{
  fprintf(stderr, "Memory allocation failed for matrixC\n");
  free(dimensions);
  return;
}

// Initialize variables for final population
scratch = 0;
counter = 1;

for (i = 0; i < subDivision; ++i)
{
  // Update rowPtr for each sub-division
  for (int j = 0; j < rowDistribute.array[i + rowDistribute.rows] - rowDistribute.array[i] + 1; ++j)
  {
    matrixC->rowPtr[counter] = matrixC->rowPtr[counter - 1] + dimensions[i].rowPtr[j];
    counter += 1;
  }
  
  // Update colInd and values for each sub-division
  if (dimensions[i].colInd != NULL && dimensions[i].values != NULL)
  {
    for (j = 0; j < dimensions[i].nzmax; ++j)
    {
      matrixC->colInd[scratch] = dimensions[i].colInd[j];
      matrixC->values[scratch] = dimensions[i].values[j];
      scratch += 1;
    }
  }
  // Clean up the thread-local sparse matrix
  destroy_sparsemat(&dimensions[i]);
}

destroy_iarray(&rowDistribute);
free(dimensions);
dimensions = NULL;
}

void sparse_nosym(const struct sparsemat* const matrixA, const struct sparsemat* const matrixB, struct sparsemat* const matrixC, int imemSize)
{
  int i, j, scratch, counter, threads, subDivision;
  struct iarray rowDistribute = { 0 };
  struct sparsemat constant = { 0, 0, 0, NULL, NULL, NULL };
  struct sparsemat* dimensions = NULL;
  
  if (matrixA->nzmax == 0 || matrixB->nzmax == 0)
  {
    sparseMatrix0(matrixC, matrixA->rows, matrixB->cols);
    return;
  }
  
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
#else
  threads = 1;
#endif
  subDivision = threads;
  limits(matrixA->rows, subDivision,&rowDistribute);
  subDivision = rowDistribute.rows;
  
  matrixC->rows = matrixA->rows;
  matrixC->cols = matrixB->cols;
  
  if (imemSize > 10 || imemSize < 0)
  {
    imemSize = (int)ceil((double)(matrixA->rows * (size_t)matrixB->cols * 0.1) / threads);
  }
  else
  {
    imemSize = (int)ceil((double)(matrixA->rows * (size_t)matrixB->cols * ((double)imemSize / 100.0)) / threads);
  }
  if (imemSize < matrixB->cols)
  {
    imemSize = matrixB->cols;
  }
  
  dimensions = (struct sparsemat*)calloc(subDivision, sizeof(struct sparsemat));
  if (dimensions == NULL)
  {
    fprintf(stderr, "Memory allocation failed for dimensions\n");
    return;
  }
  
  // Parallel block for computation
#ifdef USE_OPENMP
#pragma omp parallel firstprivate(imemSize) private(i)
#endif
{
#ifdef USE_OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (i = 0; i < subDivision; ++i)
  {
    dimensions[i] = constant;
    
    
    sparsework_nosym(matrixA, matrixB, &dimensions[i], rowDistribute.array[i], rowDistribute.array[i + rowDistribute.rows], imemSize);
    
    
    // Update nzmax atomically
#ifdef USE_OPENMP
#pragma omp atomic
#endif
    matrixC->nzmax += dimensions[i].nzmax;
  }
}
// Allocate memory for matrixC based on the total nzmax
matrixC->rowPtr = (int*)calloc(((size_t)matrixC->rows + 1), sizeof(int));
matrixC->colInd = (int*)calloc((size_t)matrixC->nzmax, sizeof(int));
matrixC->values = (double*)calloc((size_t)matrixC->nzmax, sizeof(double));

// Memory allocation failure check
if (matrixC->rowPtr == NULL || matrixC->colInd == NULL || matrixC->values == NULL)
{
  fprintf(stderr, "Memory allocation failed for matrixC\n");
  free(dimensions);
  return;
}

// Initialize variables for final population
scratch = 0;
counter = 1;

for (i = 0; i < subDivision; ++i)
{
  // Update rowPtr for each sub-division
  for (int j = 0; j < rowDistribute.array[i + rowDistribute.rows] - rowDistribute.array[i] + 1; ++j)
  {
    matrixC->rowPtr[counter] = matrixC->rowPtr[counter - 1] + dimensions[i].rowPtr[j];
    counter += 1;
  }
  
  // Update colInd and values for each sub-division
  if (dimensions[i].colInd != NULL && dimensions[i].values != NULL)
  {
    for (j = 0; j < dimensions[i].nzmax; ++j)
    {
      matrixC->colInd[scratch] = dimensions[i].colInd[j];
      matrixC->values[scratch] = dimensions[i].values[j];
      scratch += 1;
    }
  }
  // Clean up the thread-local sparse matrix
  destroy_sparsemat(&dimensions[i]);
}

destroy_iarray(&rowDistribute);
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

void sparse_combined_optimized(const struct sparsemat* const matrixA, const struct sparsemat* const matrixB, struct sparsemat* const matrixC, int imemSize, bool isSymmetric)
{
	int i, j, scratch, counter, threads, subDivision;
	struct iarray rowDistribute = { 0 };
	struct sparsemat constant = { 0, 0, 0, NULL, NULL, NULL };
	struct sparsemat* dimensions = NULL;

	if (matrixA->nzmax == 0 || matrixB->nzmax == 0)
	{
		sparseMatrix0(matrixC, matrixA->rows, matrixB->cols);
		return;
	}

	threads = omp_get_max_threads();
	subDivision = threads;
	rowDistribute = limits(matrixA->rows, subDivision);
	subDivision = rowDistribute.rows;

	matrixC->rows = matrixA->rows;
	matrixC->cols = matrixB->cols;

	if (imemSize > 10 || imemSize < 0)
	{
		imemSize = (int)ceil((double)(matrixA->rows * (size_t)matrixB->cols * 0.1) / threads);
	}
	else
	{
		imemSize = (int)ceil((double)(matrixA->rows * (size_t)matrixB->cols * ((double)imemSize / 100.0)) / threads);
	}
	if (imemSize < matrixB->cols)
	{
		imemSize = matrixB->cols;
	}

	dimensions = (struct sparsemat*)calloc(subDivision, sizeof(struct sparsemat));
	if (dimensions == NULL)
	{
		fprintf(stderr, "Memory allocation failed for dimensions\n");
		return;
	}

	// Parallel block for computation
#pragma omp parallel firstprivate(imemSize) private(i)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < subDivision; ++i)
		{
			dimensions[i] = constant;

			if (isSymmetric)
			{
				sparsework_sym(matrixA, matrixB, &dimensions[i], rowDistribute.array[i], rowDistribute.array[i + rowDistribute.rows], imemSize);
			}
			else
			{
				sparsework_nosym(matrixA, matrixB, &dimensions[i], rowDistribute.array[i], rowDistribute.array[i + rowDistribute.rows], imemSize);
			}

			// Update nzmax atomically
#pragma omp atomic
			matrixC->nzmax += dimensions[i].nzmax;
		}
	}

	// Allocate memory for matrixC based on the total nzmax
	matrixC->rowPtr = (int*)calloc(((size_t)matrixC->rows + 1), sizeof(int));
	matrixC->colInd = (int*)calloc((size_t)matrixC->nzmax, sizeof(int));
	matrixC->values = (double*)calloc((size_t)matrixC->nzmax, sizeof(double));

	// Memory allocation failure check
	if (matrixC->rowPtr == NULL || matrixC->colInd == NULL || matrixC->values == NULL)
	{
		fprintf(stderr, "Memory allocation failed for matrixC\n");
		free(dimensions);
		return;
	}

	// Initialize variables for final population
	scratch = 0;
	counter = 1;

	for (i = 0; i < subDivision; ++i)
	{
		// Update rowPtr for each sub-division
		for (int j = 0; j < rowDistribute.array[i + rowDistribute.rows] - rowDistribute.array[i] + 1; ++j)
		{
			matrixC->rowPtr[counter] = matrixC->rowPtr[counter - 1] + dimensions[i].rowPtr[j];
			counter += 1;
		}

		// Update colInd and values for each sub-division
		if (dimensions[i].colInd != NULL && dimensions[i].values != NULL)
		{
			for (j = 0; j < dimensions[i].nzmax; ++j)
			{
				matrixC->colInd[scratch] = dimensions[i].colInd[j];
				matrixC->values[scratch] = dimensions[i].values[j];
				scratch += 1;
			}
		}
		// Clean up the thread-local sparse matrix
		destroy(&dimensions[i]);
	}

	destroy(&rowDistribute);
	free(dimensions);
	dimensions = NULL;
}
