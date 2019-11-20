#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "matrix_def.h"
#include "functions.h"

void dense_sym(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct darray * const matrixc)
/* This routine multiplies two sparse matrices in CSR format
with zero based indexing and produces a dense upper triangular
matrix based on the assumption that multiplication of two sparse
matrices results in a symmetric matrix*/
{
     int i, j, k, col_num_a, col_num_b, konstant;
     double value;
	 matrixc->rows = matrixa->rows;
	 matrixc->cols = matrixb ->cols;
	 matrixc->array = (double *)calloc(matrixa->rows*matrixb->cols, sizeof(double));
# pragma omp parallel for private (i,j,k,value,col_num_a,col_num_b,konstant)
		 for (i = 0; i < matrixa->rows; i++)
		 {
			 konstant = i*matrixa->rows;
			 for (j = matrixa->rowPtr[i]; j <= matrixa->rowPtr[i + 1] - 1; j++)
			 {
				 value = matrixa->values[j];
				 col_num_a = matrixa->colInd[j];
				 for (k = matrixb->rowPtr[col_num_a];
					 k <= matrixb->rowPtr[col_num_a + 1] - 1; k++)
				 {
					 col_num_b = matrixb->colInd[k];
					 if (i <= col_num_b)
						 matrixc->array[konstant + col_num_b] += value * matrixb->values[k];
				 }
			 }

		 }
}

void dense_nosym(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct darray * const matrixc)

/* This routine multiplies two sparse matrices in CSR format
with zero based indexing and produces a dense upper triangular
matrix based on the assumption that multiplication of two sparse
matrices results in a symmetric matrix*/
{
    int i, j, k, col_num_a, col_num_b, konstant;
    double value;
    matrixc->rows = matrixa->rows;
    matrixc->cols = matrixb->cols;
     matrixc->array = (double *)calloc(matrixa->rows*matrixb->cols, sizeof(double));
    # pragma omp parallel for private (i,j,k,value,col_num_a,col_num_b,konstant)
        for (i = 0; i < matrixa->rows; i++)
            {
                konstant = i*matrixa->rows;
                for (j = matrixa->rowPtr[i]; j <= matrixa->rowPtr[i + 1] - 1; j++)
                    {
                        value = matrixa->values[j];
                        col_num_a = matrixa->colInd[j];
                        for (k = matrixb->rowPtr[col_num_a];
                                k <= matrixb->rowPtr[col_num_a + 1] - 1; k++)
                            {
                                col_num_b = matrixb->colInd[k];
                                matrixc->array[konstant + col_num_b] += value * matrixb->values[k];
                            }
                    }

            }
}
