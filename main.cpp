#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "functions.h"
#include "matrix_def.h"


int main(void)
{

    /* Test parameters for checking code
    1 = sparse_nosym
    2 = sparse_sym
    3 = dense_nosym
    4 = dense_sym
    */
    int test = 0;


    printf("Four functions in the manuscript can be tested\n\n");
    printf("The integer codes for testing these functions are: \n");
    printf("1 = sparse_nosym, 2 = sparse_sym, 3 = dense_nosym, 4 = dense_sym \n\n");
    printf("Now enter a test code\n");
    if (scanf("%d", &test)>0)
    {
        printf("You entered: %d\n", test);
    }
    else
    {
      printf("You did not enter any number.\n");
    }


     struct sparsemat *matrixa = 0, *matrixb = 0;
	 /* These matrices are used to show the structure of matrices A and B
	 and are not used for any computations*/

	 struct darray *matrixa_p = 0, *matrixb_p = 0;

	 int nrows_a = 4;
	 int ncols_a = 4;
	 int nzmax_a = 5;

	 int nrows_b = 4;
	 int ncols_b = 4;
	 int nzmax_b = 4;

	 matrixa_p = create(nrows_a, ncols_a);
	 matrixb_p = create(nrows_b, ncols_b);

	 /* Fill matrices */
	 matrixa_p->rows = nrows_a;
	 matrixa_p->cols = ncols_a;

	 matrixa_p->array[0] = 0.64;
	 matrixa_p->array[3] = 0.72;
	 matrixa_p->array[5] = 0.67;
	 matrixa_p->array[9] = 0.32;
	 matrixa_p->array[12] = 0.10;

	 /*Print Input Matrices*/
	 printf("\nNow Printing Matrix A as 2d Array\n\n");
	 printMatForm(matrixa_p->array, matrixa_p->rows, matrixa_p->cols, 5, 4);
	 destroy(matrixa_p);
	 free(matrixa_p);

	 /* Fill matrices */
	 matrixb_p->rows = nrows_b;
	 matrixb_p->cols = ncols_b;

	 matrixb_p->array[0] = 0.72;
	 matrixb_p->array[5] = 0.72;
	 matrixb_p->array[10] = 0.72;
	 matrixb_p->array[15] = 0.72;

	 /*Print Input Matrices*/
	 printf("\nNow Printing Matrix B as 2d Array\n\n");
	 printMatForm(matrixb_p->array, matrixb_p->rows, matrixb_p->cols, 5, 4);
	 destroy(matrixb_p);
	 free(matrixb_p);

    matrixa = create(nrows_a, ncols_a, nzmax_a);
    matrixb = create(nrows_b, ncols_b, nzmax_b);

	 /* Create CSR Matrices A and B*/

    matrixa->rowPtr[0] = 0;
    matrixa->rowPtr[1] = 2;
    matrixa->rowPtr[2] = 3;
    matrixa->rowPtr[3] = 4;
    matrixa->rowPtr[4] = 5;

    matrixa->colInd[0] = 0;
    matrixa->colInd[1] = 3;
    matrixa->colInd[2] = 1;
    matrixa->colInd[3] = 1;
    matrixa->colInd[4] = 0;

    matrixa->values[0] = 0.64;
    matrixa->values[1] = 0.72;
    matrixa->values[2] = 0.67;
    matrixa->values[3] = 0.32;
    matrixa->values[4] = 0.10;

    matrixb->rowPtr[0] = 0;
    matrixb->rowPtr[1] = 1;
    matrixb->rowPtr[2] = 2;
    matrixb->rowPtr[3] = 3;
    matrixb->rowPtr[4] = 4;

    matrixb->colInd[0] = 0;
    matrixb->colInd[1] = 1;
    matrixb->colInd[2] = 2;
    matrixb->colInd[3] = 3;

    matrixb->values[0] = 0.72;
    matrixb->values[1] = 0.72;
    matrixb->values[2] = 0.72;
    matrixb->values[3] = 0.72;


    if (test == 1)
        {
            struct sparsemat *matrixc = 0;
            matrixc = (struct sparsemat*) calloc(1, sizeof(struct sparsemat));

            sparse_nosym(matrixa, matrixb, matrixc,5);

            printf("\nInvoked Function sparse_nosym: output would be a full matrix in CSR form\n\n");

            printf("Rows in Matrix C = %i\n\n", matrixc->rows);
            printf("Columns in Matrix C = %i\n\n", matrixc->cols);
            printf("Total Non-Zero Entries in Matrix C = %i\n\n", matrixc->nzmax);

            printf("Now Printing Row Pointers of Matrix C\n\n");
            print_array_1d(matrixc->rowPtr, matrixc->rows + 1);
            printf("Now Printing Column Indexes of Matrix C\n\n");
            print_array_1d(matrixc->colInd, matrixc->nzmax);
            printf("Now Printing Values of Matrix C\n\n");
            print_array_1d(matrixc->values, matrixc->nzmax);
            destroy(matrixc);
            free(matrixc);
        }
    else if (test == 2)
        {
            struct sparsemat *matrixc = 0;
            matrixc = (struct sparsemat*) calloc(1, sizeof(struct sparsemat));

            sparse_sym(matrixa, matrixb, matrixc, 5);

            printf("\nInvoked Function sparse_sym: output would be a upper triangular matrix in CSR form\n\n");

            printf("Rows in Matrix C = %i\n\n", matrixc->rows);
            printf("Columns in Matrix C = %i\n\n", matrixc->cols);
            printf("Total Non-Zero Entries in Matrix C = %i\n\n", matrixc->nzmax);

            printf("Now Printing Row Pointers of Matrix C\n\n");
            print_array_1d(matrixc->rowPtr, matrixc->rows + 1);
            printf("Now Printing Column Indexes of Matrix C\n\n");
            print_array_1d(matrixc->colInd, matrixc->nzmax);
            printf("Now Printing Values of Matrix C\n\n");
            print_array_1d(matrixc->values, matrixc->nzmax);
            destroy(matrixc);
            free(matrixc);
        }
    else if (test == 3)
        {
            struct darray *matrixc = 0;
            matrixc = (struct darray*) calloc(1, sizeof(struct darray));


            dense_nosym(matrixa, matrixb, matrixc);

            printf("\nInvoked Function dense_nosym: output would be a dense matrix\n\n");

            printf("Rows in Matrix C = %i\n\n", matrixc->rows);
            printf("Columns in Matrix C = %i\n\n", matrixc->cols);

            printf("Now Printing Matrix C as 2d Array\n\n");
            printMatForm(matrixc->array, matrixc->rows, matrixc->cols, 5, 4);
            destroy(matrixc);
            free(matrixc);
        }
    else if (test == 4)
        {
            struct darray *matrixc = 0;
            matrixc = (struct darray*) calloc(1, sizeof(struct darray));

            dense_sym(matrixa, matrixb, matrixc);

            printf("\nInvoked Function dense_nosym: output would be a dense upper triangular matrix\n\n");

            printf("Rows in Matrix C = %i\n\n", matrixc->rows);
            printf("Columns in Matrix C = %i\n\n", matrixc->cols);

            printf("Now Printing Matrix C as 2d Array\n\n");
            printMatForm(matrixc->array, matrixc->rows, matrixc->cols, 5, 4);
            destroy(matrixc);
            free(matrixc);
        }
    else
        {
            printf("Not a Valid Test Value\n");
        }
    destroy(matrixa);
    free(matrixa);
    destroy(matrixb);
    free(matrixb);
    printf("Press Enter To Exit\n");
    getchar();
    getchar();
    return 0;
}
