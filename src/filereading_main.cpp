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
    /* Binary Stream Files that store matrixa and matrixb can be read by below
    given code to produce an output matrix
    Binary file should be written as:

    1. An integer telling rows in a matrix,
    2. An integer that telling number of columns in a matrix
    3. An integer that telling number of non-zero entries in a matrix
    4. Followed by Row Pointers (integer), Col Indices (integer) and Values (double)*/


    struct sparsemat *matrixa = 0, *matrixb = 0, *matrixc;
    int nrows_a = 0;
    int ncols_a = 0;
    int nzmax_a = 0;

    int nrows_b = 0;
    int ncols_b = 0;
    int nzmax_b = 0;

    double start, finish;

    FILE *fp = NULL;

    /* File names*/
    /* Input file names*/
    /* These are example file names */
    const char *filenamea = "/media/sf_D_DRIVE/sparse/Szero.bin";
    const char *filenameb = "/media/sf_D_DRIVE/sparse/SD_zero.bin";
    /* Output file name*/
    const char *filenameo = "/media/sf_D_DRIVE/sparse/output.bin";

    /* READ file A*/
    fp = fopen(filenamea, "rb");

    /*Read dimensions of A*/
    fread(&nrows_a, sizeof (int), 1, fp);
    fread(&ncols_a, sizeof (int), 1, fp);
    fread(&nzmax_a, sizeof (int), 1, fp);

    /* Assign space for storing rowpointers, colindices and values of matrix a*/
    matrixa = create(nrows_a, ncols_a, nzmax_a);

    /* Read data*/
    fread(matrixa->rowPtr, sizeof (int), nrows_a + 1, fp);
    fread(matrixa->colInd, sizeof (int), nzmax_a, fp);
    fread(matrixa->values, sizeof (double), nzmax_a, fp);

    /* Print information about matrix A*/
    printf("Rows in A = %i \n", nrows_a);
    printf("Cols in A = %i \n", ncols_a);
    printf("Non Zeros in A = %i  \n\n", nzmax_a);

    fclose(fp);
    fp = 0;

    /* Read file A*/
    fp = fopen(filenameb, "rb");

    /*Read dimensions of B*/
    fread(&nrows_b, sizeof (int), 1, fp);
    fread(&ncols_b, sizeof (int), 1, fp);
    fread(&nzmax_b, sizeof (int), 1, fp);

    /* Assign space for storing rowpointers, colindices and values of matrix B*/
    matrixb = create(nrows_b, ncols_b, nzmax_b);

    /* Read data*/
    fread(matrixb->rowPtr, sizeof (int), nrows_b + 1, fp);
    fread(matrixb->colInd, sizeof (int), nzmax_b, fp);
    fread(matrixb->values, sizeof (double), nzmax_b, fp);

    /* Print information about matrix B*/
    printf("Rows in B = %i \n", nrows_b);
    printf("Cols in B = %i \n", ncols_b);
    printf("Non Zeros in B = %i  \n\n", nzmax_b);

    fclose(fp);
    fp = 0;

    /* Assign space for storing rowpointers, colindices and values of matrix B*/
    matrixc = (struct sparsemat*) calloc(1, sizeof (struct sparsemat));

    /* Start timer*/
    start = omp_get_wtime();

    /* Compute matrix product*/
    sparse_nosym(matrixa, matrixb, matrixc, 5);

    /* End timer*/
    finish = omp_get_wtime();

    /* Print information about matrix C*/

    printf("Rows in C = %i \n", matrixc->rows);
    printf("Cols in C = %i \n", matrixc->cols);
    printf("Non Zeros in C = %i  \n", matrixc->nzmax);

    /* Print time taken in matrix multiplication*/
    printf("start = %.16g\nend = %.16g\n Time Taken in Seconds = %.16g\n", start, finish, finish - start);

    /* Write output file*/

    fp = fopen(filenameo, "wb");

    fwrite(&matrixc->rows, sizeof(int), 1, fp);
    fwrite(&matrixc->cols, sizeof(int), 1, fp);
    fwrite(&matrixc->nzmax, sizeof(int), 1, fp);
    fwrite(matrixc->rowPtr, sizeof(int), matrixc->rows+1, fp);
    fwrite(matrixc->colInd, sizeof(int), matrixc->nzmax, fp);
    fwrite(matrixc->values, sizeof(double), matrixc->nzmax, fp);
    fclose(fp);
    fp = 0;

    /* Free memory*/
    destroy(matrixa);
    destroy(matrixb);
    destroy(matrixc);
    free(matrixa);
    free(matrixb);
    free(matrixc);

    printf("Press Enter To Exit\n");
    getchar();
    return 0;

}
