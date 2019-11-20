#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct sparsemat {
	/*
	nzmax = number of non zeros in a sparse matrix
	rows  = number of rows in a sparse matrix
	cols  = number of cols in a sparse matrix
	rowPtr = number of non zero entries in each row of a sparse matrix
	in a cumulative format
	colInd = column number for non zero entries in a sparse matrix
	values = non zero entries in the sparse matrix
	Note: for more details see structure of a sparse matrix in a CSR
	format
	*/
	int nzmax;
	int rows;
	int cols;
	int *rowPtr;
	int *colInd;
	double *values;
};

struct darray
{
	double * array;
	int rows;
	int cols;
};

//Structure of an integer array that also contains information about
//the dimensions of the array
struct iarray
{
	int * array;
	int rows;
	int cols;
};