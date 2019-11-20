#include <stdio.h>
#include <stdlib.h>

void print_array_1d(double * __restrict Array, int length)
/*------------------------------------------------------------
Purpose: Print arrays
Inputs:
(1) double *Array : A pointer to a double array
(2) int &rows : Number of rows in the array
(3) int &cols : Number of cols in the array
--------------------------------------------------------------*/
{
    for (int i = 0; i <length; ++i)
        printf("Array[%i] = %f\n", i, Array[i]);
}

void print_array_1d(int * __restrict Array, int length)
/*------------------------------------------------------------
Purpose: Print arrays
Inputs:
(1) int *Array : A pointer to a int array
(2) int &rows : Number of rows in the array
(3) int &cols : Number of cols in the array
--------------------------------------------------------------*/
{
    for (int i = 0; i <length; ++i)
        printf("Array[%i] = %i\n", i, Array[i]);

}

void printMatForm(double * __restrict Array, const int rows, const int cols, int width, const int precision)
/*------------------------------------------------------------
Purpose: Print arrays
Inputs:
(1) double *Array : A pointer to a double array
(2) int &rows : Number of rows in the array
(3) int &cols : Number of cols in the array
--------------------------------------------------------------*/
{
	if (rows == 1 && cols == 1)
		printf("Array[0][0] = %f\n", Array[0]);
	else
	{
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				printf(" %*.*f ", width, precision, Array[i*cols + j]);
			}
			printf("\n");
		}
	}

}

void printMatForm(int * __restrict Array, const int rows, const int cols, int width, const int precision)
/*------------------------------------------------------------
Purpose: Print arrays
Inputs:
(1) int *Array : A pointer to a double array
(2) int &rows : Number of rows in the array
(3) int &cols : Number of cols in the array
--------------------------------------------------------------*/
{
	if (rows == 1 && cols == 1)
		printf("Array[0][0] = %i\n", Array[0]);
	else
	{
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				printf(" %*.*i ", width, precision, Array[i*cols + j]);
			}
			printf("\n");
		}
	}

}
