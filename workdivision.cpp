#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "matrix_def.h"
#include "memory.h"


struct iarray limits(int tcov_rows, int numprocs)
{

	/*
	This function distributes rows as evenly as possible among different processors
	based on number of rows to be processed in an OpenMp or/and MPI environment.

	Inputs:
	tcov_rows: no of rows that would be divided across number of processors as evenly
	as possible
	numprocs: how many processors or nodes that would be used for subDivision of work

	Output:
	Returns a structure that gives the start and the end index of the
	rows to be processed by different processors or threads
	*/

	struct iarray rowDistribute = { 0 };
	int remainder, division, index, i;
	int pset = 0;

	if (numprocs <= 0)
	{
		printf("Error! Number of Processors Cannot be Zero. Program Exiting");
		exit(0);
	}

	if (numprocs > tcov_rows)
	{
		numprocs = tcov_rows;
		printf("Number of processors greater than number of rows\n");
		printf("Number of processors are set to number of rows\n");
	}

	rowDistribute.rows = numprocs;
	rowDistribute.cols = 2;
	rowDistribute.array = (int*)calloc(numprocs * 2, sizeof(int));

	if (rowDistribute.array == NULL)
	{
		printf("Error! memory not allocated for start vector, Program Exiting");
		exit(0);
	}

	remainder = tcov_rows % numprocs;
	if (remainder == 0)
		pset = 0;
	else if (remainder > 0)
		pset = 1;

	division = (tcov_rows - remainder) / (numprocs);
	index = 0;

	for (i = 0; i < numprocs; ++i)
	{
		/* Use this branch if rows can be evenly divided across processors*/
		if (pset == 0)
		{
			rowDistribute.array[i] = index;
			rowDistribute.array[i + rowDistribute.rows] = rowDistribute.array[i]
				+ division - 1;
			index = rowDistribute.array[i + rowDistribute.rows] + 1;
		}
		/*
		Use this branch if rows cannot be evenly divided across processors.
		First evenly divide rows across processors and then add the remaining
		rows to the start and finish arrays as evenly as possible by adding one
		row at a time
		*/
		else if (pset == 1)
		{
			if (remainder > 0)
			{
				rowDistribute.array[i] = index;
				rowDistribute.array[i + rowDistribute.rows] = index + (division - 1)
					+ (remainder - (remainder - 1));
				remainder = remainder - 1;
				index = rowDistribute.array[i + rowDistribute.rows] + 1;
			}
			else if (remainder == 0)
			{
				rowDistribute.array[i] = index;
				rowDistribute.array[i + rowDistribute.rows] = index + division - 1;
				index = rowDistribute.array[i + rowDistribute.rows] + 1;
			}
		}
	}
	return rowDistribute;
}
