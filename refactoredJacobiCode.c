/*code was taken from and edited upon from https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html*/

/* Collin Campbell ~ A Simple Jacobi Iteration ~ COMP 233A */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* This example handles a 12 x 12 mesh, on 4 processors only. */
#define maxn 12

int main(argc, argv) int argc; //command line input arguments
char **argv;
{
	int rank, value, size, errCnt, totErr, row, col, itCnt; /* never uses errcnt, toterr, and value, other than that all other vars are used */
	int row_first, row_last;								/* The first index that the process gets, then the row_last is the last indes the process gets */
	MPI_Status status;										/* creates a status var for signaling whether or not messaging was sucessful or not */
	double diffNorm, gDiffNorm;								/* The difference between the interations, diffNorm is the local difference and the gDiffNorm is the devation - https://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/*/
	double xLocal[(12 / 4) + 2][12];						/* makes an array for the local array for four processes */
	double xNew[(12 / 3) + 2][12];							/*new array for after the calculating the new change*/
	double epsilon = 0.0;									/* setting up the difference for the iterations to stop */

	MPI_Init(&argc, &argv);			 //takes in two command line arguments
	epsilon = strtod(argv[1], NULL); /* takes one command line argument at sets it to the epsilon value */

	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /*grabs the processes rank*/
	MPI_Comm_size(MPI_COMM_WORLD, &size); /* how many processes there are */

	if (size != 4)
		MPI_Abort(MPI_COMM_WORLD, 1); /* about if there are anything else than 4 processes */

	if (rank == 0)
	{
		printf("\nCollin Campbell ~ A Simple Jacobi Iteration ~ COMP 233A\n");
	}

	/* xLocal[][0] is lower ghostpoints, xLocal[][maxn+2] is upper */

	/* Note that top and bottom processes have one less row of interior
       points */
	row_first = 1;
	row_last = maxn / size; /* cuts the into how many rows each process will have*/
	if (rank == 0)
		row_first++; /* if master interate the row_first */
	if (rank == size - 1)
		row_last--; /* otherwise the last process go ahead and  de-incremate row_last*/

	/* 
    Representation of 34 - 38 lines of code
    
    master -
        will have ifirst = 2 and ilast = 3 with only 2 rows it has responsibility for
	slave 1 -
		will have rows ifirst = 1 and ilast = 3 with 3 rows it has responsibility for
	slave 2 -
		will have rows ifirst = 1 and ilast = 3 with 3 rows it has responsibility for
	slave 3 - 
		will have rows ifirst = 1 and ilast = 2 with 2 rows it has responsibility for

	the top edge and lower edge are always the same value and will never change only giving 10 rows of reponsibility
    */

	/* setting bounds to the conditions in the project North = 100, East = 75, South = 100, West = 10 */

	/* Fill the data as specified */
	for (row = 1; row <= maxn / size; row++) /* filling the array with data*/
		for (col = 0; col < maxn; col++)
			xLocal[row][col] = -1;
	for (col = 0; col < maxn; col++)
	{									  /* setting certain values in the array to a null value */
		xLocal[row_first - 1][col] = 100; /* setting the left side to be -1 */
		xLocal[row_last + 1][col] = 100;  /* setting the right side to be -1 */
	}
	for (row = row_first; row <= row_last; row++)
	{
		xLocal[row][0] = 75;  /* setting the east side to 75 */
		xLocal[row][11] = 10; /* setting the west side to 10 */
	}

	/* Fill the data as specified */
	// for (row=1; row<=maxn/size; row++) /* filling the array with data*/
	// for (col=0; col<maxn; col++)
	//     xLocal[row][col] = rank;
	// for (col=0; col<maxn; col++) { /* setting certain values in the array to a null value */
	// xLocal[row_first-1][col] = -1;
	// xLocal[row_last+1][col] = -1;
	// }

	/* 
	Representation of code line 55 to 61 with data point either -1 or rank of the xLocal array
	
	-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 
	3  3  3  3  3  3  3  3  3  3  3  3
	3  3  3  3  3  3  3  3  3  3  3  3
	2  2  2  2  2  2  2  2  2  2  2  2
	2  2  2  2  2  2  2  2  2  2  2  2
	2  2  2  2  2  2  2  2  2  2  2  2
	1  1  1  1  1  1  1  1  1  1  1  1
	1  1  1  1  1  1  1  1  1  1  1  1
	1  1  1  1  1  1  1  1  1  1  1  1
	0  0  0  0  0  0  0  0  0  0  0  0
	0  0  0  0  0  0  0  0  0  0  0  0
	-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 
		
	*/

	printf("\n\n");
	itCnt = 0; /* init the interations to 0 */
	do
	{ /* does this loop untill it either hits 100 interation or the difference between two iterations is less that 1.0e-2 */
		/* Send up unless I'm at the top, then receive from below */
		/* Note the use of xLocal[i] for &xLocal[i][0] */
		if (rank < size - 1) /*if not last process send your array to the next process */
			MPI_Send(xLocal[maxn / size], maxn, MPI_DOUBLE, rank + 1, 0,
					 MPI_COMM_WORLD);
		if (rank > 0) /* listening for the incoming rows and places them into xLocal[0] */
			MPI_Recv(xLocal[0], maxn, MPI_DOUBLE, rank - 1, 0,
					 MPI_COMM_WORLD, &status);
		/* Send down unless I'm at the bottom */
		if (rank > 0) /* all process except the master send your xLocal[1] and send it to the process behind you */
			MPI_Send(xLocal[1], maxn, MPI_DOUBLE, rank - 1, 1,
					 MPI_COMM_WORLD);
		if (rank < size - 1) /* all the process execpt the last process listen for incoming information and place it xLocal[maxn/size+1] */
			MPI_Recv(xLocal[maxn / size + 1], maxn, MPI_DOUBLE, rank + 1, 1,
					 MPI_COMM_WORLD, &status);

		/* Compute new values (but not on boundary) */
		itCnt++;
		diffNorm = 0.0;
		for (row = row_first; row <= row_last; row++)
			for (col = 1; col < maxn - 1; col++)
			{
				xNew[row][col] = (xLocal[row][col + 1] + xLocal[row][col - 1] + /* logic for caluculating the new value and the diff norm */
								  xLocal[row + 1][col] + xLocal[row - 1][col]) /
								 4.0;
				diffNorm += (xNew[row][col] - xLocal[row][col]) *
							(xNew[row][col] - xLocal[row][col]);
			}
		/* Only transfer the interior points */
		for (row = row_first; row <= row_last; row++)
			for (col = 1; col < maxn - 1; col++)
				xLocal[row][col] = xNew[row][col];

		MPI_Allreduce(&diffNorm, &gDiffNorm, 1, MPI_DOUBLE, MPI_SUM, /* takes all the diffNorm and sum reduces it to gDiffNorm */
					  MPI_COMM_WORLD);
		gDiffNorm = sqrt(gDiffNorm); /* get the average from all the processes  */
		if (rank == 0)
			printf("At iteration %d, diff is %e\n", itCnt,
				   gDiffNorm);
	} while (gDiffNorm > epsilon && itCnt < 100);

	/* printing the final result of the plate itself, to lazy to figure out how it would send it to master and put it into the xLocal correctly*/
	int turnToPrint = 0;

	if (rank == 0)
	{
		printf("\nThe finilized array of the plate is\n");
	}

	for (turnToPrint = 0; turnToPrint < size; turnToPrint++)
	{
		if (rank == turnToPrint)
		{
			for (int row = 1; row <= maxn / size; row++)
			{
				printf("\n");
				for (int col = 0; col < maxn; col++)
				{
					printf("%.2f	", xLocal[row][col]);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD); /* makes sure that each process stays together for the printing so no one process gets ahead */
	}

	MPI_Finalize(); /* wrapping up */
	if (rank == 0)
	{
		printf("\n\n\n****Finished Without Issue****\n");
	}
	return 0;
}