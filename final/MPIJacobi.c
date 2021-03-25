/*Some code was taken from https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html*/
/* Collin Campbell ~ MPI Jacobi Iteration ~ COMP 233A */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

void printArray(double **array, int row, int col); /* declaration of the print array function */

int main(int argc, char *argv[])
{
    int rows = 400;         /* amount of rows in the array */
    int cols = 800;         /* amount of cols in the array */
    int iterationCount = 0; /* init the iterationCount for later in the code, limits the interation to a set value */
    double diffNorm = 0.0;  /* the difference between new and old */
    double **tmp;           /* this will be used as an holder when switching the arrays */
    FILE *f;                /* makes a pointer for the reference to the file */
    double epsilon = 0.0;   /* takes an epislon value in order to stop interation */
    int rank = 0;           /* get the rank of the process */
    int size = 0;           /* get the size of the work */
    double gDiffNorm = 0.0; /* overall gdiffnorm */
    double t1 = 0.0;        /* timing vars */
    double t2 = 0.0;
    int row_first = 0;      /* row first for the local array */
    int row_last = 0;       /* the last row for all */
    int row = 0;            /* rest are interation varibles */
    int col = 0;    
    int r = 0;
    int c = 0;
    int i;
    int interation_stop = 0;

    MPI_Init(&argc, &argv);          //takes in two command line arguments
    epsilon = strtod(argv[1], NULL); /* takes one command line argument at sets it to the epsilon value */
    interation_stop = strtod(argv[2], NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /*grabs the processes rank*/
    MPI_Comm_size(MPI_COMM_WORLD, &size); /* how many processes there are */

    if (rank == 0)
    {
        printf("Collin Campbell ~ Serial Jacobi Iteration ~ COMP 233A\n");
        printf("~The purpose of this code is to run a jacobi iteration ~\n");
        printf("~in a parallel manner\n");
    }

    /* Note that top and bottom processes have one less row of interior
       points */
    row_first = 1;
    row_last = rows / size; /* cuts the into how many rows each process will have*/
    if (rank == 0)
        row_first++; /* if master interate the row_first */
    if (rank == size - 1)
        row_last--; /* otherwise the last process go ahead and  de-incremate row_last*/

    double **orgin = (double **)malloc(((rows / size) + 2) * sizeof(double *));        /* sets up the array for it to be populated by the init data */
    double *orgintemp = (double *)malloc(((rows / size) + 2) * cols * sizeof(double)); /* grabs the **orgin and turns it into a continous block of memory so MPI shuts up */
    for (i = 0; i < ((rows / size) + 2); i++)
    {
        orgin[i] = orgintemp + (i * cols);
    }

    double **updatedArray = (double **)malloc(((rows / size) + 2) * sizeof(double *));        /* sets up an array for it to be populated by new updated data from the orgin array*/
    double *updatedArraytemp = (double *)malloc(((rows / size) + 2) * cols * sizeof(double)); /* grabs the **updatedArray and turns it into a continous block of memory so MPI shuts up */
    for (i = 0; i < ((rows / size) + 2); i++)
    {
        updatedArray[i] = updatedArraytemp + (i * cols);
    }

    double **final = (double **)malloc(rows * sizeof(double *));          /* sets up the array for it to be populated by the init data */
    double *finaltemp = (double *)malloc((rows * cols) * sizeof(double)); /* grabs the **orgin and turns it into a continous block of memory so MPI shuts up */
    for (i = 0; i < rows; i++)
    {
        final[i] = finaltemp + (i * cols);
    }

    double **result = (double **)malloc(rows * sizeof(double *));                   /* sets up the array for it to be populated by the init data */
    double *resulttemp = (double *)malloc(((rows / size)) * cols * sizeof(double)); /* grabs the **orgin and turns it into a continous block of memory so MPI shuts up */
    for (i = 0; i < rows / size; i++)
    {
        result[i] = resulttemp + (i * cols);
    }

    /* fill the orgin array with data to start with the specs of 100 north, 100 east, 0 south, 0 west, all other cells are set to an init -1*/
    for (r = 0; r < (rows / size) + 2; r++)
    {
        for (c = 0; c < cols; c++)
        {
            if (r == row_first - 1)
            { /* if on the row 0 fill with north specified value */
                orgin[r][c] = 100;
                updatedArray[r][c] = 100;
            }
            else if (r == row_last + 1)
            { /* if on the last row then fill south specified value */
                orgin[r][c] = 0;
                updatedArray[r][c] = 0;
            }
            else if (c == 0)
            { /* if on the first column then specify with west specified value */
                orgin[r][c] = 0;
                updatedArray[r][c] = 0;
            }
            else if (c == cols - 1)
            { /* if of the last column then fill with east specified value */
                orgin[r][c] = 100;
                updatedArray[r][c] = 100;
            }
            else
            {
                orgin[r][c] = -1;
                updatedArray[r][c] = -1;
            }
        }
    }

    t1 = MPI_Wtime();
    iterationCount = 0; /* init the interations to 0 */
    do
    { /* does this loop untill it either hits 100 interation or the difference between two iterations is less that 1.0e-2 */
        /* Send up unless I'm at the top, then receive from below */
        /* Note the use of orgin[i] for &orgin[i][0] */
        if (rank < size - 1) /*if not last process send your array to the next process */
            MPI_Send(orgin[rows / size], cols, MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD);
        if (rank > 0) /* listening for the incoming rows and places them into orgin[0] */
            MPI_Recv(orgin[0], cols, MPI_DOUBLE, rank - 1, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /* Send down unless I'm at the bottom */
        if (rank > 0) /* all process except the master send your orgin[1] and send it to the process behind you */
            MPI_Send(orgin[1], cols, MPI_DOUBLE, rank - 1, 1,
                     MPI_COMM_WORLD);
        if (rank < size - 1) /* all the process execpt the last process listen for incoming information and place it orgin[maxn/size+1] */
            MPI_Recv(orgin[rows / size + 1], cols, MPI_DOUBLE, rank + 1, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /* Compute new values (but not on boundary) */
        iterationCount++;
        diffNorm = 0.0;

        for (row = row_first; row <= row_last; row++)
            for (col = 1; col < cols - 1; col++)
            {
                updatedArray[row][col] = (orgin[row][col + 1] + orgin[row][col - 1] + /* logic for caluculating the new value and the diff norm */
                                          orgin[row + 1][col] + orgin[row - 1][col]) /
                                         4.0;
                diffNorm += (updatedArray[row][col] - orgin[row][col]) *
                            (updatedArray[row][col] - orgin[row][col]);
            }
        /* Only transfer the interior points */
        tmp = orgin;          // originprt => tmpprt
        orgin = updatedArray; // updatedArrayprt => orgin effectivly taking updatedArray data and making its orgins data
        updatedArray = tmp;   // tmpprt => updatedArrayprt effectivly taking orgins data and making its updatedArrays data

        MPI_Allreduce(&diffNorm, &gDiffNorm, 1, MPI_DOUBLE, MPI_SUM, /* takes all the diffNorm and sum reduces it to gDiffNorm */
                      MPI_COMM_WORLD);
        gDiffNorm = sqrt(gDiffNorm); /* get the average from all the processes  */
        if (rank == 0 && iterationCount % 1000 == 0)
        {
            printf("At iteration %d, diff is %e\n", iterationCount,
                   gDiffNorm);
        }
    } while (gDiffNorm > epsilon && iterationCount < interation_stop);

    if (rank == 0)
    {
        for (r = 0; r < rows / size; r++)
        {
            for (c = 0; c < cols; c++)
            {
                final[r][c] = orgin[r + 1][c]; //this is the master copying its array to final to be printed
            }
        }
        int source = 0; //who's telling me stuff
        for (source = 1; source < size; source++)
        {
            MPI_Recv(result[0], ((rows / size) * cols), MPI_DOUBLE, source, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (r = 0; r < rows / size; r++)
            {
                for (c = 0; c < cols; c++)
                {
                    final[r + (rows / size * source)][c] = result[r][c]; //this is the master copying its array to final to be printed
                }
            }
        }
    }
    else
    {
        for (r = 0; r < rows / size; r++)
        {
            for (c = 0; c < cols; c++)
            {
                result[r][c] = orgin[r + 1][c];
            }
        }
        MPI_Send(result[0], ((rows / size) * cols), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }
    t2 = MPI_Wtime();
    if (rank == 0)
        printf("%d took %.2f seconds\n", size, (t2 - t1));

    MPI_Finalize(); /* wrapping up */

    if (rank == 0)
    {
        if (argv[3] != NULL)
        {
            /* Creating an PPM image */
            int const W = cols; //widht of the picture
            int const H = rows; //height of the picture
            int arrayTemp = 0;  //array value for that pixel
            int R = 0;          //color value red
            int G = 0;          //color value green always will be zero
            int B = 0;          //color value blue
            FILE *image;        //creating the file pointer
            image = fopen(argv[3], "w");
            fprintf(image, "P3\n#Collin Campbell ~ Serial Jacobi Iternation ~ COMP233\n#Jacobi Iteration, with EPSILON as 1e-2\n%d %d\n255\n", W, H);
            //Calculating the RGB values for the file
            for (r = 0; r < rows; r++)
            {
                for (c = 0; c < cols; c++)
                {
                    arrayTemp = final[r][c];
                    R = (255 * (arrayTemp) / 100);
                    B = 255 - R;
                    fprintf(image, "%d %d %d ", R, G, B);
                }
            }
            fclose(image); //closes the image file
        }

        /* Outputting to the timing file */
        FILE *f = fopen("./timingInfo.txt", "a");                    //opens the file
        fprintf(f, "The timing for %d was %.2f\n", size, (t2 - t1)); //writes to the file
        fclose(f);                                                   //closes the file
    }

    free(orgin); /* cleaning up the garbarge */
    free(orgintemp);
    free(final);
    free(finaltemp);
    free(updatedArray);
    free(updatedArraytemp);
}

/* 
    input - a pointer pointer referencing the array, how many rows the array has and how many col | output - NULL |
    functions as an easy way to print an array
*/
void printArray(double **array, int row, int col)
{
    int r; //interation int
    int c;
    for (r = 0; r < row; r++)
    {
        printf("\n");
        for (c = 0; c < col; c++)
        {
            printf("%6.2f|", array[r][c]);
        }
    }
}