/*Some code was taken from https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html*/
/* Collin Campbell ~ Serail Jacobi Iteration ~ COMP 233A */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void printArray(double **array, int row, int col); /* declaration of the print array function */

void garbageCollection(double **array, int row); /* declaration of the garbage collection function */

int main(int argc, char *argv[])
{
    printf("Collin Campbell ~ Serial Jacobi Iteration ~ COMP 233A\n");
    printf("~The purpose of this code is to run a jacobi iteration ~\n");
    printf("~in a serial manner\n");

    int row = 400;                               /* amount of rows in the array */
    int col = 800;                               /* amount of cols in the array */
    double epsilon = strtod(argv[1], NULL);      /* the value for the process to stop once it reaches a diff of less than current value, taken from the second cmd line arg*/
    int interation_stop = strtod(argv[2], NULL); /* init the iterationCount for later in the code, limits the interation to a set value */
    double diffNorm = 0.0;                       /* the difference between new and old */
    double **tmp;                                /* this will be used as an holder when switching the arrays */
    FILE *f;                                     /* makes a pointer for the reference to the file */
    clock_t t1 = 0.0;                            /* inits the clock time */
    int iterationCount = 0;                      /* init the iterationCount */
    int k;                                       /* iteration vars */
    int i;
    int r;
    int c;

    double **orgin = (double **)malloc(row * sizeof(double *)); /* sets up the array for it to be populated by the init data */
    for (i = 0; i < row; i++)
    {
        orgin[i] = (double *)malloc(col * sizeof(double));
    }

    double **updatedArray = (double **)malloc(row * sizeof(double *)); /* sets up an array for it to be populated by new updated data from the orgin array*/
    for (k = 0; k < row; k++)
    {
        updatedArray[k] = (double *)malloc(col * sizeof(double));
    }

    /* fill the orgin array with data to start with the specs of 100 north, 100 east, 0 south, 0 west, all other cells are set to an init -1*/
    for (r = 0; r < row; r++)
    {
        for (c = 0; c < col; c++)
        {
            if (r == 0)
            { /* if on the row 0 fill with north specified value */
                orgin[r][c] = 100;
                updatedArray[r][c] = 100;
            }
            else if (r == row - 1)
            { /* if on the last row then fill south specified value */
                orgin[r][c] = 0;
                updatedArray[r][c] = 0;
            }
            else if (c == 0)
            { /* if on the first column then specify with west specified value */
                orgin[r][c] = 0;
                updatedArray[r][c] = 0;
            }
            else if (c == col - 1)
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
    t1 = clock();
    /* do the calc on the orgin plate and update it to the new plate and keep switching the pointers until max iterations or epsilon is true */
    do
    {
        iterationCount++; /* interates the iteration count */
        diffNorm = 0.0;   /* need to reset the diffNorm everytime or else we'll get the previous interations also included into the calc */
        for (r = 1; r < row - 1; r++)
        { /* set up the iteration so it doesnt touch the edges that are suppose to stay static */
            for (c = 1; c < col - 1; c++)
            {
                updatedArray[r][c] = (orgin[r][c + 1] + orgin[r][c - 1] + orgin[r + 1][c] + orgin[r - 1][c]) / 4; /* grabs the north, south, east, and west cell and update the array */
                diffNorm = diffNorm + (updatedArray[r][c] - orgin[r][c]) * (updatedArray[r][c] - orgin[r][c]);    /* grab the diffnorm between the two plates */
            }
        }
        diffNorm = sqrt(diffNorm); /* get the average */
        /* make the new plate the old by switching the pointers of the array effectivly making the new the old and the old forgotten in a way since it's going to be overwritten */
        tmp = orgin;          // originprt => tmpprt
        orgin = updatedArray; // updatedArrayprt => orgin effectivly taking updatedArray data and making its orgins data
        updatedArray = tmp;   // tmpprt => updatedArrayprt effectivly taking orgins data and making its updatedArrays data

        // for (int r = 1; r < row - 1; r++)
        // {
        //     for (int c = 1; c < col - 1; c++)
        //     {
        //         orgin[r][c] = updatedArray[r][c];
        //     }
        // }
        if (iterationCount % 1000 == 0)
        {
            printf("At interation %d, the diff is %e\n", iterationCount, diffNorm);
        }
    } while (epsilon < diffNorm && iterationCount < interation_stop);
    t1 = clock() - t1; //timing for the clock

    printf("Serial time took %.2f seconds\n", ((double)(t1) / CLOCKS_PER_SEC));
    f = fopen("./timingInfo.txt", "a");
    fprintf(f, "Serial time took %.2f seconds\n", ((double)(t1) / CLOCKS_PER_SEC));
    fclose(f); //closes the file for the timing

    /* Creating an PPM image */
    if (argv[3] != NULL)
    {
        int const W = col; //widht of the picture
        int const H = row; //height of the picture
        int arrayTemp = 0; //array value for that pixel
        int R = 0;         //color value red
        int G = 0;         //color value green always will be zero
        int B = 0;         //color value blue
        FILE *image;       //creating the file pointer
        image = fopen(argv[3], "w");
        fprintf(image, "P3\n#Collin Campbell ~ Serial Jacobi Iternation ~ COMP233\n#Jacobi Iteration, with EPSILON as 1e-2\n%d %d\n255\n", W, H);
        //Calculating the RGB values for the file
        for (r = 1; r < row - 1; r++)
        {
            for (c = 1; c < col - 1; c++)
            {
                arrayTemp = orgin[r][c];
                R = (255 * (arrayTemp) / 100);
                B = 255 - R;
                fprintf(image, "%d %d %d ", R, G, B);
            }
        }
        fclose(image); //closes the image file
    }
    garbageCollection(orgin, row);        //cleaning out the orgin array
    garbageCollection(updatedArray, row); //cleaning out the updatedArray
}

/* 
    input - a pointer pointer referencing the array, how many rows the array has and how many col | output - NULL |
    functions as an easy way to print an array
*/
void printArray(double **array, int row, int col)
{
    int r; //iteration vars
    int c;
    for (r = 0; r < row; r++)
    {
        printf("\n");
        for (c = 0; c < col; c++)
        {
            printf("%6.2f ", array[r][c]);
        }
    }
}

/*  
    input - a pointer pointer referencing the array and the amount of rows the array has| output - NULL |
    functions as an easy way to clean out the array from memory, as long as the array is a prtprt array
*/
void garbageCollection(double **array, int row)
{
    int i;
    for (i = 0; i < row; i++)
    {
        free(array[i]);
    }
    free(array);
}