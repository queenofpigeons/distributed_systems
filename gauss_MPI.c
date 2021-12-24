#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.0000000000001

double
gauss(double **matrix, int n);

void
swap(double *first, double *second);

int rank;
int size;

int
main(int argc, char *argv[])
{
    double **matrix;
    double determinant;
    int n;
    FILE *fd;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    #ifndef DEBUG
    if (rank == 0) {
        fd = fopen("input.txt", "r");
    }
    #else
    fd = stdin;
    #endif
    if (rank == 0) {
        if (fd == NULL) {
            fprintf(stderr, "Can not open input file.\n");
            exit(1);
        }
    }

    // Enter size of matrix
    if (rank == 0) {
        fscanf(fd, "%d", &n);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    matrix = calloc(n, sizeof(double *));
    if (matrix == NULL) {
        fprintf(stderr, "Out of memory, can't allocate matrix!");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = calloc(n, sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Out of memory, can't allocate matrix!");
            exit(1);
        }
    }

    // Enter matrix
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fscanf(fd, "%lf", &matrix[i][j]);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        MPI_Bcast(matrix[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    determinant = gauss(matrix, n);
    if (rank == 0) {
        printf("\nDeterminant of your matrix is %.3lf\n", determinant);

        fclose(fd);
    }

    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);


    MPI_Finalize();

    return 0;
}

double
gauss(double **matrix, int n)
{

    double start_time, finish_time;

    double lead;
    double determinant;

    int i, j, k;
    int row, sel;

    if (rank == 0) {
        start_time = MPI_Wtime();
    }
    determinant = 1;

    for (i = 0, row = 0; i < n - 1 && row < n; i++) {
        // Search for max element in column
        sel = i;
        for (j = i; j < n; j++) {
            if ( fabs(matrix[j][i]) > fabs(matrix[sel][i]) ) {
                sel = j;
            }
        }
        if (fabs(matrix[sel][i]) < EPS) {
            continue;
        }
        MPI_Bcast(&sel, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (sel != row) {
            if (rank == 0) {
                for (j = i; j < n; j++) {
                    swap(&matrix[sel][j], &matrix[row][j]);

                    if (j == i) {
                        #if DEBUG
                        printf("===>sel = %d\trow = %d\tj = %d\n", sel, row, j);
                        #endif
                        determinant *= -1;
                    }
                }
                if (sel % size != 0) {
                    MPI_Send(matrix[sel], n, MPI_DOUBLE, sel % size, 0, MPI_COMM_WORLD);
                    MPI_Send(matrix[row], n, MPI_DOUBLE, sel % size, 0, MPI_COMM_WORLD);
                }

                if (row % size != 0) {
                    MPI_Send(matrix[sel], n, MPI_DOUBLE, row % size, 0, MPI_COMM_WORLD);
                    MPI_Send(matrix[row], n, MPI_DOUBLE, row % size, 0, MPI_COMM_WORLD);
                }
            }
            if ((rank > 0) && (sel % size == rank)) {
                MPI_Recv(matrix[sel], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(matrix[row], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if ((rank > 0) && (row % size == rank)) {
                MPI_Recv(matrix[sel], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(matrix[row], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        for (j = i + 1; j < n; j++) {
            if (rank == 0) {
                lead = matrix[j][i] / matrix[i][i];
            }
            MPI_Bcast(&lead, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (j % size == rank) {
                if (j != row) {
                    for (k = i; k < n; k++) {
                        matrix[j][k] -= lead * matrix[i][k];
                    }
                }
            }
        }
        for (j = i + 1; j < n; j++) {
            if ((rank == 0) && (j % size != 0)) {
                MPI_Recv(matrix[j], n, MPI_DOUBLE, j % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if ((j % size == rank) && (j % size != 0)) {
                MPI_Send(matrix[j], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
        row++;
    }

    if (rank == 0) {
        finish_time = MPI_Wtime();
        printf("\n=========================\n");
        printf("Calculation time is %lf\n", finish_time - start_time);
        printf("=========================\n");

        for (i = 0; i < n; i++) {
            if ( isnan(matrix[i][i]) ) {
                determinant = 0;
                break;
            }
            determinant *= matrix[i][i];
        }
    }

    #if DEBUG
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
    #endif

    return determinant;
}

void
swap(double *first, double *second)
{
    double tmp = *first;
    *first = *second;
    *second = tmp;
}
