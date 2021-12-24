#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi-ext.h"
#include <signal.h>
#include <sys/time.h>

#define EPS 0.0000000000001

double
gauss(double **matrix, int n);

void
swap(double *first, double *second);

MPI_Comm main_comm = MPI_COMM_WORLD;
int rank;
int size;
double **matrix;
bool err_happens = false;
int n;

static void save_checkpoint()
{
    if (rank == 0) {
        FILE* f = fopen("gauss_data.txt", "w");
        for(int i = 0; i < n; i++) {
            for (int j =0; j < n; j++) {
                fprintf(f, "%lf", matrix[i][j]);
            }
        }
        fclose(f);
    }
    MPI_Barrier(main_comm);
}

static void load_checkpoint()
{
    FILE* f = fopen("gauss_data.txt", "r");

    for(int i = 0; i < n; i++) {
        for (int j =0; j < n; j++) {
            fscanf(f, "%lf", &matrix[i][j]);
        }
    }
    fclose(f);
    printf("Proc %d loaded checkpoint\n", rank);
    MPI_Barrier(main_comm);
}

static void verbose_errhandler(MPI_Comm* pcomm, int* perr, ...)
{
    MPI_Comm comm = *pcomm;
    int err = *perr;
    char errstr[MPI_MAX_ERROR_STRING];
    int i, err_rank, err_size, nf, len, eclass;
    MPI_Group group_c, group_f;
    int *ranks_gc, *ranks_gf;

    MPI_Error_class(err, &eclass);
    if( MPIX_ERR_PROC_FAILED != eclass ) {
        MPI_Abort(comm, err);
    }

    MPI_Comm_rank(comm, &err_rank);
    MPI_Comm_size(comm, &err_size);

    MPIX_Comm_failure_ack(comm);
    MPIX_Comm_failure_get_acked(comm, &group_f);
    MPI_Group_size(group_f, &nf);
    MPI_Error_string(err, errstr, &len);

    printf("Rank %d / %d: Notified of error %s. %d found dead: { ",
           err_rank, err_size, errstr, nf);

    ranks_gf = (int*)malloc(nf * sizeof(int));
    ranks_gc = (int*)malloc(nf * sizeof(int));
    MPI_Comm_group(comm, &group_c);
    for(i = 0; i < nf; i++)
        ranks_gf[i] = i;
    MPI_Group_translate_ranks(group_f, nf, ranks_gf,
                              group_c, ranks_gc);
    for(i = 0; i < nf; i++)
        printf("%d ", ranks_gc[i]);
    printf("}\n");

    MPIX_Comm_shrink(comm, &main_comm);
    MPI_Comm_rank(main_comm, &rank);
    MPI_Comm_size(main_comm, &size);
    load_checkpoint();

    err_happens = true;
    free(ranks_gc);
    free(ranks_gf);
}

int
main(int argc, char *argv[])
{
    double determinant;
    FILE *fd;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(main_comm, &size);
    MPI_Comm_rank(main_comm, &rank);

    MPI_Errhandler errh;
    MPI_Comm_create_errhandler(verbose_errhandler, &errh);
    MPI_Comm_set_errhandler(main_comm, errh);

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

    MPI_Bcast(&n, 1, MPI_INT, 0, main_comm);

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
        MPI_Bcast(matrix[i], n, MPI_DOUBLE, 0, main_comm);
    }
    MPI_Barrier(main_comm);

    save_checkpoint();
    MPI_Barrier(main_comm);

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
    bool first_itter = false;

    double lead;
    double determinant;

    int i, j, k;
    int row, sel;

    if (rank == 0) {
        start_time = MPI_Wtime();
    }
    determinant = 1;


    for (i = 0, row = 0; i < n - 1 && row < n; i++) {
        first_itter = true;
        while (err_happens || first_itter) {
            err_happens = false;
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
            MPI_Bcast(&sel, 1, MPI_INT, 0, main_comm);
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
                        MPI_Send(matrix[sel], n, MPI_DOUBLE, sel % size, 0, main_comm);
                        MPI_Send(matrix[row], n, MPI_DOUBLE, sel % size, 0, main_comm);
                    }

                    if (row % size != 0) {
                        MPI_Send(matrix[sel], n, MPI_DOUBLE, row % size, 0, main_comm);
                        MPI_Send(matrix[row], n, MPI_DOUBLE, row % size, 0, main_comm);
                    }
                }
                if ((rank > 0) && (sel % size == rank)) {
                    MPI_Recv(matrix[sel], n, MPI_DOUBLE, 0, 0, main_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(matrix[row], n, MPI_DOUBLE, 0, 0, main_comm, MPI_STATUS_IGNORE);
                }
                if ((rank > 0) && (row % size == rank)) {
                    MPI_Recv(matrix[sel], n, MPI_DOUBLE, 0, 0, main_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(matrix[row], n, MPI_DOUBLE, 0, 0, main_comm, MPI_STATUS_IGNORE);
                }
            }

            if (!err_happens)
                MPI_Barrier(main_comm);

            for (j = i + 1; j < n; j++) {
                if (rank == 0) {
                    lead = matrix[j][i] / matrix[i][i];
                }
                MPI_Bcast(&lead, 1, MPI_DOUBLE, 0, main_comm);
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
                    MPI_Recv(matrix[j], n, MPI_DOUBLE, j % size, 0, main_comm, MPI_STATUS_IGNORE);
                }
                if ((j % size == rank) && (j % size != 0)) {
                    MPI_Send(matrix[j], n, MPI_DOUBLE, 0, 0, main_comm);
                }
            }
            row++;

            first_itter = false;
            if (!err_happens)
                MPI_Barrier(main_comm);
        }
        save_checkpoint();
        MPI_Barrier(main_comm);
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
