#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 4
#define ARRLEN 16

int rank;
int tasks;
MPI_Comm comm;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);

    int size[2] = {SIZE, SIZE};
    int periodic[2] = {0};
    // создание транспьютерной матрицы
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 0, &comm);
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    srand(rank + 4);
    int arr[ARRLEN];
    for (int i = 0; i < ARRLEN; i++) {
        arr[i] = rand() % 1000;
        printf("proc [%d][%d]; arr[%d] = %d\n", coords[0], coords[1], i, arr[i]);
    }
    printf("Coordinates for process %d: (%d, %d)\n", rank, coords[0], coords[1]);

    int arr_result[ARRLEN];
    int local_max;
    int other_coords[2];
    int other_rank = 0;
    MPI_Status status;

    // шаг 1
    other_coords[1] = coords[1];
    if (coords[0] == 0 || coords[0] == 3) {
        if (coords[0] == 0) {
            other_coords[0] = coords[0] + 1;
        } else {
            other_coords[0] = coords[0] - 1;
        }
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(arr, ARRLEN, MPI_INT, other_rank, 0, comm);
    } else {
        if (coords[0] == 1) {
            other_coords[0] = coords[0] - 1;
        } else {
            other_coords[0] = coords[0] + 1;
        }
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(arr_result, ARRLEN, MPI_INT, other_rank, 0, comm, &status);
        for (int i = 0; i < ARRLEN; i++) {
            if (arr_result[i] > arr[i]) {
                arr[i] = arr_result[i];
            }
        }
    }
    MPI_Barrier(comm);
    // шаг 2
    other_coords[0] = coords[0];
    if (coords[0] == 1 || coords[0] == 2) {
        if (coords[1] == 0 || coords[1] == 3) {
            if (coords[1] == 0) {
                other_coords[1] = coords[1] + 1;
            }
            if (coords[1] == 3) {
                other_coords[1] = coords[1] - 1;
            }
            MPI_Cart_rank(comm, other_coords, &other_rank);
            MPI_Send(arr, ARRLEN, MPI_INT, other_rank, 0, comm);
        } else {
            if (coords[1] == 1) {
                other_coords[1] = coords[1] - 1;
            }
            if (coords[1] == 2) {
                other_coords[1] = coords[1] + 1;
            }
            MPI_Cart_rank(comm, other_coords, &other_rank);
            MPI_Recv(arr_result, ARRLEN, MPI_INT, other_rank, 0, comm, &status);
            for (int i = 0; i < ARRLEN; i++) {
                if (arr_result[i] > arr[i]) {
                    arr[i] = arr_result[i];
                }
            }
        }
    }
    MPI_Barrier(comm);
    // шаг 3
    other_coords[1] = coords[1];
    if ((coords[0] == 1 || coords[0] == 2) && (coords[1] == 1 || coords[1] == 2)) {
        if (coords[0] == 1) {
            other_coords[0] = coords[0] + 1;
        }
        if (coords[0] == 2) {
            other_coords[0] = coords[0] - 1;
        }
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Sendrecv(arr, ARRLEN, MPI_INT, other_rank, 0, arr_result, ARRLEN, MPI_INT, other_rank, 0, comm, &status);
        for (int i = 0; i < ARRLEN; i++) {
            if (arr_result[i] > arr[i]) {
                arr[i] = arr_result[i];
            }
        }
    }
    MPI_Barrier(comm);

    // шаг 4
    other_coords[0] = coords[0];
    if ((coords[0] == 1 || coords[0] == 2) && (coords[1] == 1 || coords[1] == 2)) {
        if (coords[1] == 1) {
            other_coords[1] = coords[1] + 1;
        }
        if (coords[1] == 2) {
            other_coords[1] = coords[1] - 1;
        }
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Sendrecv(arr, ARRLEN, MPI_INT, other_rank, 0, arr_result, ARRLEN, MPI_INT, other_rank, 0, comm, &status);
        for (int i = 0; i < ARRLEN; i++) {
            if (arr_result[i] > arr[i]) {
                arr[i] = arr_result[i];
            }
        }
    }
    MPI_Barrier(comm);

    // шаг 5
    if (coords[0] == 1 && coords[1] == 1) {
        local_max = arr[5];
        other_coords[0] = 0;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[0], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 1;
        other_coords[1] = 0;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[4], 1, MPI_INT, other_rank, 0, comm);
    }

    if ((coords[0] == 1 && coords[1] == 0) || (coords[0] == 0 && coords[1] == 1)) {
        other_coords[0] = 1;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        if (coords[0] == 1) {
            MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
        } else {
            MPI_Recv(&arr[0], 1, MPI_INT, other_rank, 0, comm, &status);
        }
    }

    if (coords[0] == 1 && coords[1] == 2) {
        local_max = arr[6];
        other_coords[0] = 0;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[3], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 1;
        other_coords[1] = 3;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[7], 1, MPI_INT, other_rank, 0, comm);
    }

    if ((coords[0] == 1 && coords[1] == 3) || (coords[0] == 0 && coords[1] == 2)) {
        other_coords[0] = 1;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        if (coords[0] == 1) {
            MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
        } else {
            MPI_Recv(&arr[3], 1, MPI_INT, other_rank, 0, comm, &status);
        }
    }

    if (coords[0] == 2 && coords[1] == 1) {
        local_max = arr[9];
        other_coords[0] = 2;
        other_coords[1] = 0;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[8], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 3;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[12], 1, MPI_INT, other_rank, 0, comm);
    }

    if ((coords[0] == 2 && coords[1] == 0) || (coords[0] == 3 && coords[1] == 1)) {
        other_coords[0] = 2;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        if (coords[0] == 2) {
            MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
        } else {
            MPI_Recv(&arr[12], 1, MPI_INT, other_rank, 0, comm, &status);
        }
    }

    if (coords[0] == 2 && coords[1] == 2) {
        local_max = arr[10];
        other_coords[0] = 2;
        other_coords[1] = 3;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[11], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 3;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[15], 1, MPI_INT, other_rank, 0, comm);
    }

    if ((coords[0] == 2 && coords[1] == 3) || (coords[0] == 3 && coords[1] == 2)) {
        other_coords[0] = 2;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        if (coords[0] == 2) {
            MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
        } else {
            MPI_Recv(&arr[15], 1, MPI_INT, other_rank, 0, comm, &status);
        }
    }

    MPI_Barrier(comm);

    // шаг 6
    if (coords[0] == 1 && coords[1] == 1) {
        other_coords[0] = 0;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[1], 1, MPI_INT, other_rank, 0, comm);
    }

    if (coords[0] == 0 && coords[1] == 1) {
        other_coords[0] = 0;
        other_coords[1] = 0;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[0], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 1;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

    if (coords[0] == 0 && coords[1] == 0) {
        other_coords[0] = 0;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

    //--------------------


    if (coords[0] == 1 && coords[1] == 2) {
        other_coords[0] = 0;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[2], 1, MPI_INT, other_rank, 0, comm);
    }

    if (coords[0] == 0 && coords[1] == 2) {
        other_coords[0] = 0;
        other_coords[1] = 3;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[3], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 1;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

    if (coords[0] == 0 && coords[1] == 3) {
        other_coords[0] = 0;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

        //--------------------

    if (coords[0] == 2 && coords[1] == 1) {
        other_coords[0] = 3;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[13], 1, MPI_INT, other_rank, 0, comm);
    }

    if (coords[0] == 3 && coords[1] == 1) {
        other_coords[0] = 3;
        other_coords[1] = 0;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[12], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 2;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

    if (coords[0] == 3 && coords[1] == 0) {
        other_coords[0] = 3;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }
        //--------------------

    if (coords[0] == 2 && coords[1] == 2) {
        other_coords[0] = 3;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[14], 1, MPI_INT, other_rank, 0, comm);
    }

    if (coords[0] == 3 && coords[1] == 2) {
        other_coords[0] = 3;
        other_coords[1] = 3;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&arr[15], 1, MPI_INT, other_rank, 0, comm);
        other_coords[0] = 2;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

    if (coords[0] == 3 && coords[1] == 3) {
        other_coords[0] = 3;
        other_coords[1] = 2;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&local_max, 1, MPI_INT, other_rank, 0, comm, &status);
    }

    MPI_Barrier(comm);

    // вывод результата
    printf("Max result for %d: %d\n", rank, local_max);
    MPI_Finalize();
    return 0;
}