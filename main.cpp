#include <iostream>
#include <cmath>

#include <mpi.h>

using namespace std;

double func_calculation(int m, double s, double sd);

double func_calculation(int m, double x1, double x2) {
    double sum1 = 0;
    double sum2 = 0;
    double g;
    for (int i = 1; i <= m; ++i) {
        sum1 += i * cos((i + 1) * x1 + 1);
        sum2 += i * cos((i + 1) * x2 + 1);
    }

    g = -sum1 * sum2;

    return g;
}

double integration(double x0, double x, double y0, double y, int m, double pr) {
    double sum = 0;
    for (double i = x0; i <= x; i += pr) {
        for (double j = y0; j <= y; j += pr) {
            sum += func_calculation(m, i + pr / 2.0, j + pr / 2.0) * pr * pr;
        }
    }
    return sum;
}

int main() {
    int commsize, rank, len;
    char procname[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &len);
    double res;
    if (rank == 0) {
        double abs_er, rel_er, x0, x1, y0, y1;

        int m;
        abs_er = 0.001;
        rel_er = 0.001;
        x0 = 0;
        x1 = 3;
        y0 = 0;
        y1 = 1;
        m = 5;

        double pr = 1E-3;

        double integral = 0;
        double interval_x = (x1 - x0) / 2;
        double x = x0;
        cout << "  Calculating...\n" << endl;
        double step1 = 1E-3;
        double step2 = step1 / 2.0;
        double integral1 = 0, integral2 = 0;
        double j = x, l = x;

        while (j < x1) {
            integral1 += integration(j, j + step1, y0, y1, m, pr);
            j += step1;
        }

        while (l < x1) {
            integral2 += integration(l, l + step2, y0, y1, m, pr);
            l += step2;
        }

        double abs_dif = abs(integral1 - integral2);
        double rel_dif = abs((integral1 - integral2) / max(integral1, integral2));

        if (abs_dif <= abs_er)
            cout << "| Absolute error is okay\t";
        else
            cout << "| Absolute error is not okay\t";

        cout << abs_dif << " vs " << abs_er << endl;

        if (rel_dif <= rel_er)
            cout << "| Relative error is okay\t";
        else
            cout << "| Relative error is not okay\t";
        cout << rel_dif << " vs " << rel_er << endl;
        double from_to[] = {x0, x1 / commsize, y0, y1 / commsize, m / commsize, pr};
        double start_time = MPI_Wtime();
        for (int i = 1; i < commsize; ++i) {
            MPI_Send(from_to, 6, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            from_to[0] = from_to[1]; //Може треба додавати 0.0001
            from_to[1] = from_to[0] + x1 / commsize;
            from_to[2] = from_to[3];
            from_to[3] = from_to[2] + y1 / commsize;
            from_to[4] = from_to[4] + from_to[4];
        }
#ifdef PRINT_PARTS
        printf("%2i: %10i - %10i\n", 0, from_to[0], max_number);
#endif // PRINT_PARTS
        res = integration(from_to[0], from_to[1] * commsize, from_to[2], from_to[3] * commsize,
                          (int) (from_to[4] * commsize), pr);//Cast not shure
#ifdef PRINT_PARTS
        printf("Recv from %2i: %10i\n", 0, res);

#endif // PRINT_PARTS
        for (int i = 1; i < commsize; ++i) {
            double tr;
            MPI_Recv(&tr, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef PRINT_PARTS
            printf("Recv from %2i: %10i\n", i, tr);
#endif // PRINT_PARTS

            res += tr;
        }
        printf("Result: %10i\n", res);
        printf("Total time: %g \n", MPI_Wtime() - start_time);
    } else {
        double start_time = MPI_Wtime();
        double from_to[5] = {};
        MPI_Recv(from_to, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        res = integration(from_to[0], from_to[1], from_to[2], from_to[3], from_to[4], from_to[5]);
        MPI_Send(&res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        printf("Process %d/%d execution time: %.6f\n", rank, commsize, MPI_Wtime() - start_time);
    }

    MPI_Finalize();
    return 0;

}
