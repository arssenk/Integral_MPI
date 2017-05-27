#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <map>

#include <sstream>

#include <mpi.h>

using namespace std;



template<class D>

double func_calculation(int m, double x1, double x2) {
    double sum1 = 0;
    double sum2 = 0;
    double g;
    for (int i = 1; i <= m; ++i)
    {
        sum1 += i * cos((i + 1) * x1 + 1);
        sum2 += i * cos((i + 1) * x2 + 1);
    }

    g = - sum1 * sum2;

    return g;
//    double g,sum;
//    int j;
//    for (int i = -2; i <= 2; ++i)
//    {
//        j = i;
//        sum += 1 / (5 * (i + 2) + j + 3 + pow(x1 - 16* j,6) + pow(x2 - 16* i,6));
//    }
//
//    g = pow(0.002 + sum, -1);
//
//    return g;

}

double integration(double x0, double x, double y0, double y, int m, double pr) {
    assert (m >= 5);
    double sum = 0;
    for (double i = x0; i <= x; i += pr) {
        for (double j = y0; j <= y; j += pr) {
            sum += func_calculation(m, i + pr / 2.0, j + pr / 2.0) * pr * pr;
        }
    }
    return sum;
}

//double thread_integration(double x0, double x, double y0, double y, int m, double pr, ) {
//
//    auto result = integration(x0, x, y0, y, m, pr);
//    lock_guard<mutex> lg(mx);
//    *r += result;
//    return *r+result;
//}

template <class T>
T get_param(string key, map<string, string> myMap) {
    istringstream ss(myMap[key]);
    T val;
    ss >> val;
    return val;
}

map<string, string> read_config(string filename) {
    string line, delimiter;
    ifstream myfile (filename);
    map<string, string> mp;

    delimiter = "=";


    if (myfile.is_open())
    {
        while (getline(myfile,line))
        {
            int pos = line.find(delimiter);
            string key = line.substr(0, pos);
            string value = line.substr(pos + 1);
            mp[key] = value;
        }

        myfile.close();
    }
    else {
        cout << "Error with opening the file!" << endl;
    }
    return mp;

};


int main()
{
    int commsize, rank, len;
    char procname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &len);

    if (rank == 0){
        string filename;
        filename = "config.txt";
        map<string, string> mp = read_config(filename);
        double abs_er, rel_er, x0, x1, y0, y1;

        int m;
        if (mp.size() != 0) {
            abs_er = get_param<double>("absol_er", mp);
            rel_er = get_param<double>("rel_er", mp);
            x0 = get_param<double>("x0", mp);
            x1 = get_param<double>("x1", mp);
            y0 = get_param<double>("y0", mp);
            y1 = get_param<double>("y1", mp);
            m = get_param<int>("m", mp);
            //
            // num_of_threads = get_param<int>("threads", mp);

            //thread threads[num_of_threads];

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
            unsigned double from_to[] = {x0, x1/commsize, y0, y1/commsize, m/commsize, pr};
            double start_time = MPI_Wtime();
            for(int i = 1; i<commsize; ++i)
            {
                printf("%2i: X:%10i - %10i; Y:%10i - %10i. M/commsize: %10i\n", i, from_to[0], from_to[1], from_to[2], from_to[3], from_to[4]);
                MPI_Send(from_to, 6, MPI_UNSIGNED_DOUBLE, i, 0, MPI_COMM_WORLD);
                from_to[0] = from_to[1]; //Може треба додавати 0.0001
                from_to[1] = from_to[0] + x1/commsize;
                from_to[2] = from_to[3];
                from_to[3] = from_to[2] + y1/commsize;
                from_to[4] = from_to[4] + from_to[4];
            }
#ifdef PRINT_PARTS
            printf("%2i: %10i - %10i\n", 0, from_to[0], max_number);
#endif // PRINT_PARTS
            double  res = integration(from_to[0], from_to[1]*commsize, from_to[2], from_to[3] * commsize,
                                      (int) (from_to[4] * commsize), pr);//Cast not shure
#ifdef PRINT_PARTS
            printf("Recv from %2i: %10i\n", 0, res);

#endif // PRINT_PARTS
            for(int i = 1; i<commsize; ++i)
            {
                unsigned double tr;
                MPI_Recv(&tr, 1, MPI_UNSIGNED_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef PRINT_PARTS
                printf("Recv from %2i: %10i\n", i, tr);
#endif // PRINT_PARTS

                res += tr;
            }
            printf("Result: %10i\n", res);
            printf("Total time: %g \n",   MPI_Wtime() - start_time);
    }else
        {
            double start_time = MPI_Wtime();
            unsigned double from_to[5]= {};
            MPI_Recv(from_to, 6, MPI_UNSIGNED_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            res = integration(from_to[0], from_to[1], from_to[2], from_to[3], from_to[4], from_to[5]);
            MPI_Send(&res, 1, MPI_UNSIGNED_DOUBLE, 0, 0, MPI_COMM_WORLD);
            printf("Process %d/%d execution time: %.6f\n", rank, commsize, MPI_Wtime() - start_time);
        }

        MPI_Finalize();
        return 0;
    }
}
