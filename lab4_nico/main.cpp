#include <iostream>
#include <math.h>
#include <fstream>

#include "lab4.cpp"

using namespace std;

int main()
{
    double x_min = 0.;
    double x_max = 1.;
    int M = 10; // number of grid points in x
    double h = (x_max-x_min)/(M-1);

    double x[M];
    for (int j=0; j<M; j++) {
        x[j] = x_min + j*h;
    }

    int N = 10000; // number of timesteps
    double t_max = 0.001;
    double tau = t_max/N;

    double t[N+1];
    for (int n=0; n<N; n++) {
        t[n] = n*tau;
    }

    double* u = new double[M];

    // initial conditions
    for (int j=0; j<M; j++) {
        u[j] = 1;
    }

    // boundary conditions
    u[0] = 0;
    u[M-1] = 0;

    // specify the vectors for the Thomas algorithm
    // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    double a[M], b[M], c[M], z[M], d[M];
    for (int j=0; j<M; j++) {
        a[j] = 1/(h*h);
        b[j] = 1/tau - 2/(h*h);
        c[j] = 1/(h*h);
    }

    // animation constants
    int frames = 100;
    int frame_stride = N/frames;

    // write file header and initial condition
    cout << "iteration " << 0 << "/" << N;
    FILE* file = fopen("data/1d_simulation.csv", "w");
    fprintf(file, "{\"x_min\":%.5e,\"x_max\":%.5e,\"M\":%d}\n", x_min, x_max, M);
    fprintf(file, "time, u[0...M=%d]\n", M);
    fprintf(file, "%.5e, ", t[0]);
    for (int j=0; j<M; j++) {
        if (j==M-1) fprintf(file, "%.5e\n", u[j]);
        else fprintf(file, "%.5e, ", u[j]);
    }
    fclose(file);

    // step through the simulation
    for (int n=1; n<N+1; n++) {
        // calculate the right side of the system of equations
        for (int j=0; j<M; j++) {
            d[j] = u[j]/tau;
        }

        for (int j=0; j<M; j++) {
            a[j] = 1/(h*h);
            b[j] = 1/tau - 2/(h*h);
            c[j] = 1/(h*h);
        }

        // solve for u_(n) using u_(n-1)
        thomasAlgorithm(a, b, c, d, u, M);

        // apply boundary conditions
        u[0] = 0;
        u[M-1] = 0;

        // save frames to a files
        if (n % frame_stride == 0) {
            cout << "\riteration " << n << "/" << N;
            FILE* file = fopen("data/1d_simulation.csv", "a");
            fprintf(file, "%.5e, ", t[n]);
            for (int j=0; j<M; j++) {
                if (j==M-1) fprintf(file, "%.5e\n", u[j]);
                else fprintf(file, "%.5e, ", u[j]);
            }
            fclose(file);
        }
    }
    cout << endl;


    delete[] u;
    return 0;
}