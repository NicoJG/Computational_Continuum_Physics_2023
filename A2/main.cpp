#include <complex.h>
#include "fftw3.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// complex number i = sqrt(-1)
const complex<double> i(0.0, 1.0);

double V(double x, double y) {
    return - 5.0 / pow(1 + (x/5)*(x/5) + (y/4)*(y/4), 4);
}

int main()
{

    // Space grid
    double    L = 10.0;         // Domain is -L < x < L,  -L < y < L
    const int M = 201;
    double    h = 2*L / (M-1);

    // Fourier space grid
    double    delta_k = 2 * M_PI / (2 * L);
    double    k_max   = 2 * M_PI / h;

    // Time grid
    double t_min = 0;
    double t_max = 100;
    int    N     = 1001;
    double tau   = (t_max - t_min) / (N-1);

    // space linspaces
    double* x = new double[M];
    double* y = new double[M];
    for (int m=0; m<M; m++) {
        x[m] = -L + h*m;
        y[m] = -L + h*m;
    }

    // wave vector linspaces
    double* k_x = new double[M];
    double* k_y = new double[M];
    for (int m=0; m<M; m++) {
        k_x[m] = -k_max + delta_k*m;
        k_y[m] = -k_max + delta_k*m;
    }
    
    auto psi_r = new complex<double>[M][M];
    auto psi_k = new complex<double>[M][M];

    // initial wave function
    for (int m_x=0; m_x<M; m_x++)
    for (int m_y=0; m_y<M; m_y++)
        psi_r[m_x][m_y] = exp(-(x[m_x]-1)*(x[m_x]-1) - (y[m_y]-1)*(y[m_y]-1)) / sqrt(M_PI);
    
    // stepping through the simulation
    for (int n=0; n<N; n++) {

        // Propagate in space by tau/2
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_r[m_x][m_y] *= exp(-i*V(x[m_x], y[m_y]) * tau/2.);

        // Compute fourier transform
        fftw_plan plan = fftw_plan_dft_2d(
            M, 
            M, 
            reinterpret_cast<fftw_complex*>(psi_r), 
            reinterpret_cast<fftw_complex*>(psi_k),
            1,
            0
        );
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    
        // Propagate in fourier domain by tau
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_k[m_x][m_y] *= exp(-i*(k_x[m_x]*k_x[m_x]+k_y[m_y]*k_y[m_y]) * tau);

        // Compute inverse Fourier transform
        plan = fftw_plan_dft_2d(
            M, 
            M, 
            reinterpret_cast<fftw_complex*>(psi_k), 
            reinterpret_cast<fftw_complex*>(psi_r),
            -1,
            0
        );
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        // Propagate in space by tau/2
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_r[m_x][m_y] *= exp(-i*V(x[m_x], y[m_y]) * tau/2.);

        cout << "iteration " << n << endl;
    }
    
    return 0;
}