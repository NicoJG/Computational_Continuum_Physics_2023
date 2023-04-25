#include <complex.h>
#include "fftw3.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// complex number i = sqrt(-1)
const complex<double> i(0.0, 1.0);

double V(double x, double y) {
    return - 5. / pow(1. + (x/5.)*(x/5.) + (y/4.)*(y/4.), 4.);
}

//double d = 1. / M_PI;
//double V(double x, double y) {
//    return (x*x+y*y)/(d*d);
//}

int main()
{
    // Grid resolution
    const int M = 64;

    // Real space grid
    double    L  = 10.0;         // Domain is -L < x < L,  -L < y < L
    double    dx = 2.0 * L / M;

    // Fourier space grid
    double    dk    = 2 * M_PI     / (2.0 * L);
    double    k_max = 2 * M_PI * M / (2.0 * L);

    // Time grid
    double t_min = 0;
    double t_max = 100;
    int    N     = 1001;
    double tau   = (t_max - t_min) / (N-1);

    // For saving frames
    int frames = 20;
    int frame_stride = N/frames;

    // Real space coordinates
    double x[M], y[M];
    for (int m=0; m<M; m++) {
        x[m] = -L + dx * m;
        y[m] = -L + dx * m;
    }

    // Fourier space coordinates
    double k_y[M], k_x[M];
    for (int m=0; m<M; m++) {
        if (m<M/2) {
            k_x[m] = dk * m;
            k_y[m] = dk * m;
        }
        else {
            k_x[m] = dk * m - k_max;
            k_y[m] = dk * m - k_max;
        }
    }
    
    // Grids in real and fourier space with plans to transform between them
    complex<double> psi_r[M][M];
    complex<double> psi_k[M][M];
    fftw_plan transform_r_to_k = fftw_plan_dft_2d(M, M, 
        reinterpret_cast<fftw_complex*>(&psi_r[0][0]), 
        reinterpret_cast<fftw_complex*>(&psi_k[0][0]),
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan transform_k_to_r = fftw_plan_dft_2d(M, M, 
        reinterpret_cast<fftw_complex*>(&psi_k[0][0]), 
        reinterpret_cast<fftw_complex*>(&psi_r[0][0]),
        FFTW_BACKWARD, FFTW_ESTIMATE);

    // initial wave function
    for (int m_x=0; m_x<M; m_x++)
    for (int m_y=0; m_y<M; m_y++)
        //psi_r[m_x][m_y] = exp(-(x[m_x]*x[m_x]+y[m_y]*y[m_y])/(2*d));
        psi_r[m_x][m_y] = exp(-(x[m_x]-0)*(x[m_x]-0) - (y[m_y]-0)*(y[m_y]-0)) / sqrt(M_PI);
    
    // stepping through the simulation
    for (int n=0; n<N; n++) {

        // Save wavefunction for plotting
        if (n % frame_stride == 0) {
            cout << "iteration " << n << "/" << N-1 << endl;
            char filename[100];
            sprintf(filename, "data/psi_r_mag_t=%.2f.csv", n * tau);
            FILE* f = fopen(filename, "w");
            for (int m_y = 0; m_y < M; m_y++)
            for (int m_x = 0; m_x < M; m_x++) {
                if (m_x != M-1) fprintf(f, "%.5e,",  real(psi_r[m_x][m_y]));
                else            fprintf(f, "%.5e\n", real(psi_r[m_x][m_y]));
            }
            fclose(f);
        }


        // Propagate in space by tau/2
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_r[m_x][m_y] *= exp(-i * V(x[m_x], y[m_y]) * tau/2.);

        // Compute fourier transform
        fftw_execute(transform_r_to_k);
    
        // Propagate in fourier domain by tau
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_k[m_x][m_y] *= exp(-i * (k_x[m_x]*k_x[m_x] + k_y[m_y]*k_y[m_y]) * tau);
        

        // Compute inverse Fourier transform
        fftw_execute(transform_k_to_r);

        // Propagate in space by tau/2
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_r[m_x][m_y] *= exp(-i * V(x[m_x], y[m_y]) * tau/2.);

        // Normalize (fftw forward then backward will scale by M*M)
        for (int m_x=0; m_x<M; m_x++)
        for (int m_y=0; m_y<M; m_y++)
            psi_r[m_x][m_y] /= (M*M);
    }

    fftw_destroy_plan(transform_r_to_k);
    fftw_destroy_plan(transform_k_to_r);
    
    return 0;
}