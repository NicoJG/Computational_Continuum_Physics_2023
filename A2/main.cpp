#include <complex.h>
#include "fftw3.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// complex number i = sqrt(-1)
const complex<double> i(0.0, 1.0);

// Given potential and wavefunction for task
double V(double x, double y)    {
    return - 5. / pow(1. + (x/5.)*(x/5.) + (y/4.)*(y/4.), 4.);
}

// Test with harmonic oscillator
//double V(double x, double y) { 
//    return (x*x + y*y) / 4; 
//}

complex<double> psi0(double x, double y) {
    return exp(-(x-1.)*(x-1.) - (y-1.)*(y-1.)) / sqrt(M_PI);
}


int main()
{
    // Grid resolution
    const int M = 200;

    // Real space grid
    double    L  = 10.0;         // Domain is -L < x < L,  -L < y < L
    double    dx = 2.0 * L / M;

    // Fourier space grid
    double    dk    = 2 * M_PI     / (2.0 * L);
    double    k_max = 2 * M_PI * M / (2.0 * L);

    // Time grid
    double t_min = 0;
    double t_max = 100;
    const int N  = 1000;
    double tau   = (t_max - t_min) / N;

    // For saving snapshots of wavefunction to plotting
    int frames = 10; // Number of frames excluding initial
    int frame_stride = N/frames;

    // Real space coordinates
    double x[M], y[M];
    for (int m=0; m<M; m++) {
        x[m] = -L + dx * (m + 0.5);
        y[m] = -L + dx * (m + 0.5);
    }

    // Fourier space coordinates
    double k_y[M], k_x[M];
    for (int m=0; m<M; m++) {
        if (m<M/2) {
            k_x[m] = dk * m;
            k_y[m] = dk * m;
        } else {
            k_x[m] = dk * (m-M);
            k_y[m] = dk * (m-M);
        }
    }

    // For saving value of wave function over time
    int m_x_0p1 = round((L+0.1)/(2*L) * M)-1;
    int m_y_0p0 = round((L+0.0)/(2*L) * M)-1;
    printf("Storing wavefunction at (x,y) = (%.5f,%.5f)\n", 
            x[m_x_0p1], y[m_y_0p0]);
    complex<double> *psi_over_time     = new complex<double>[N];
    complex<double> *psi_time_spectrum = new complex<double>[N];

    // Create linspace for time and frequency
    double time[N], omega[N];
    for (int n=0; n<N; n++) {
        time[n] = t_min + tau * (n+1);
        omega[n] = (n<N/2) ? n     * (2. * M_PI) / (t_max - t_min) :
                             (n-N) * (2. * M_PI) / (t_max - t_min);
    }
    
    // Grids in real and fourier space with plans to transform between them
    auto psi_r = new complex<double>[M][M];
    auto psi_k = new complex<double>[M][M];
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
        psi_r[m_x][m_y] = psi0(x[m_x], y[m_y]);

    // Save initial wavefunction
    char filename[100];
    sprintf(filename, "data/psi_r_real_n=%05d_t=%.2f.csv", 0, 0.0);
    FILE* f = fopen(filename, "w");
    for (int m_x = 0; m_x < M; m_x++)
    for (int m_y = 0; m_y < M; m_y++) {
        if (m_y != M-1) fprintf(f, "%.5e,",  real(psi_r[m_x][m_y]));
        else            fprintf(f, "%.5e\n", real(psi_r[m_x][m_y]));
    }
    fclose(f);

    // stepping through the simulation
    for (int n=1; n<N+1; n++) {

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

        // Save updated wavefunction for plotting
        if (n % frame_stride == 0) {
            cout << "iteration " << n << "/" << N << endl;
            char filename[100];
            sprintf(filename, "data/psi_r_real_n=%05d_t=%.2f.csv", n, n * tau);
            FILE* f = fopen(filename, "w");
            for (int m_x = 0; m_x < M; m_x++)
            for (int m_y = 0; m_y < M; m_y++) {
                if (m_y != M-1) fprintf(f, "%.5e,",  real(psi_r[m_x][m_y]));
                else            fprintf(f, "%.5e\n", real(psi_r[m_x][m_y]));
            }
            fclose(f);
        }

        // Save value of wavefuntion at (x,y) = (0.1, 0)
        psi_over_time[n-1] = psi_r[m_x_0p1][m_y_0p0];
    }

    // Compute fourier transform in time
    printf("Fourier transforming\n");
    fftw_plan transform_t_to_omega = fftw_plan_dft_1d(N,
        reinterpret_cast<fftw_complex*>(&psi_over_time[0]), 
        reinterpret_cast<fftw_complex*>(&psi_time_spectrum[0]),
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(transform_t_to_omega);
    fftw_destroy_plan(transform_t_to_omega);

    // Save wavefunction at a point and spectrum
    printf("Saving to file\n");
    sprintf(filename, "data/psi_time_and_frequency.csv");
    f = fopen(filename, "w");
    fprintf(f, "time, omega, abs(psi(t)), arg(psi(t)), abs(psi(omega)), arg(psi(omega))\n");
    for (int n = 0; n < N; n++) {
        fprintf(f, "%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n", 
            time[n], 
            omega[n],
            abs(psi_over_time[n]), 
            arg(psi_over_time[n]),
            abs(psi_time_spectrum[n]), 
            arg(psi_time_spectrum[n])
        );
    }
    fclose(f);

    fftw_destroy_plan(transform_r_to_k);
    fftw_destroy_plan(transform_k_to_r);
    
    return 0;
}