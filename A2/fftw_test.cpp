#include <complex>
#include "fftw3.h"
#include <iostream>
#include <cmath>

using namespace std;

const complex<double> i = complex<double>(0, 1);

int main() {

    const int N = 16;
    double    x_max = 1.0;
    double    dx    = x_max/N;
    double    k_max = 2.0 * M_PI * N / x_max;
    double    dk    = 2.0 * M_PI / x_max;

    double x[N], k[N];
    for (int n=0; n<N; n++) x[n] = n * dx;
    for (int n=0; n<N; n++) k[n] = n * dk;

    complex<double> phi_x[N][N], phi_k[N][N];

    for (int m=0; m<N; m++)
    for (int n=0; n<N; n++)
        phi_x[m][n] = exp( i * (3. * dk * x[m] + 2. * dk * x[n]));

    printf("\nphi(r) = \n");
    for (int n=0; n<N; n++) {
        for (int m=0; m<N; m++) {
            printf("%.2f ,", real(phi_x[m][n]));
        }
        printf("\n");
    }

    fftw_plan plan = fftw_plan_dft_2d(N, N,
        reinterpret_cast<fftw_complex*>(&phi_x[0][0]), 
        reinterpret_cast<fftw_complex*>(&phi_k[0][0]), 
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    printf("\nphi(k) = \n");
    for (int n=0; n<N; n++) {
        for (int m=0; m<N; m++) {
            printf("%.2f ,", abs(phi_k[m][n]));
        }
        printf("\n");
    }

    plan = fftw_plan_dft_2d(
        N, N,
        reinterpret_cast<fftw_complex*>(&phi_k[0][0]), 
        reinterpret_cast<fftw_complex*>(&phi_x[0][0]), 
        FFTW_BACKWARD, 
        FFTW_ESTIMATE
    );
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    

    // Normalize
    for (int m=0; m<N; m++)
    for (int n=0; n<N; n++) 
        phi_x[m][n] /= (N*N);

    printf("\nphi(r) = \n");
    for (int n=0; n<N; n++) {
        for (int m=0; m<N; m++) {
            printf("%.2f ,", real(phi_x[m][n]));
        }
        printf("\n");
    }
}