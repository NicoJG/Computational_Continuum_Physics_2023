#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

// CONSTANTS
const double a = sqrt(3); // parameter of the problem

const int N_steps = 15; // number of steps to iterate
const int N_bases = 1000; // number of basis functions
const int N_x = 10000; // number of space points for numerical integration


double scalar_product(double beta_l, double* y_prime, double* x) {
    double scalar_prod = 0;

    // perform numerical integration using the trapezoidal rule
    double integrant_left = cos(beta_l * x[0]) * sqrt(1 + a*(y_prime[0])*(y_prime[0]));
    for (int i=0; i<(N_x-1); i++) {
        double integrant_right = cos(beta_l * x[i+1]) * sqrt(1 + a*(y_prime[i+1])*(y_prime[i+1]));
        scalar_prod += (x[i+1]-x[i])*(integrant_left + integrant_right)/2.;
        integrant_left = integrant_right;
    }

    return scalar_prod;
}

int main()
{
    double w[N_bases]; // vector of the weights

    double x[N_x]; // linspace of the space points for numerical integration
    for (int i=0; i<N_x; i++) {
        x[i] = -1 + 2./(N_x-1) * i;
    }

    double beta[N_bases]; // factor inside the cosine of each basis function
    for (int k=0; k<N_bases; k++) {
        beta[k] = (2.*k + 1)/2. * M_PI;
    }

    double y_at_0 = 0;
    double y_prime[N_x];
    double y[N_x];
    // initial guess y^0 = c*(x-1)*(x+1)
    double c = 1;
    for (int i=0; i<N_x; i++) {
        y[i] = c*(x[i]-1)*(x[i]+1); 
        y_prime[i] = 2*c*x[i];
    }
    y_at_0 = -c;
    cout << "y^(n=" << 0 << ")(x=0) = " << y_at_0 << endl;

    // step through the iterations
    for (int n=1; n<N_steps; n++) {
        // calculate the new weights
        for (int l=0; l<N_bases; l++) {
            w[l] = -1./(beta[l]*beta[l]) * scalar_product(beta[l], y_prime, x);
        }
        
        // calculate y(x) and y_prime(x)
        for (int i=0; i<N_x; i++) {
            y[i] = 0;
            y_prime[i] = 0;
            for (int l=0; l<N_bases; l++) {
                y[i] += w[l]*cos(beta[l]*x[i]);
                y_prime[i] += -beta[l]*w[l]*sin(beta[l]*x[i]);
            }
        }
        y_at_0 = 0;
        for (int l=0; l<N_bases; l++) {
            y_at_0 += w[l];
        }

        cout << "y^(n=" << n << ")(x=0) = " << y_at_0 << endl;
    }


    return 0;
}