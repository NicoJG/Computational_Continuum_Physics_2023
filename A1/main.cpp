#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

using namespace std;
 
using Eigen::MatrixXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
 
int main()
{
    double x_min = -1.;
    double x_max = 1.;

    int N = 4000; // number of timesteps
    int M = 1000; // number of spacesteps
    
    double h = (x_max - x_min) / (M - 0.5);
    double tau = 0.5*h;

    ArrayXd E(M);
    ArrayXd B(M);

    ArrayXXd E_arr(N,M);
    ArrayXXd B_arr(N,M);

    ArrayXd x_E = ArrayXd::LinSpaced(M, x_min,       x_min+(M-1.0)*h);
    ArrayXd x_B = ArrayXd::LinSpaced(M, x_min+0.5*h, x_min+(M-0.5)*h);

    // initial conditions:
    E = (abs(x_E) < 0.3).select(exp(-x_E*x_E/(2*0.1*0.1))*sin(20 * M_PI * x_E), 0 * x_E);
    B = (abs(x_B) < 0.3).select(exp(-x_B*x_B/(2*0.1*0.1))*sin(20 * M_PI * x_B), 0 * x_B);

    // save initial conditions to the output matrix
    E_arr.row(0) = E;
    B_arr.row(0) = B;

    for (int n=1; n<N; n++) {
        // FD step
        E.tail(M-1) = E.tail(M-1) - (tau/h) * (B.tail(M-1) - B.head(M-1));
        B.head(M-1) = B.head(M-1) - (tau/h) * (E.tail(M-1) - E.head(M-1));

        // reflective boundary conditions
        E.head(0) = 0;
        B.tail(0) = 0;
        
        // dampening boundary conditions (Z = 1)
        // E.head(0) =  B.head(0);
        // B.tail(0) = -E.tail(0);

        // save the current timestep
        E_arr.row(n) = E;
        B_arr.row(n) = B;
    }

    ofstream file("data/output_E.csv");
    if (file.is_open())
    {
        file << "# E field\n";
        file << E_arr << '\n';
    }
    file.close();

    file = ofstream("data/output_B.csv");
    if (file.is_open())
    {
        file << "# B field\n";
        file << B_arr << '\n';
    }
    file.close();
}