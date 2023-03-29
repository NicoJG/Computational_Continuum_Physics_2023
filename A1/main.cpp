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

    double t_min = 0.;
    double t_max = 5.;

    int N = 3000; // number of timesteps
    int M = 1000; // number of spacesteps
    
    double h = (x_max - x_min) / (M - 0.5);
    double tau = (t_max - t_min) / (N - 1);

    // check if stability condition (sigma = tau/h < 1) is true
    if (tau/h > 1) {
        cerr << "ERROR! Method is instable: tau = " << tau << " h = " << h << "\n";
        exit(1);
    }

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

        // FD step E
        double E_old = E[1];
        E.tail(M-1) = E.tail(M-1) - (tau/h) * (B.tail(M-1) - B.head(M-1));
        
        // Boundary conditions for E field (left)
        //E[0] = E_old+(tau-h)/(tau + h)*(E[1]-E[0]);     // absorbing
        //E[0] = 0;                                       // reflective
        //E[0] = E[M-1];                                  // periodic

        // Antenna
        double A = 0.5;
        double omega = 10; 
        E[0] = A*sin(omega*n*tau);

        // FD step B
        double B_old = B[M-2];
        B.head(M-1) = B.head(M-1) - (tau/h) * (E.tail(M-1) - E.head(M-1));
        
        // Boundary conditions for B field (right)
        B[M-1] = B_old+(tau-h)/(tau+h)*(B[M-2]-B[M-1]); // absorbing
        //B[M-1] = 0;                                     // reflective        
        //B[M-1] = B[0];                                  // periodic        

        // save the current timestep
        E_arr.row(n) = E;
        B_arr.row(n) = B;
        cout << "\rSimulating... " << n+1 << "/" << N;
    }
    cout << "\n";

    cout << "Writing files...\n";
    ofstream file("data/output_E.csv", ios::trunc);
    if (file.is_open())
    {
        file << "# E field (N rows, M columns)\n";
        file << E_arr << '\n';
    }
    file.close();

    file = ofstream("data/output_B.csv", ios::trunc);
    if (file.is_open())
    {
        file << "# B field\n";
        file << B_arr << '\n';
    }
    file.close();

    // print metadata for plotting
    file = ofstream("data/output_constants.json", ios::trunc);
    if (file.is_open())
    {
        file << "{\n";
        file << "\"M\": " << M << ",\n";
        file << "\"x_min\": " << x_min << ",\n";
        file << "\"x_max\": " << x_max << ",\n";
        file << "\"dx\": " << h << ",\n";
        file << "\"N\": " << N << ",\n";
        file << "\"t_min\": " << t_min << ",\n";
        file << "\"t_max\": " << t_max << ",\n";
        file << "\"dt\": " << tau << "\n";
        file << "}\n";
    }
    file.close();

    cout << "Done.\n";
}
