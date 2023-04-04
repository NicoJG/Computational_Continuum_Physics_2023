#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

int main()
{
    // constants of the simulation
    const double c = 0.1;
    const double k = M_PI;
    const double omega_0 = sqrt(2);
    const double omega_1 = sqrt(3);

    const double x_min = 0;
    const double x_max = 1;
    const double t_min = 0;
    const double t_max = 10;

    const double sigma = 0.45;
    const double dx = 0.005;
    const double dt = sigma*dx*dx/c;

    const int N = (int) ((t_max - t_min) / dt) + 1; // number of timesteps
    const int M = (int) ((x_max - x_min) / dx) + 1; // number of spacesteps

    cout << "dx = " << dx << " dt = " << dt << "\n";
    cout << "M = " << M << " N = " << N << "\n";
    cout << "sigma = " << sigma << "\n";

    // check if stability condition is satisfied
    if (sigma>0.5) {
        cout << "WARNING: Instable FDM\n";
    }

    // create matrix with each row representing one timestep
    auto u = new double[N][M];

    // initial conditions
    for (int m=0; m<M; m++) {
        u[0][m] = sin(k*(x_min + m*dx))*sin(k*(x_min + m*dx));
    }

    // loop through the timesteps
    cout << "\rSimulating... " << 1 << "/" << N;
    for (int n=1; n<N; n++) {
        // implement the FTCS scheme:
        for (int m=1; m<M-1; m++) {
            u[n][m] = u[n-1][m] + sigma*(u[n-1][m+1] - 2*u[n-1][m] + u[n-1][m-1]);
        }

        // boundary conditions:
        u[n][0] = u[n][1] - dx*sin(omega_0*n*dt);
        u[n][M-1] = u[n][M-2] + dx*sin(omega_1*n*dt);

        cout << "\rSimulating... " << n+1 << "/" << N;
    }
    cout << "\n";

    // write simulation to a file
    cout << "Writing file...\n";
    ofstream file = ofstream("data/output.csv", ios::trunc);
    file << "# u values of shape (N,M), with rows as timesteps\n";
    for (int n=0; n<N; n++) {
        for (int m=0; m<M; m++) {
            file << u[n][m] << ", ";
        }
        file << "\b\b\n";
    }
    file.close();
    cout << "Done.\n";


    // find u(x=0.3, t=10.0)
    int n = N-1;
    int m = round((0.3-x_min)/dx);

    double x = x_min + m*dx;
    double t = t_min + n*dt;

    printf("u[%i][%i] = u(x = %.2f, t = %.1f) = %.5f \n", n,m,x,t, u[n][m]);
}