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

    const double sigma = 0.5;
    const double dx = 0.0005;
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
    auto u = new double[M];
    auto u_next = new double[M];

    // initial conditions
    for (int m=0; m<M; m++) {
        u[m] = sin(k*(x_min + m*dx))*sin(k*(x_min + m*dx));
    }

    // loop through the timesteps
    cout << "\rSimulating... " << 1 << "/" << N;
    for (int n=1; n<N; n++) {
        // implement the FTCS scheme:
        for (int m=1; m<M-1; m++) {
            u_next[m] = u[m] + sigma*(u[m+1] - 2*u[m] + u[m-1]);
        }

        // boundary conditions:
        u_next[0] = u[1] - dx*sin(omega_0*n*dt);
        u_next[M-1] = u[M-2] + dx*sin(omega_1*n*dt);

        for (int m=0; m<M; m++) {
            u[m] = u_next[m];
        }

        if (n % 1000 == 0) {
            cout << "\rSimulating... " << n+1 << "/" << N;
        }
    }
    cout << "\n";

    // find u(x=0.3, t=10.0)
    int n = N-1;
    int m = floor((0.3-x_min)/dx);
    double x = x_min + m*dx;
    double t = t_min + n*dt;
    printf("u[%i][%i] = u(x = %.4f, t = %.1f) = %.5f \n", n,m,x,t, u[m]);

    n = N-1;
    m = ceil((0.3-x_min)/dx);
    x = x_min + m*dx;
    t = t_min + n*dt;
    printf("u[%i][%i] = u(x = %.4f, t = %.1f) = %.5f \n", n,m,x,t, u[m]);
}