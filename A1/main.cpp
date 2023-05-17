#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
 
int main()
{

    double sigma = 1.0;  // Courant number 

    // Space discretization
    double x_min = -1.0;
    double x_max = 1.0;
    const int M = 500; // Number of grid points for E
    double h = (x_max - x_min) / M;

    // Time discretization
    double t_min = 0.0;
    double t_max = 2.5;
    double tau   = sigma * h; // Magic time step
    const int N  = ceil((t_max - t_min) / tau); // number of timesteps


    // For storing snapshots of fields
    int frame_count  = 25;
    int frame_stride = N / frame_count;

    // check if stability condition (sigma = tau/h < 1) is true
    if (tau/h > 1) {
        cerr << "ERROR! Method is instable: tau = " << tau << " h = " << h << "\n";
        exit(1);
    }

    // For storing fields
    double E[M], B[M-1], E_old[M], B_old[M-1];

    // linspaces for the x values
    double x_E[M], x_B[M-1]; 
    for (int m=0; m<M; m++) {
        x_E[m] = x_min + m*h;
        x_B[m] = x_min + (m+0.5)*h;
    }

    // initial conditions:
    for (int m=0; m<M; m++) {

        // Calculate the shift in the initial condition for the B field for the first halfstep.
        // We are simulating the wave equation with an initial state traveling to the right
        // (E, B > 0 --> k = E x B > 0), so the shift is given by:
        double c = 1.0;
        double x_shift = c * tau / 2;

        // Gaussian wave packet:
        //E[m] = abs(x_E[m])<0.3 ? exp(-x_E[m]*x_E[m]/(2*0.1*0.1))*sin(20 * M_PI * x_E[m]) : 0;
        //B[m] = abs(x_B[m]+x_shift)<0.3 ? exp(-x_B[m]*(x_B[m]+x_shift)/(2*0.1*0.1))*sin(20 * M_PI * (x_B[m]+x_shift)) : 0;
        
        // Sine wave packet
        //E[m] = abs(x_E[m])<0.1 ? sin(20 * M_PI * x_E[m]) : 0;
        //B[m] = abs(x_B[m]+x_shift)<0.1 ? sin(20 * M_PI * (x_B[m]+x_shift)) : 0;
        
        // Nothing
        E[m] = 0;
        B[m] = 0;
    }

    // write header and initial conditions to the output file
    ofstream file = ofstream("data/output.csv", ios::trunc);
    file << "# n, t, E[0:(M-1)], B[0:(M-1)]\n";
    file << 0 << ", " << t_min;
    for (int m=0; m<M; m++) {
        file << ", " << E[m];
    }
    for (int m=0; m<M-1; m++) {
        file << ", " << B[m];
    }
    file << "\n";
    file.close();
    std::cout << "\rSimulating... " << 1 << "/" << N;

    for (int n=1; n<N; n++) {

        // copy arrays
        for (int m=0; m<M; m++)   E_old[m] = E[m];
        for (int m=0; m<M-1; m++) B_old[m] = B[m];

        // FD step B
        for (int m=0; m<M-1; m++)
            B[m] = B_old[m] - (tau/h) * (E_old[m+1] - E_old[m]);

        // FD step E
        for (int m=1; m<M-1; m++)
            E[m] = E_old[m] - (tau/h) * (B[m] - B[m-1]);
        
        // Boundary conditions on the left
        //E[0] = E_old[1] + (tau-h)/(tau + h) * (E[1]-E[0]);    // absorbing
        E[0] = 0;                                             // reflective
        //E[0] = E_old[0] - (tau/h) * (B[0] - B[M-2]);          // periodic
        //E[0] = 0.5 * sin(10 * n * tau);                       // antenna

        // Boundary conditions on the right
        //E[M-1] = E_old[M-2] + (tau-h)/(tau+h) * (E[M-2]-E[M-1]);  // absorbing
        E[M-1] = 0;                                               // reflective        
        //E[M-1] = E[0];                                            // periodic
        //E[M-1] = 0.5 * sin(10 * n * tau);                         // antenna     

        // antenna in the middle
        //E[M/2] = sin(20 * n * tau);

        // save the current timestep in regular intervals
        if (((n+1) % frame_stride == 0) || (n==(N-1))) {   
            file = ofstream("data/output.csv", ios::app);
            file << n << ", " << t_min+n*tau;
            for (int m=0; m<M; m++) {
                file << ", " << E[m];
            }
            for (int m=0; m<M-1; m++) {
                file << ", " << B[m];
            }
            file << "\n";
            file.close();
            std::cout << "\rSimulating... " << n+1 << "/" << N;
        }
    }
    std::cout << "\n";

    // write metadata to a file for plotting 
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

    std::cout << "Done.\n";
}
