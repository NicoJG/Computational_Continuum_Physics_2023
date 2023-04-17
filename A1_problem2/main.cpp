#include <stdio.h>
#include <math.h>

using namespace std;
 
int main() {

    // --- SIMULATION PARAMETERS --- //

    double c     = 1.0; // Propagation speed
    double sigma = 1.0; // Courant number
    double T     = 2.0; // Final time
    double L     = 2.0; // Width of the centered simulation box

    // Stability condition for sigma = 1: 
    // 0.25 <= xi_x + xi_y <= 0.5
    double xi_x = 0.125;
    double xi_y = 0.125;


    // Space step
    const int M = 100;
    double h = L / M;

    // Time step
    double tau = 1 * sqrt(sigma) * h / c; // sigma = c2 tau2 / h2
    int N      = ceil(T / tau);

    // For saving frames
    int frames = 20;
    int frame_stride = N/frames;


    // --- GRID SETUP --- //

    //   -   -   -        -  E_x
    // | o | o | o |      |  E_y
    //   -   -   -        o  B_z
    // | o | o | o |      
    //   -   -   -        
    // | o | o | o | 
    //   -   -   -  
    //
    // |<--------->| x domain

    double E_x[M][M+1],     E_y[M+1][M],     B_z[M][M];
    double E_x_old[M][M+1], E_y_old[M+1][M], B_z_old[M][M];

    // --- INITIAL CONDITIONS --- //

    // E_x
    for (int m_x = 0; m_x < M;   m_x++)
    for (int m_y = 0; m_y < M+1; m_y++) {
        double x = (m_x + 0.5) * h - L/2;
        double y = (m_y      ) * h - L/2;
        E_x[m_x][m_y] = 0 * exp(-(x*x + y*y)/(0.5 * 0.01));
    }

    // E_y
    for (int m_x = 0; m_x < M+1; m_x++)
    for (int m_y = 0; m_y < M  ; m_y++) {
        double x = (m_x      ) * h - L/2;
        double y = (m_y + 0.5) * h - L/2;
        E_y[m_x][m_y] = 0 * sin(2 * 3.14 * 5 * x) * exp(-(x*x + y*y)/(2 * 0.01));
    }

    // B_z
    for (int m_x = 0; m_x < M; m_x++)
    for (int m_y = 0; m_y < M; m_y++) {
        double x = (m_x + 0.5) * h - L/2;
        double y = (m_y + 0.5) * h - L/2;
        B_z[m_x][m_y] = 0 * sin(2 * 3.14 * 5 * x) * exp(-(x*x + y*y)/(2 * 0.01));
    }


    // --- FINITE DIFFERENCE METHOD --- //

    for (int n = 0; n < N+1; n++) {
        
        // SAVE TO FILE 

        if (n % frame_stride == 0) {
            char filename[100];
            sprintf(filename, "data/B_field_t=%.2f.csv", n * tau);
            FILE* f = fopen(filename, "w");
            for (int m_y = 0; m_y < M; m_y++) {
                for (int m_x = 0; m_x < M; m_x++) {
                    if (m_x != M-1) fprintf(f, "%.5e,",  B_z[m_x][m_y]);
                    else            fprintf(f, "%.5e\n", B_z[m_x][m_y]);
                }
            }  
            fclose(f);
        }

        if (n == 0) continue;

        // COPY FIELDS

        // E_x
        for (int m_x = 0; m_x < M;   m_x++)
        for (int m_y = 0; m_y < M+1; m_y++)
            E_x_old[m_x][m_y] = E_x[m_x][m_y];

        // E_y
        for (int m_x = 0; m_x < M+1; m_x++)
        for (int m_y = 0; m_y < M  ; m_y++)
            E_y_old[m_x][m_y] = E_y[m_x][m_y];

        // B_z
        for (int m_x = 0; m_x < M; m_x++)
        for (int m_y = 0; m_y < M; m_y++)
            B_z_old[m_x][m_y] = B_z[m_x][m_y];


        // FD STEP

        // dE_x/dt = c dB_z/dy
        for (int m_x = 0; m_x < M; m_x++) {

            // Inside domain 
            for (int m_y = 1; m_y < M; m_y++)
                E_x[m_x][m_y] = 
                    E_x_old[m_x][m_y] + (tau*c/h) * 
                    (B_z[m_x][m_y] - B_z[m_x][m_y-1]);

            
            // On boundaries
            E_x[m_x][0] = 0;
            E_x[m_x][M] = 0;
        }

        // dE_y/dt = -c *dB_z/dz
        for (int m_y = 0; m_y < M; m_y++) {

            // Inside domain
            for (int m_x = 1; m_x < M  ; m_x++)
                E_y[m_x][m_y] = 
                    E_y_old[m_x][m_y] - (tau*c/h) * 
                    (B_z[m_x][m_y] - B_z[m_x-1][m_y]);
            
            // On boundaries
            E_y[0][m_y] = (m_y < M/2) ? sin(2 * 3.14159 * 14 * n * tau) : 
                                        sin(2 * 3.14159 * 7  * n * tau); // antenna on left edge
            E_y[M][m_y] = 0;
        }
        
        // Wall in the middle
        for (int m_x = 0; m_x < M; m_x++) {
            E_x[m_x][M/2] = 0;
            E_x[m_x][M/2+1] = 0;
        }

        // Wall in the middle
        // for (int m_x = 0; m_x < M; m_x++)
        //    E_y[m_x][M/2] = 0;
            
        // dB_z/dt = c (dE_x/dy - dE_y/dx)
        // Using special stencil!
        for (int m_x = 0; m_x < M; m_x++)
        for (int m_y = 0; m_y < M; m_y++) {
            
            if (m_x > 0  ) B_z[m_x][m_y] += (c*tau/h) * xi_x       * (E_x[m_x-1][m_y+1] - E_x[m_x-1][m_y]);
                           B_z[m_x][m_y] += (c*tau/h) * (1-2*xi_x) * (E_x[m_x  ][m_y+1] - E_x[m_x  ][m_y]);
            if (m_x < M-1) B_z[m_x][m_y] += (c*tau/h) * xi_x       * (E_x[m_x+1][m_y+1] - E_x[m_x+1][m_y]);

            if (m_y > 0)   B_z[m_x][m_y] -= (c*tau/h) * xi_y       * (E_y[m_x+1][m_y-1] - E_y[m_x][m_y-1]);
                           B_z[m_x][m_y] -= (c*tau/h) * (1-2*xi_y) * (E_y[m_x+1][m_y  ] - E_y[m_x][m_y  ]);
            if (m_y < M-1) B_z[m_x][m_y] -= (c*tau/h) * xi_y       * (E_y[m_x+1][m_y+1] - E_y[m_x][m_y+1]);

        }

        // Wall in the middle
        // for (int m_x = 0; m_x < M; m_x++)
        //    B_z[m_x][M/2] = 0;

    }

    return 0;
}
