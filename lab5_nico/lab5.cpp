#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <stdio.h>

using namespace std;

void pyPlot2D(double** X, int Nx, int Ny, string fileName = "", string kwarg = "", bool noFrame = false, bool positive = false) // X[iy][ix]
{
    ofstream file; file.open("plot2D_tmp.dat", ios::out); file.close();
    FILE * dataFile = fopen ("plot2D_tmp.dat", "ab");
    fwrite (&Ny, sizeof(int), 1, dataFile);
    fwrite (&Nx, sizeof(int), 1, dataFile);
    for(int y = 0; y < Ny; y++) for(int x = 0; x < Nx; x++) fwrite (&(X[y][x]), sizeof(double), 1, dataFile);
    fclose (dataFile);
    ofstream pyFile; pyFile.open("plot2D_tmp.py", ios::out); 
    pyFile << "import matplotlib.pyplot as plt" << endl;
    pyFile << "import numpy as np" << endl; 
    pyFile << "import struct" << endl;
    pyFile << "f = open('plot2D_tmp.dat', 'rb')" << endl;
    pyFile << "Ny = struct.unpack('I', f.read(4))[0]" << endl;
    pyFile << "Nx = struct.unpack('I', f.read(4))[0]" << endl;
    pyFile << "v = np.empty(shape=(Ny, Nx))" << endl;
    pyFile << "for y in range(Ny):" << endl;
    pyFile << "    for x in range(Nx):" << endl;
    pyFile << "        v[y][x] = struct.unpack('d', f.read(8))[0]" << endl;
    pyFile << "fig, ax = plt.subplots()" << endl;
    
    pyFile << "vAbsMax = np.maximum(abs(np.amin(v)), abs(np.amax(v)))" << endl;
    if(!positive) kwarg += ", vmin = -vAbsMax, vmax = vAbsMax";
        else kwarg += ", vmin = 0, vmax = vAbsMax";

    pyFile << "plot = ax.imshow(v, interpolation='none', extent=[0,1,0,1], origin='lower'" << kwarg << ")" << endl;
    pyFile << "ax.set_aspect('equal')" << endl;
    pyFile << "fig.colorbar(plot, ax=ax, location='right')" << endl;
    //pyFile << "ax.contour(v, np.array([0]), colors='k', extent=[0,1,0,1])" << endl;
    //if(fileName == "") pyFile << "plt.show()" << endl;
    //    else pyFile << "plt.savefig('" << fileName << ".png')" << endl;
    if(noFrame){
        pyFile << "plt.gca().set_axis_off()" << endl;
        pyFile << "plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)" << endl;
        pyFile << "plt.margins(0,0)" << endl;
        pyFile << "plt.gca().xaxis.set_major_locator(plt.NullLocator())" << endl;
        pyFile << "plt.gca().yaxis.set_major_locator(plt.NullLocator())" << endl;
        pyFile << "plt.savefig('" << fileName << ".png', bbox_inches = 'tight', pad_inches = 0)" << endl;
    } else 
        pyFile << "plt.savefig('" << fileName << ".png')" << endl;
    
    pyFile.close();
    system("python plot2D_tmp.py");
};

void test_plot2D()
{
    cout << "test_plot2D" << endl;
    int Nx = 128, Ny = 64;
    double **X = new double*[Ny]; 
    for(int i = 0; i < Ny; i++) X[i] = new double[Nx];
    for(int iy = 0; iy < Ny; iy++)
    for(int ix = 0; ix < Nx; ix++){
        double y = iy/double(Ny), x = ix/double(Nx);
        X[iy][ix] = sin(2*M_PI*y) + 0.0*sin(2*M_PI*(x + 2*y));
    }
    pyPlot2D(X, Nx, Ny, "figures/test_plot2d", ", cmap='RdBu'");
    for(int i = 0; i < Ny; i++) delete[](X[i]);
    delete[](X);
};

void plot1D(vector<double> X, string fileName = "", double yRange = 0) // function for plotting a vector<double>
{
    ofstream file; file.open("plot1D_tmp.dat", ios::out); file.close();
    FILE * dataFile = fopen ("plot1D_tmp.dat", "ab");
    int N = X.size(); fwrite (&N, sizeof(int), 1, dataFile);
    for(int i = 0; i < X.size(); i++) fwrite (&(X[i]), sizeof(double), 1, dataFile);
    fclose (dataFile);
    ofstream pyFile; pyFile.open("plot1D_tmp.py", ios::out); 
    pyFile << "import matplotlib.pyplot as plt" << endl;
    pyFile << "import struct" << endl;
    pyFile << "f = open('plot1D_tmp.dat', 'rb')" << endl;
    pyFile << "N = struct.unpack('I', f.read(4))[0]" << endl;
    pyFile << "y = [0]*N" << endl;
    pyFile << "for i in range(N):" << endl;
    pyFile << "    y[i] = struct.unpack('d', f.read(8))[0]" << endl;
    pyFile << "fig, ax = plt.subplots()" << endl;
    pyFile << "ax.plot(y)" << endl;
    if(yRange != 0) pyFile << "ax.set_ylim(" << -yRange << ", "<< yRange << ")" << endl;
    if(fileName == "") pyFile << "plt.show()" << endl;
        else pyFile << "plt.savefig('" << fileName << ".png')" << endl;
    pyFile.close();
    system("python plot1D_tmp.py");
};
double randNormal(double sigma2) // generation of a pseudorandom number following normal distribution.
{ 
    return cos(2*M_PI*(rand()/double(RAND_MAX)))*sqrt(-2*sigma2*log(1 - (1 + rand())/double((long)RAND_MAX + 2)));
}
//CGS units are used
const double electronCharge = -4.80320427e-10;
const double electronMass = 9.10938215e-28;
const double lightVelocity = 29979245800.0;

const double k_B = 1.380649e-16;
double sqr(double x){return x*x;};

// Simulation constants
const double T = 1e-6*electronMass*lightVelocity*lightVelocity;
const double N_e = 1e18;
const double lambda_D = sqrt(T/(4*M_PI*N_e*electronCharge*electronCharge));
const double T_p = sqrt(M_PI*electronMass/(N_e*electronCharge*electronCharge));


const double L = 100*lambda_D;
const int time_steps_per_T_p = 1e3;
const double dt = T_p/time_steps_per_T_p;
const double dx = lambda_D;
const int particles_per_cell = 1e3;
const int nCells = ceil(L/dx);
const int nParticles = nCells*particles_per_cell;

const double a = 0.1;

struct electrostaticPIC1D{
    vector<double> Ex, Jx;
    vector<double> px, x;
    double L;
    electrostaticPIC1D(double L, double dt, int nParticles, int nCells): L(L)
    {
        Ex.resize(nCells);
        Jx.resize(nCells);
        px.resize(nParticles);
        x.resize(nParticles);

        double sigma2 = 2*T*electronMass;

        for(int i = 0; i < px.size();i++)
        {

            x[i] = L*(rand() + 1)/double((long)(RAND_MAX) + 2);
            px[i] = randNormal(sigma2); // here we need to generate a random value to get the normal distribution with a given temperature; fucntion randNormal() can be helpful 
        }
    }
    void advance()// this function is to advance the state of both particles and fields
    {
        //compute Jx
        for(int ix = 0; ix < Ex.size(); ix++) {
            Jx[ix] = 0;
        }
        
        double weight = N_e*L/nParticles;
        for(int ip = 0; ip < px.size(); ip++){
            // here we need to deposite the contribution of ip-th particle to the grid value of Jx at its location
            int ix = floor(x[ip]/dx);
            Jx[ix] += px[ip]*electronCharge/electronMass * weight / dx;
        }

        //advance Ex
        for(int ix = 0; ix < Ex.size(); ix++){
            // here we need to do one step for Ex[ix] to advance its state in time
            Ex[ix] += -4*M_PI*Jx[ix]*dt;
        }

        //advance particles
        for(int ip = 0; ip < px.size(); ip++){
            //here we need to advance the state of particles
            x[ip] += px[ip]/electronMass*dt;
            // periodic boundary conditions
            while (x[ip] >= L) x[ip] -= L;
            while (x[ip] < 0) x[ip] += L;

            int ix = floor(x[ip]/dx);
            px[ip] += electronCharge*Ex[ix]*dt;
        }
    }
    void plot_xpx(double temperature, int imageNumber)
    {
        int nx = 100;
        int npx = 100;
        double pxMax = 30*sqrt(temperature*electronMass);
        double **d = new double*[npx];
        for(int i = 0; i < npx; i++){
            d[i] = new double[nx];
            for(int ix = 0; ix < nx; ix++)d[i][ix] = 0;
        }

        for(int ip = 0; ip < px.size(); ip++){
            int ix = int(nx*x[ip]/L + 0.5) % nx;
            int ipx = int(npx*(0.5 + px[ip]/pxMax));
            if((ipx > 0)&&(ipx < npx)){
                d[ipx][ix] += 1;
            }
        }
        char filename[100];
        sprintf(filename, "figures/xpx_%05d", imageNumber);
        pyPlot2D(d, nx, npx, filename, ", cmap='plasma'", false, true);
        for(int i = 0; i < npx; i++) delete []d[i];
        delete []d;
    }
};

int main()
{
    cout << "lab5" << endl;
    // here we need to create the instance of sruct "electrostaticPIC1D", set up initial conditions and call advance() to advance its state.
    electrostaticPIC1D pic(L,dt,nParticles,nCells);

    // initial conditions:
    double E_amplitude = a * 4*M_PI*L*electronCharge*N_e;
    for(int ix = 0; ix < pic.Ex.size(); ix++){
        pic.Ex[ix] = E_amplitude * sin(2*M_PI*ix*dx/L);
    }

    char filename[100];
    sprintf(filename, "figures/Ex_%05d",0);
    plot1D(pic.Ex, filename, abs(E_amplitude));
    pic.plot_xpx(T,0);

    int Nt = 10*time_steps_per_T_p;
    int frames = 20;
    int frame_stride = Nt/frames;

    double initial_energy = 0;
    for (int ix = 0; ix < pic.Ex.size(); ix++) {
        initial_energy += pic.Ex[ix]*pic.Ex[ix];
    }

    cout << "iteration " << 0 << "/" << Nt;
    for (int it=0; it < Nt; it++) {
        pic.advance();
        cout << "\riteration " << it+1 << "/" << Nt << flush;
        if ((it == Nt-1) || (it % frame_stride == 0)) {
            sprintf(filename, "figures/Ex_%05d",it+1);
            plot1D(pic.Ex, filename, abs(E_amplitude));
            pic.plot_xpx(T,it+1);
        }
    }
    cout << endl;

    
    double final_energy = 0;
    for (int ix = 0; ix < pic.Ex.size(); ix++) {
        final_energy += pic.Ex[ix]*pic.Ex[ix];
    }

    cout << "field energy ratio: " << initial_energy/final_energy << endl;

}