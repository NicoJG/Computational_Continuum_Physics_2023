#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <stdio.h>

using namespace std;

void pyPlot2D(double** X, int Nx, int Ny, string fileName = "", string kwarg = "", bool noFrame = false) // X[iy][ix]
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
    kwarg += ", vmin = -vAbsMax, vmax = vAbsMax";

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

void test_plot2D_dots()
{
    cout << "test_plot2D_dots" << endl;
    int Nx = 128, Ny = 64;
    double **X = new double*[Ny]; for(int i = 0; i < Ny; i++) X[i] = new double[Nx];
    for(int iy = 0; iy < Ny; iy++)
    for(int ix = 0; ix < Nx; ix++){
        double y = iy/double(Ny), x = ix/double(Nx);
        X[iy][ix] = sin(2*M_PI*y) + 0.0*sin(2*M_PI*(x + 2*y));
    }
    vector<vector<double>> dots;
    dots.push_back({0.25, 0.5});
    dots.push_back({0.27, 0.5});
    pyPlot2D(X, Nx, Ny, "test_plot2d", ", cmap='RdBu'");
};