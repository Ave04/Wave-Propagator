#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>

#include "ShallowWater.h"

#include "cblas.h"
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]) {
    
    // Parse command line arguments
    // define parameters
    double dt;
    double T;
    int Nx;
    int Ny;
    int ic;
    
    // read from command line and store in variable
    for (int i = 1; i < argc; i +=2) {
        
        string arg  = argv[i];
        
        if (arg == "--dt") {dt = atof(argv[i+1]);}
        else if (arg == "--T") {T = atof(argv[i+1]);}
        else if (arg == "--Nx") {Nx = atoi(argv[i+1]);}
        else if (arg == "--Ny") {Ny = atoi(argv[i+1]);}
        else if (arg == "--ic") {ic = atoi(argv[i+1]);}
        
        // Invalid condition
        else {
            cout << "invalid input: " << arg << endl;
            return 1;
        }
        
    }
    
    // code tester
   // dt = 0.1;
   // T = 80;
   // Nx = 100;
   // Ny = 100;
   // ic = 3;
    
    // initialise solution fields   
    double* u = new double[Nx*Ny];
    double* v = new double[Nx*Ny];
    double* h = new double[Nx*Ny];
    
    // display to test
    cout << "dt: " << dt << endl;
    cout << "T: " << T << endl;
    cout << "Nx: " << Nx << endl;
    cout << "Ny: " << Ny << endl;
    cout << "ic: " << ic << endl;

    // Initialise ShallowWater class
    ShallowWater sw(Nx, Ny, dt, T, ic, u, v, h); 
    
    // Set Initial conditions
    sw.SetInitialConditions(sw.ic, sw.Ny, sw.Nx, sw.dx, sw.dy, sw.h, sw.u, sw.v); 
    
    // Time integrate for T time domain
    sw.TimeIntegrate(dt, T, sw.u, sw.v, sw.h, sw.Nx, sw.Ny, sw.g, sw.dx, sw.dy); 

    // Write output to file
    sw.WriteOutput(sw.u, sw.v, sw.h, sw.dx, sw.dy); 
    
    // Deallocate dynamic memory
    delete[] u;
	delete[] v;
	delete[] h;
    
    return 0;
}