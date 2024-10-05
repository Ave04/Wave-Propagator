#ifndef ShallowWater_H
#define ShallowWater_H

class ShallowWater {
    
public:
    ShallowWater(int _Nx, int _Ny, double _dt, double _T, int _ic, double* _u, double* _v, double* _h); // constructor
    
    // member functions
    void SpatialDerivativeY(double* u, double* v,double* h,double* du_dy,double* dv_dy,double* dh_dy, double dy, int Nx, int Ny);
    void SpatialDerivativeX(double* u, double* v,double* h,double* du_dx,double* dv_dx,double* dh_dx, double dx, int Nx, int Ny);
    void PrintMatrix(int Nx, int Ny, double* _u);
    void SetInitialConditions(int ic, int Ny, int Nx, double dx, double dy, double* h, double* u, double* v);
    void TimeIntegrate(double dt, double T, double* u, double* v, double* h, int Nx, int Ny, double g, double dx, double dy);
    void WriteOutput(double* u, double* v, double* h, double dx, double dy);
    void RHS(double* R1, double* R2, double* R3, double* u, double* v, double* h, double* dudx, double* dudy, double* dhdx, double* dvdx, double* dvdy, double* dhdy, double g);
    void BLASversionRHS(double* R1, double* R2, double* R3, double* u, double* v, double* h, double* dudx, double* dudy, double* dhdx, double* dvdx, double* dvdy, double* dhdy, double g);
;

    // PARAMETERS
    int Nx;
    int Ny;
    double dt;
    double T;
    int ic;
    
    double g = 9.81; // acceleration due to gravity
    double Lx = 100.0, Ly = 100.0; // dimensions of domain
    double dx = Lx/Nx, dy = Ly/Ny; // grid spacing
 
    //solution fields
    double* u ; 
    double* v ;
    double* h ;

};

#endif