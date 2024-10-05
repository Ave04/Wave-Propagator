#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <cstdlib>
#include "ShallowWater.h"

#include "cblas.h"
#include <omp.h>


using namespace std;

// constructor
ShallowWater::ShallowWater(int _Nx, int _Ny, double _dt, double _T, int _ic, double* _u, double* _v, double* _h)
    :Nx(_Nx), Ny(_Ny), dt(_dt), T(_T), ic(_ic), u(_u), v(_v), h(_h)
{}

// Print Matrix
void ShallowWater::PrintMatrix(int Nx, int Ny, double* u) {
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int k = j*Nx+i;
            cout << setprecision(4) << u[k] << " ";
        }
        cout << endl;
    }
}

// Function to calculate RHS of equation ( f functions )
void ShallowWater::RHS(double* R1, double* R2, double* R3, double* u, double* v, double* h, double* dudx, double* dudy, double* dhdx, double* dvdx, double* dvdy, double* dhdy, double g) {
    #pragma omp parallel for  
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int k = j*Nx+i;
            R1[k] = -1*(u[k]*dudx[k] + v[k]*dudy[k] + g*dhdx[k]);
            R2[k] = -1*(u[k]*dvdx[k] + v[k]*dvdy[k] + g*dhdy[k]);
            R3[k] = -1*(u[k]*dhdx[k] + v[k]*dhdy[k] + h[k]*dudx[k] + h[k]*dvdy[k]);
        }
    }
}

// FDM STENCIL IMPLEMENTATION

// x spatial derivatives
void  ShallowWater::SpatialDerivativeX(double* u, double* v,double* h, double* dudx, double* dvdx,double* dhdx, double dx, int Nx, int Ny) {
    #pragma omp parallel for
    for ( int j = 0; j<Ny; j++) {
			for ( int i = 0; i<Nx; i++){
            
				int k = j*Nx+i ;
				if (i >2){
					dudx[k] = 1/dx*(-1.0/60*u[k-3] 
								+ 3.0/20*u[k-2] 
								- 3.0/4*u[k-1] 
								+ 3.0/4*u[j*Nx +(k+1)%Nx] 
								- 3.0/20*u[j*Nx +(k+2)%Nx] 
								+ 1.0/60*u[j*Nx +(k+3)%Nx]);
							  
					dvdx[k] = 1/dx*(-1.0/60*v[k-3] 
								+ 3.0/20*v[k-2] 
								- 3.0/4*v[k-1] 
								+ 3.0/4*v[j*Nx +(k+1)%Nx] 
								- 3.0/20*v[j*Nx +(k+2)%Nx] 
								+ 1.0/60*v[j*Nx +(k+3)%Nx]);
							  
					dhdx[k] = 1/dx*(-1.0/60*h[k-3] 
								+ 3.0/20*h[k-2] 
								- 3.0/4*h[k-1] 
								+ 3.0/4*h[j*Nx +(k+1)%Nx] 
								- 3.0/20*h[j*Nx +(k+2)%Nx] 
								+ 1.0/60*h[j*Nx +(k+3)%Nx]);
					
                }
                    
                if (i==2){
					dudx[k] = 1/dx*(-1.0/60*u[Nx*(j+1)-1] 
								+ 3.0/20*u[k-2] 
								- 3.0/4*u[k-1] 
								+ 3.0/4*u[j*Nx +(k+1)%Nx] 
								- 3.0/20*u[j*Nx +(k+2)%Nx] 
								+ 1.0/60*u[j*Nx +(k+3)%Nx]);
							  
					dvdx[k] = 1/dx*(-1.0/60*v[Nx*(j+1)-1]  
								+ 3.0/20*v[k-2] 
								- 3.0/4*v[k-1] 
								+ 3.0/4*v[j*Nx +(k+1)%Nx] 
								- 3.0/20*v[j*Nx +(k+2)%Nx] 
								+ 1.0/60*v[j*Nx +(k+3)%Nx]);
							  
					dhdx[k] = 1/dx*(-1.0/60*h[Nx*(j+1)-1] 
								+ 3.0/20*h[k-2] 
								- 3.0/4*h[k-1] 
								+ 3.0/4*h[j*Nx +(k+1)%Nx] 
								- 3.0/20*h[j*Nx +(k+2)%Nx] 
								+ 1.0/60*h[j*Nx +(k+3)%Nx]); 
                }
				 
                if (i ==1){
					dudx[k] = 1/dx*(-1.0/60*u[Nx*(j+1)-2]
								+ 3.0/20*u[Nx*(j+1)-1] 
								- 3.0/4*u[k-1] 
								+ 3.0/4*u[j*Nx +(k+1)%Nx] 
								- 3.0/20*u[j*Nx +(k+2)%Nx] 
								+ 1.0/60*u[j*Nx +(k+3)%Nx]);
							  
					dvdx[k] = 1/dx*(-1.0/60*v[Nx*(j+1)-2]
								+ 3.0/20*v[Nx*(j+1)-1] 
								- 3.0/4*v[k-1] 
								+ 3.0/4*v[j*Nx +(k+1)%Nx] 
								- 3.0/20*v[j*Nx +(k+2)%Nx] 
								+ 1.0/60*v[j*Nx +(k+3)%Nx]);
							  
					dhdx[k] = 1/dx*(-1.0/60*h[Nx*(j+1)-2]
								+ 3.0/20*h[Nx*(j+1)-1] 
								- 3.0/4*h[k-1] 
								+ 3.0/4*h[j*Nx +(k+1)%Nx] 
								- 3.0/20*h[j*Nx +(k+2)%Nx] 
								+ 1.0/60*h[j*Nx +(k+3)%Nx]); 
                }
				 
                if (i==0){
					dudx[k] = 1/dx*(-1.0/60*u[Nx*(j+1)-3] 
								+ 3.0/20*u[Nx*(j+1)-2] 
								- 3.0/4*u[Nx*(j+1)-1] 
								+ 3.0/4*u[j*Nx +(k+1)%Nx] 
								- 3.0/20*u[j*Nx +(k+2)%Nx] 
								+ 1.0/60*u[j*Nx +(k+3)%Nx]);
							  
					dvdx[k] = 1/dx*(-1.0/60*v[Nx*(j+1)-3] 
								+ 3.0/20*v[Nx*(j+1)-2] 
								- 3.0/4*v[Nx*(j+1)-1] 
								+ 3.0/4*v[j*Nx +(k+1)%Nx] 
								- 3.0/20*v[j*Nx +(k+2)%Nx] 
								+ 1.0/60*v[j*Nx +(k+3)%Nx]);
							  
					dhdx[k] =  1/dx*(-1.0/60*h[Nx*(j+1)-3] 
								+ 3.0/20*h[Nx*(j+1)-2] 
								- 3.0/4*h[Nx*(j+1)-1] 
								+ 3.0/4*h[j*Nx +(k+1)%Nx] 
								- 3.0/20*h[j*Nx +(k+2)%Nx] 
								+ 1.0/60*h[j*Nx +(k+3)%Nx]);
                }

		}
	}
}

// y spatial derivatives
void  ShallowWater::SpatialDerivativeY(double* u, double* v,double* h,double* dudy, double* dvdy, double* dhdy, double dy, int Nx, int Ny) {
    #pragma omp parallel for
    for ( int j = 0; j<Ny; j++) {
			for ( int i = 0; i<Nx; i++){
				
				int k = j*Nx+i ;
				if (j>2){
					dudy[k] = -1/dy*(-1.0/60*u[k-3*Nx]
								+ 3.0/20*u[k-2*Nx] 
								- 3.0/4*u[k-1*Nx] 
								+ 3.0/4*u[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*u[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*u[(k+3*Nx)%(Nx*Ny)]);
							
					dvdy[k] = -1/dy*(-1.0/60*v[k-3*Nx]
								+ 3.0/20*v[k-2*Nx] 
								- 3.0/4*v[k-1*Nx] 
								+ 3.0/4*v[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*v[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*v[(k+3*Nx)%(Nx*Ny)]);
							  
					dhdy[k] = -1/dy*(-1.0/60*h[k-3*Nx]
								+ 3.0/20*h[k-2*Nx] 
								- 3.0/4*h[k-1*Nx] 
								+ 3.0/4*h[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*h[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*h[(k+3*Nx)%(Nx*Ny)]);
				}
				 
				 
                if (j==2){
					dudy[k] = -1/dy*(-1.0/60*u[k+ (Ny-j-1)*Nx]
									+ 3.0/20*u[k-2*Nx] 
									- 3.0/4*u[k-1*Nx] 
									+ 3.0/4*u[(k+1*Nx)%(Nx*Ny)] 
									- 3.0/20*u[(k+2*Nx)%(Nx*Ny)]
									+ 1.0/60*u[(k+3*Nx)%(Nx*Ny)]);
							
					dvdy[k] = -1/dy*(-1.0/60*v[k+ (Ny-j-1)*Nx]
									+ 3.0/20*v[k-2*Nx] 
									- 3.0/4*v[k-1*Nx] 
									+ 3.0/4*v[(k+1*Nx)%(Nx*Ny)] 
									- 3.0/20*v[(k+2*Nx)%(Nx*Ny)]
									+ 1.0/60*v[(k+3*Nx)%(Nx*Ny)]);
							  
					dhdy[k] = -1/dy*(-1.0/60*h[k+ (Ny-j-1)*Nx]
									+ 3.0/20*h[k-2*Nx] 
									- 3.0/4*h[k-1*Nx] 
									+ 3.0/4*h[(k+1*Nx)%(Nx*Ny)] 
									- 3.0/20*h[(k+2*Nx)%(Nx*Ny)]
									+ 1.0/60*h[(k+3*Nx)%(Nx*Ny)]);
                }
				 
				 
                if (j==1){
					dudy[k] = -1/dy*(-1.0/60*u[k+ (Ny-j-2)*Nx]
								+ 3.0/20*u[k+ (Ny-j-1)*Nx] 
								- 3.0/4*u[k-1*Nx] 
								+ 3.0/4*u[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*u[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*u[(k+3*Nx)%(Nx*Ny)]);
							
					dvdy[k] = -1/dy*(-1.0/60*v[k+ (Ny-j-2)*Nx]
								+ 3.0/20*v[k+ (Ny-j-1)*Nx] 
								- 3.0/4*v[k-1*Nx] 
								+ 3.0/4*v[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*v[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*v[(k+3*Nx)%(Nx*Ny)]);
							  
					dhdy[k] = -1/dy*(-1.0/60*h[k+ (Ny-j-2)*Nx]
								+ 3.0/20*h[k+ (Ny-j-1)*Nx] 
								- 3.0/4*h[k-1*Nx] 
								+ 3.0/4*h[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*h[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*h[(k+3*Nx)%(Nx*Ny)]);
                }
				 
                if (j==0){
					 dudy[k] = -1/dy*(-1.0/60*u[k+ (Ny-j-3)*Nx]
								+ 3.0/20*u[k+ (Ny-j-2)*Nx]
								- 3.0/4*u[k+ (Ny-j-1)*Nx] 
								+ 3.0/4*u[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*u[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*u[(k+3*Nx)%(Nx*Ny)]);
							
					dvdy[k] = -1/dy*(-1.0/60*v[k+ (Ny-j-3)*Nx]
								+ 3.0/20*v[k+ (Ny-j-2)*Nx] 
								- 3.0/4*v[k+ (Ny-j-1)*Nx] 
								+ 3.0/4*v[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*v[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*v[(k+3*Nx)%(Nx*Ny)]);
							  
					dhdy[k] = -1/dy*(-1.0/60*h[k+ (Ny-j-3)*Nx]
								+ 3.0/20*h[k+ (Ny-j-2)*Nx] 
								- 3.0/4*h[k+ (Ny-j-1)*Nx] 
								+ 3.0/4*h[(k+1*Nx)%(Nx*Ny)] 
								- 3.0/20*h[(k+2*Nx)%(Nx*Ny)]
								+ 1.0/60*h[(k+3*Nx)%(Nx*Ny)]);
                }

			}
	}

}

// TimeIntegration Function
void ShallowWater::TimeIntegrate(double dt, double T, double* u, double* v, double* h, int Nx, int Ny, double g, double dx, double dy) {
    
    // Initialise time
    double t = 0;
    
    while (t <= T) {
        
    
        // derivatives
        double* dudx = new double[Nx*Ny];
        double* dudy = new double[Nx*Ny];
        double* dvdx = new double[Nx*Ny];
        double* dvdy = new double[Nx*Ny];
        double* dhdx = new double[Nx*Ny];
        double* dhdy = new double[Nx*Ny];
        
        // Intiatialise temporary variables
        double* utemp = new double[Nx*Ny];
        double* vtemp = new double[Nx*Ny];
        double* htemp = new double[Nx*Ny];
        
        // Initialise RK4 intermediates
        double* k1u = new double[Nx*Ny];
        double* k2u = new double[Nx*Ny]; 
        double* k3u = new double[Nx*Ny]; 
        double* k4u = new double[Nx*Ny];
        
        double* k1v = new double[Nx*Ny];
        double* k2v = new double[Nx*Ny]; 
        double* k3v = new double[Nx*Ny]; 
        double* k4v = new double[Nx*Ny];
        
        double* k1h = new double[Nx*Ny];
        double* k2h = new double[Nx*Ny]; 
        double* k3h = new double[Nx*Ny]; 
        double* k4h = new double[Nx*Ny];
    
        
        // calculate k intermediates
        SpatialDerivativeX(u, v, h, dudx, dvdx, dhdx,  dx,  Nx,  Ny);
        SpatialDerivativeY(u, v, h, dudy, dvdy, dhdy,  dy,  Nx,  Ny);
        RHS(k1u, k1v, k1h, u, v, h, dudx, dudy, dhdx, dvdx, dvdy, dhdy, g);
        
        #pragma omp parallel for
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int k = j*Nx+i;
                utemp[k] = u[k] + 0.5*dt*k1u[k];
                vtemp[k] = v[k] + 0.5*dt*k1v[k];
                htemp[k] = h[k] + 0.5*dt*k1h[k];
            }
        }
                
        SpatialDerivativeX(utemp, vtemp, htemp, dudx, dvdx, dhdx,  dx,  Nx,  Ny);
        SpatialDerivativeY(utemp, vtemp, htemp, dudy, dvdy, dhdy,  dy,  Nx,  Ny);
        RHS(k2u, k2v, k2h, utemp, vtemp, htemp, dudx, dudy, dhdx, dvdx, dvdy, dhdy, g);
        
        #pragma omp parallel for
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int k = j*Nx+i;
                utemp[k] = u[k] + 0.5*dt*k2u[k];
                vtemp[k] = v[k] + 0.5*dt*k2v[k];
                htemp[k] = h[k] + 0.5*dt*k2h[k];
                
            }
        }
        
        SpatialDerivativeX(utemp, vtemp, htemp, dudx, dvdx, dhdx,  dx,  Nx,  Ny);
        SpatialDerivativeY(utemp, vtemp, htemp, dudy, dvdy, dhdy,  dy,  Nx,  Ny);
        RHS(k3u, k3v, k3h, utemp, vtemp, htemp, dudx, dudy, dhdx, dvdx, dvdy, dhdy, g);
        
        #pragma omp parallel for
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int k = j*Nx+i;
                utemp[k] = u[k] + dt*k3u[k];
                vtemp[k] = v[k] + dt*k3v[k];
                htemp[k] = h[k] + dt*k3h[k];
            }
        }
    
        SpatialDerivativeX(utemp, vtemp, htemp, dudx, dvdx, dhdx,  dx,  Nx,  Ny);
        SpatialDerivativeY(utemp, vtemp, htemp, dudy, dvdy, dhdy,  dy,  Nx,  Ny);
        RHS(k4u, k4v, k4h, utemp, vtemp, htemp, dudx, dudy, dhdx, dvdx, dvdy, dhdy, g);
                
        // fill u, v and h
        #pragma omp parallel for
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int k = j*Nx+i;
                u[k] = u[k] + (k1u[k] + 2*k2u[k] +2*k3u[k] +k4u[k])*dt/6;
                v[k] = v[k] + (k1v[k] + 2*k2v[k] +2*k3v[k] +k4v[k])*dt/6;
                h[k] = h[k] + (k1h[k] + 2*k2h[k] +2*k3h[k] +k4h[k])*dt/6;
            }
        }
        
        // update time
        t = t + dt;
                
    }
    
}

void ShallowWater::SetInitialConditions(int ic, int Ny, int Nx, double dx, double dy, double* h, double* u, double* v) {
    // Set initial conditions for u and v
    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            u[j*Nx+i] = 0;
            v[j*Nx+i] = 0;
        }
    }
    
    // Set initial conditions based on chosen ic for h
    switch (ic) {
        case 1:
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    double x = i * dx;
                    // double y = j * dy;
                    h[j*Nx+i] = 10.0 + exp(-(x - 50.0) * (x - 50.0) / 25.0);
                }
            }
            break;
        case 2:
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    // double x = i * dx;
                    double y = j * dy;
                    h[j*Nx+i] = 10.0 + exp(-(y - 50.0) * (y - 50.0) / 25.0);
                }
            }
            break;
       case 3:
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    double x = i * dx;
                    double y = j * dy;
                    h[j*Nx+i] = 10.0 + exp(-((y - 50.0) * (y - 50.0) + (x - 50.0) * (x - 50.0)) / 25.0);
                }
            }
            break;
      case 4:
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    double x = i * dx;
                    double y = j * dy;
                    h[j*Nx+i] = 10.0 + exp(-((y - 25.0) * (y - 25.0) + (x - 25.0) * (x - 25.0)) / 25.0) + exp(-((y - 75.0) * (y - 75.0) + (x - 75.0) * (x - 75.0)) / 25.0);
                }
            }
            break;
    }
}

void ShallowWater::WriteOutput(double* u, double* v, double* h, double dx, double dy) {
    
    // create file
    ofstream outfile("output.txt");
    outfile << fixed << setprecision(6);
    
    // write to file
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            double x = i*dx;
            double y = j*dy;
            int k = j*Nx+i; //col.major index
            
            outfile << x << " " << y << " " << u[k] << " " << v[k] << " " << h[k] << endl;
        }
        
        outfile << endl;
    }
    
    outfile.close();
    
}

// BLAS calculation of RHS (f)
void ShallowWater::BLASversionRHS(double* R1, double* R2, double* R3, double* u, double* v, double* h, double* dudx, double* dudy, double* dhdx, double* dvdx, double* dvdy, double* dhdy, double g) {

    double* ududx = new double[Nx*Ny];
    double* vdudy = new double[Nx*Ny];
    double* udvdx = new double[Nx*Ny];
    double* vdvdy = new double[Nx*Ny];
    double* udhdx = new double[Nx*Ny];
    double* vdhdy = new double[Nx*Ny];
    double* hdudx = new double[Nx*Ny];
    double* hdvdy = new double[Nx*Ny];
    
    double* gdhdx = new double[Nx*Ny];
    double* gdhdy = new double[Nx*Ny];
    
    double* sum1 = new double[Nx*Ny];
    double* sum2 = new double[Nx*Ny];
    double* sum3 = new double[Nx*Ny];
    double* sum4 = new double[Nx*Ny];
    double* sum5 = new double[Nx*Ny];
    double* sum6 = new double[Nx*Ny];
    double* sum7 = new double[Nx*Ny];
    
    // u * dudx
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, u, Nx, dudx, Nx, 0.0, ududx, Nx);
    // v * dudy
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, v, Nx, dudy, Nx, 0.0, vdudy, Nx);
    // u * dvdx
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, u, Nx, dvdx, Nx, 0.0, udvdx, Nx);
    // v * dvdy
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, v, Nx, dvdy, Nx, 0.0, vdvdy, Nx);
    // u * dhdx
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, u, Nx, dhdx, Nx, 0.0, udhdx, Nx);
    // v * dhdy
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, v, Nx, dhdy, Nx, 0.0, vdhdy, Nx);
    // h * dudx
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, h, Nx, dudx, Nx, 0.0, hdudx, Nx);
    // h * dvdy
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Nx, Ny, 1.0, h, Nx, dvdy, Nx, 0.0, hdvdy, Nx);
    // g * dhdx
    cblas_dscal(Nx, g, dhdx, 1);
    // g * dhdy
    cblas_dscal(Nx, g, dhdy, 1);
    
    // R1
    cblas_daxpy(Nx, 1.0, ududx, 1, vdudy, 1);
    cblas_daxpy(Nx, 1.0, sum1, 1, gdhdx, 1);
    cblas_dscal(Nx, -1.0, sum2, 1);
    
    // R2
    cblas_daxpy(Nx, 1.0, udvdx, 1, vdvdy, 1);
    cblas_daxpy(Nx, 1.0, sum3, 1, gdhdy, 1);
    cblas_dscal(Nx, -1.0, sum4, 1);
        
    // R3
    cblas_daxpy(Nx, 1.0, udhdx, 1, vdhdy, 1);
    cblas_daxpy(Nx, 1.0, hdudx, 1, hdvdy, 1);
    cblas_daxpy(Nx, 1.0, sum7, 1, sum6, 1);
    cblas_dscal(Nx, -1.0, sum7, 1);
    
    
}
