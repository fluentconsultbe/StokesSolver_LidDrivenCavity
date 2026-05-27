// StokesSolver.cpp : Defines the entry point for the application.
//

#include "StokesSolver.h"

#include <vector>
#include <cmath>
#include <iostream>

// A Simple Squared Field.
struct Field 
{
    int nx, ny;
    double dx, dy;
    std::vector<double> data;

    Field(
        int nx_, 
        int ny_, 
        double Lx, 
        double Ly
        )
        : 
        nx( nx_), 
        ny( ny_), 
        dx( Lx / (nx_ - 1)), 
        dy( Ly / (ny_ - 1)), 
        data( nx_* ny_, 0.0)
    {
    
    }

    double& operator()(int i, int j) { return data[j * nx + i]; }
    const double& operator()(int i, int j) const { return data[j * nx + i]; }
};

void applyBC(
    Field & u, 
    Field & v
    ) 
{
    int nx = u.nx, ny = u.ny;

    // No-slip on walls
    for (int i = 0; i < nx; i++) 
    {
        u(i, 0) = v(i, 0) = 0.0;
        u(i, ny - 1) = 1.0;  // lid velocity
        v(i, ny - 1) = 0.0;
    }

    for (int j = 0; j < ny; j++)
    {
        u(0, j) = v(0, j) = 0.0;
        u(nx - 1, j) = v(nx - 1, j) = 0.0;
    }
}

// Normalizes the pressure on the field.
void normalizePressure(
    Field & p
    ) 
{
    double sum = 0.0;

    // Compute the total pressure.
    for (double val : p.data) 
        sum += val;

    // Compute the mean pressure.
    double mean = sum / p.data.size();
    
    // Subtract mean from every cell.
    for ( double & val : p.data) 
        val -= mean;
}

int main() {

    // 
    int nx = 64, ny = 64;
    double Lx = 1.0, Ly = 1.0;

    Field u(nx, ny, Lx, Ly);
    Field v(nx, ny, Lx, Ly);
    Field p(nx, ny, Lx, Ly);

    applyBoundaryConditions(u, v);

    // TODO:
    // 1) Sæt en kraft f (fx driven cavity eller body force)
    // 2) Diskretisér Stokes-ligningerne
    // 3) Løs for u, v, p iterativt
    // 4) Enforce ∫ p dΩ = 0 (fx træk middelværdien fra)

    
    return 0;
}