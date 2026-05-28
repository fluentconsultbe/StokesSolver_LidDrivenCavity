# Stokes Solver for  Lid-Driven Cavity

NOTICE!!! Work in progress...

The aim of this project is to develop a C++ Stoke's solver for a Lid-Driven Cavity Flow.

The FEM simulation is developed in C++ using only the standard template library.

The code is provided in the _StokesSolver_ subfolder in the form of a CMake buildfile.

# 1. Problem Formulation

The lid‑driven cavity is a classical benchmark problem in incompressible fluid dynamics.
It consists of a square domain filled with a viscous, incompressible fluid, where the top boundary (the “lid”) moves with a constant horizontal velocity while all other walls remain stationary. 

The motion of the lid induces a recirculating flow inside the cavity, generating characteristic vortical structures.

<img width="400" alt="Lid-Driven cavity" src="Figures/problem formulation.png" />

# 2. Assumptions

Domain is assumed to be a unit-rectangle:  

$$\Omega = [0,1] \times [0,1]$$
Boundary conditions:

- Top:  $u = 1, v = 0$ ( The Driven Lid )
- No-slip on walls: $u = v = 0$

### Pressure‑normalisation

As the pressure in incompressible Stokes flow is only defined up to an arbitrary constant, we enforce uniqueness by requiring the mean pressure to vanish:

$$\int_\Omega p , d\Omega = 0$$

# 3. Governing Equations

### Stokes‑equations:

$$-\nu \nabla^2 \mathbf{u} + \nabla p = 0$$$$\nabla \cdot \mathbf{u} = 0$$
### Pressure Poisson Equation


To enforce incompressibility we solve a pressure Poisson equation derived from the projection method. 

$$\nabla^2 p = \frac{1}{\Delta t},\nabla \cdot u^*$$
 in which $u^*$ is the __intermediate velocity field__ obtained after the viscous diffusion step, used to compute the pressure correction that enforces incompressibility.

This equation ensures that the corrected velocity field becomes divergence‑free.

I discretize the Laplacian using a standard 5‑point finite‑difference stencil on a uniform grid, and solve the resulting linear system iteratively using Gauss–Seidel relaxation. Once the pressure field is updated, the velocity is projected onto a solenoidal field:

$$u = u^* - \Delta t,\frac{\partial p}{\partial x}, \qquad v = v^* - \Delta t,\frac{\partial p}{\partial y}$$

This projection step enforces the incompressibility condition $\nabla \cdot u = 0$ and yields a physically consistent Stokes flow solution.

# 4. Discretization

We will construct a 2D‑grid with similar element distances, $h$, in both directions.

The Laplacianen in  2D is:

$$\nabla^2 u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}$$

We will use a 5-points Laplace-operator with central difference given by: 

$$\nabla^2 u_{i,j} \approx \frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}}{h^2}$$

This operator uses the center point of the current grid elements plus its four neighbours, aka left, right, up and down.

The formula for the field, $v_{i,j}$ is similar.
