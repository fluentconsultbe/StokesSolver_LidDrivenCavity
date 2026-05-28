# Stokes Solver for  Lid-Driven Cavity

The aim of this project is to develop a C++ Stoke's solver for a Lid-Driven Cavity Flow.

The FEM simulation is developed in C++ using only the standard template library.

The code is provided in the _StokesSolver_ subfolder in the form of a CMake buildfile.

# 1. Problem Formulation

The lid‑driven cavity is a classical benchmark problem in incompressible fluid dynamics.
It consists of a square domain filled with a viscous, incompressible fluid, where the top boundary (the “lid”) moves with a constant horizontal velocity while all other walls remain stationary. 

The motion of the lid induces a recirculating flow inside the cavity, generating characteristic vortical structures.

<img width="400" alt="Lid-Driven cavity" src="Figures/problem formulation.png" />

# 2. Assumptions

The domain is assumed to be a unit-square:  

$$\Omega = [0,1] \times [0,1]$$

We are imposing the following boundary conditions:

- Top: $u = 1, v = 0$ ( The Driven Lid )
- Walls: $u = v = 0$ ( No-slip )

### Pressure‑normalisation

As the pressure in incompressible Stokes flow is only defined up to an arbitrary constant, we enforce uniqueness by requiring the mean pressure to vanish:

$$\int_\Omega p \hspace{3pt} d\Omega = 0$$

# 3. Governing Equations

The lid‑driven cavity flow is governed by the **steady incompressible Stokes equations**, which describe slow, viscosity‑dominated motion of a Newtonian fluid. The system consists of a momentum balance and the incompressibility constraint:

### Momentum equation
$$-\nu \nabla^2 \mathbf{u} + \nabla \mathbf{p} = \mathbf{0}$$

### Continuity equation ( incompressibility )
$$\nabla \cdot \mathbf{u} = 0$$

Here:

- $\mathbf{u}=(u,v)$ is the velocity field
    
- $\mathbf{p}$ is the pressure
    
- ν is the kinematic viscosity
---
### Pressure Poisson Equation

To enforce incompressibility numerically, the solver uses a **projection method**.  
First, an intermediate velocity field $\mathbf{u}^*$ is computed by applying only the viscous diffusion step:

$$\mathbf{u}^* = \mathbf{u}^n + \Delta t \hspace{2pt} \nu \nabla^2 \mathbf{u}^n$$

This intermediate field generally **does not** satisfy the divergence‑free condition.  
To correct this, we solve the **pressure Poisson equation**:

$$\nabla^2 p = \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^*$$

Solving this equation yields a pressure field that enforces incompressibility.

---
## Velocity Projection

Once the pressure field is obtained, the velocity is projected onto a **solenoidal** (divergence‑free) field:

$$u = u^* - \Delta t \frac{\partial p}{\partial x}, \qquad v = v^* - \Delta t \frac{\partial p}{\partial y}$$

This projection step ensures:

$$\nabla \cdot \mathbf{u} = 0$$

and produces a physically consistent Stokes flow solution.

---

### Numerical Solution of the Poisson Equation

The pressure Poisson equation forms a large sparse linear system.  
In this solver, it is iteratively approximated using **Gauss–Seidel relaxation**, which updates each pressure value based on its neighbours:

$$p_{i,j}^{(k+1)} = \frac{1}{4} \left( p_{i+1,j}^{(k)} + p_{i-1,j}^{(k)} + p_{i,j+1}^{(k)} + p_{i,j-1}^{(k)}

- h^2\hspace{2pt} \text{rhs}_{i,j} \right)$$

where:

$$\text{rhs}_{i,j} = \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^*_{i,j}$$

Repeated sweeps drive the solution toward the correct pressure field.

# 4. Discretization

We construct a uniform 2D grid with equal spacing, $h$, in both the x- and y-directions.

The Laplacian in  2D is:

$$\nabla^2 u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}$$

We approximate this operator using the **five‑point central‑difference stencil**:

$$\nabla^2 u_{i,j} \approx \frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}}{h^2}$$

This stencil uses the value at the current grid point $(i,j)$ and its four direct neighbours (left, right, top, bottom). 

The same discretization is applied to the second component of the velocity field, $v_{i,j}$.
