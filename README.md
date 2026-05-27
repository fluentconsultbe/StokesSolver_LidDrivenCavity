# Stokes Solver for  Lid-Driven Cavity

Work in progress...

The aim of this project is to develop a C++ Stoke's solver for a Lid-Driven Cavity Flow.

The FEM simulation is developed in C++ using only the standard template library.

The code is provided in the _notebook_ subfolder of this repository, in the form of a CMake buildfile.

# 1. Problem Formulation

Draw Lid driven cavity bla bla.

1. # 1. 
![[Pasted image 20260527122321.png]]

# 2. Assumptions

Domæne:  

$$\Omega = [0,1] \times [0,1]$$
Randbetingelser:

- Top: (u = 1, ; v = 0) (den “drivende” låg)
- Øvrige vægge: (u = v = 0) (no‑slip)


### Tryk‑normalisering

As the pressure in incompressible Stokes flow is only defined up to an arbitrary constant, we enforce uniqueness by requiring the mean pressure to vanish:

$$\int_\Omega p , d\Omega = 0$$

# 3. Governing Equations

Stokes‑ligninger:

$$-\nu \nabla^2 \mathbf{u} + \nabla p = 0$$$$\nabla \cdot \mathbf{u} = 0$$

