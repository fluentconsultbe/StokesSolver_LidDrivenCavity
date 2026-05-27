# Stokes Solver for  Lid-Driven Cavity

Work in progress...

The aim of this project is to develop a C++ Stoke's solver for a the laminar velocity profile in a circular pipe at low Reynold's number using the __finite element method__.

The FEM simulation is developed in C++ using only the standard template library.

https://scikit-fem.readthedocs.io/en/latest/index.html

The code is provided in the _notebook_ subfolder of this repository, in the jupyter notebook _flow.ipynb_. It is assumed that the user has a standard Anaconda installation.
A stokes solver in C++

# 1. Problem Formulation

Draw Lid driven cavity bla bla.

1. # 1. 


# 2. Assumptions

Domæne:  

$$\Omega = [0,1] \times [0,1]$$
Randbetingelser:

- Top: (u = 1, ; v = 0) (den “drivende” låg)
- Øvrige vægge: (u = v = 0) (no‑slip)


# 3. Governing Equations

Stokes‑ligninger:

$$-\nu \nabla^2 \mathbf{u} + \nabla p = 0$$$$\nabla \cdot \mathbf{u} = 0$$

Tryk‑normalisering:

$$\int_\Omega p \hspace{3pt}, d\Omega = 0$$
