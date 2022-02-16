# Runge_Kutta repository

### Purpose

This project hosts programs demonstrating two integration algorithms used in molecular dynamics : velocity-Verlet and Runge-Kutta to forth-order (RK4). I developed those codes during my PhD thesis using the code developed by Jofre. The context is that of molecular dynamics simulation of trapped of ions in electric fields.

### Codes

Several codes are presented in this repository. Two codes are intended to demonstrate integration algorithms in the case of a single particle in 1D. A more advanced program is dedicated to the molecular dynamics of $N>1$ particles in a mix of radio-frequency and static fields, with Coulomb interaction and friction, using RK4.

##### 1D codes

```Verlet_Pendulum.ipynb```  is a short Python code demonstrating the velocity-Verlet algorithm for a simple pendulum. Created with jupyter notebook and Python 3.7.

```Runge_Kutta_ions.ipynb```  is a Python code demonstrating the velocity-Verlet algorithm and Runge-Kutta fourth order (RK4) algorithm for a single trapped ion in 1D. The motion equation is a second-order ODE ( $\ddot{x} = f(t,x,\dot{x})$ ). This notebook compares the dynamics obtained with velocity-Verlet and the dynamics obtained Runge-Kutta at fourth order. Runge-Kutta is implemented in two different ways. This code also read data output from the last Fortran code for another comparison. Created with jupyter notebook and Python 3.7.

```Runge_Kutta_1ion.F90```  is a Fortran code demonstrating the Runge-Kutta at fourth order algorithm with the same equation of motion of a single ion in 1D as in ```Runge_Kutta_ions.ipynb```. The data from the dynamics is generated with this code for later analysis.

Open ```Runge_Kutta_ions.ipynb``` to see comparison of all methods with graphs. Those programs can be compiled with ```gfortran```.

##### 3D Codes

```Runge_Kutta_Nions.F90``` is a functional program designed to render the molecular dynamics of a collection of $N>1$ charged particles in a radio-frequency and static field, with Coulomb interaction and friction, using RK4 integration algorithm.

This program requires ```ifort``` compiler and relies on ```OpenMP``` for parallelization.

### Licence

Licence Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).