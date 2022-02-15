# Runge_Kutta repository

### Purpose

This project is intended to host programs demonstrating two integration algorithms used in molecular dynamics : velocity-Verlet and Runge-Kutta 4. I developped those codes during my PhD thesis. They are intended to be part of a code designed for molecular dynamics simulation (trapping of ions in electric fields).

The present code is  a simple version rendering a 1D single particle dynamics in a radio-frequency field. It is the support for further implementation of poly-ion dynamics, which requires to consider Coulomb repulsion. The transition to 3D is also a requirement in order to provide an actual functional code.



### Codes

```Verlet_Pendulum.ipynb```  is a short Python code demonstrating the velocity-Verlet algorithm for a simple pendulum. Created with jupyter notebook and Python 3.7.

```Runge_Kutta_ions.ipynb```  is a Python code demonstrating the velocity-Verlet algorithm and Runge-Kutta fourth order (RK4) algorithm for a trapped ion whose motion equation is a second-order ODE ($\ddot{x} = f(t,x,\dot{x})$). It compares the dynamics obtained with velocity-Verlet and the dynamics obtained Runge-Kutta at fourth order. Runge-Kutta is implemented in two differejnt ways. This code also read data output from the last code for another comparison. Created with jupyter notebook and Python 3.7.

```Runge_Kutta_ions.F90```  is a Fortran code demonstrating the Runge-Kutta at fourth order algorithm the same equation of motion of a single ion in 1D as in ```Runge_Kutta_ions.ipynb```. The data from the dynamics is generated with this code for later analysis.



Open ```Runge_Kutta_ions.ipynb``` to see comparison of all methods with graphs.



### Licence

Licence Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).