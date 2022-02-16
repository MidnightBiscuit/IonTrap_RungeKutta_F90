program main
implicit none

! Trap dimensions
double precision, parameter :: r_0   = 2.865d-3/1.14511d0
double precision, parameter :: L     = 0.0140608827209d0
double precision, parameter :: d_0   = 4d-3 ! longueur du piege
! Trapping voltages
double precision, parameter :: Udc   = 1
double precision, parameter :: V_rf  = 4d1
double precision, parameter :: Omega = pi2 * 2.00d+06
! Integration constants:
integer         , parameter :: n_dt = 100 ! number of time steps per RF period
double precision, parameter :: dt_n = pi2 / n_dt
double precision, parameter :: dt   = pi2 / (n_dt*Omega) ! simu time step

! RK4 parameters
double precision, parameter :: gam = 1e-20 ! friction parameter
double precision, parameter :: h = dt
integer, parameter :: n_simu = n_dt*50 ! ! length of simulation
integer, parameter :: order = 2
integer :: n
double precision, allocatable, dimension(:) :: k1, k2, k3, k4
double precision :: dydx_temp=0, dydx2_temp=0

double precision, allocatable, dimension(:) :: t, r, v, a ! variables for ion dynamics

allocate(t(n_simu))
allocate(r(n_simu))
allocate(v(n_simu))
allocate(a(n_simu))
allocate(k1(order))
allocate(k2(order))
allocate(k3(order))
allocate(k4(order))
allocate(k_save(n_simu,4,2))

t(0) = 0
r(0) = r_0/4
v(0) = 0
a(0) = qe*2*V_rf /(mass*r_0**2)*cos(Omega*0)*r(0) - gam/mass*v(0)

do n = 0 , n_simu-1
	call update_rk4(derivs)
	t(n+1) = t(n) + h
	r(n+1) = r(n) + (k1(0)+2*k2(0)+2*k3(0)+k4(0))/6	! n is the time step
	v(n+1) = v(n) + (k1(1)+2*k2(1)+2*k3(1)+k4(1))/6	! n is the time step
	call derivs()
	a(n+1) = dydx2_temp
enddo
print*, 'dt', h

call save_data_to_file()

deallocate(t)
deallocate(r)
deallocate(v)
deallocate(a)
deallocate(k1)
deallocate(k2)
deallocate(k3)
deallocate(k4)
deallocate(k_save)

contains

! https://github.com/MidnightBiscuit/Runge_Kutta/blob/main/Runge_Kutta_ions.ipynb
! https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od

! Subroutine that provides the derivate for the physical problem
SUBROUTINE derivs()
implicit none
	dydx_temp  = v(n)
	dydx2_temp = qe*2*V_rf /(mass*r_0**2)*cos(Omega*t(n))*r(n) - gam/mass*v(n)
END

SUBROUTINE update_rk4(f) ! f is the function that derives
						 ! it provides derivate for each variable 
						 ! It should be changed so y = (pos, vel)
						 ! and then dydx = (dpos, dvel)
						 ! for each time step
						 ! for each particle
implicit none
EXTERNAL f
	call f(t(n),r(n),v(n))	
	k1(0) = h*dydx_temp
	k1(1) = h*dydx2_temp
	
	call f(t(n)+h/2,r(n)+0.5*k1(0),v(n)+0.5*k1(1))
	k2(0) = h*dydx_temp
	k2(1) = h*dydx2_temp
	
	call f(t(n)+h/2,r(n)+0.5*k2(0),v(n)+0.5*k2(1))
	k3(0) = h*dydx_temp
	k3(1) = h*dydx2_temp
	
	call f(t(n)+h,r(n)+k3(0),v(n)+k3(1))
	k4(0) = h*dydx_temp
	k4(1) = h*dydx2_temp
END

subroutine save_data_to_file()
implicit none
integer i

print*, 'n final',n
open(9, file = 'xva.dat', status='replace', access='sequential', action='write')
    do i = 0, n_simu-1
        write(9,221) t(i),&
                     r(i), &
                     v(i), &
                     a(i)
    enddo

221     format(4(1X,e24.17e3))
close(9)
end
	
end program main
