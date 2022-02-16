program main
implicit none

integer          , parameter :: n_ions    = 1024
real    :: rand_num(1)
integer :: rand_seed
brng=VSL_BRNG_MCG31
method=VSL_RNG_METHOD_GAUSSIAN_ICDF

CALL RANDOM_SEED()
call random_number(rand_num)
rand_seed = int(rand_num(1)*1e4)
errcode=vslnewstream( stream, brng, rand_seed )
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
integer, parameter :: i_simu = n_dt*50 ! length of simulation
integer, parameter :: order = 2 ! order of ODE
integer :: jj, nn ! loop variable
double precision, allocatable, dimension(:,:,:) :: k1, k2, k3, k4
double precision :: dydx_temp=0, dydx2_temp=0

double precision, allocatable, dimension(:) :: t
double precision, allocatable, dimension(:,:,:) :: r, v, a ! variables for ion dynamics
double precision   , dimension(n_ions)   ::  r_z, r_x, r_r, r_y, &
                                             v_z, v_x, v_r, v_y,   &
                                             a1_z, a1_x, a1_r,a1_y, &
                                             a2_z, a2_x, a2_r,a2_y, &

allocate(t(i_simu))
allocate(r(i_simu,n_ions,3))
allocate(v(i_simu,n_ions,3))
allocate(a(i_simu,n_ions,3))
allocate(k1(order,3))
allocate(k2(order,3))
allocate(k3(order,3))
allocate(k4(order,3))
allocate(k_save(i_simu,4,2))

t(0) = 0
v(0,:,:) = 0
a(0,:,:) = 0 ! qe*2*V_rf /(mass*r_0**2)*cos(Omega*0)*r(0) - gam/mass*v(0)
r_x = r(0,:,0)
r_y = r(0,:,1)
r_z = r(0,:,2)


do jj = 1 , i_simu-1    ! jj is the time step
	do nn = 0, n_ions-1 ! nn is the ion index
		call update_rk4(derivs)
		t(jj+1) = t(jj) + h
		
		r_x(nn) = r_x(nn) + (k1(0,0)+2*k2(0,0)+2*k3(0,0)+k4(0,0))/6
		r_y(nn) = r_y(nn) + (k1(0,1)+2*k2(0,1)+2*k3(0,1)+k4(0,1))/6
		r_z(nn) = r_z(nn) + (k1(0,2)+2*k2(0,2)+2*k3(0,2)+k4(0,2))/6
			
		v_x(nn) = v(nn) + (k1(1,0)+2*k2(1,0)+2*k3(1,0)+k4(1,0))/6
		v_y(nn) = v(nn) + (k1(1,1)+2*k2(1,1)+2*k3(1,1)+k4(1,1))/6
		v_z(nn) = v(nn) + (k1(1,2)+2*k2(1,2)+2*k3(1,2)+k4(1,2))/6
		
		r_temp = (/ r_x(nn),r_y(nn),r_z(nn) /)
		v_temp = (/ v_x(nn),v_y(nn),v_z(nn) /)
		call derivs(t(n),r_temp,v_temp))
		a_x(nn) = dydx2_x_temp
		a_y(nn) = dydx2_y_temp
		a_z(nn) = dydx2_z_temp
	enddo
	r(jj,:,0) = r_x(nn); r(jj,:,1) = r_y(nn); r(jj,:,2) = r_z(nn)
	v(jj,:,0) = v_x(nn); v(jj,:,1) = v_y(nn); v(jj,:,2) = v_z(nn)
	a(jj,:,0) = a_x(nn); a(jj,:,1) = a_y(nn); a(jj,:,2) = v_z(nn)
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
	dydx_x_temp  = v_temp(0)
	dydx_y_temp  = v_temp(1)
	dydx_z_temp  = v_temp(2)
	dydx2_x_temp = qe*2*V_rf /(mass*r_0**2)*cos(Omega*t(n))*r_temp(0) - gam/mass*v_temp(0)
	dydx2_y_temp = qe*2*V_rf /(mass*r_0**2)*cos(Omega*t(n))*r_temp(1) - gam/mass*v_temp(1)
	dydx2_z_temp = qe*2*V_rf /(mass*r_0**2)*cos(Omega*t(n))*r_temp(2) - gam/mass*v_temp(2)
END

SUBROUTINE update_rk4(f) ! f is the function that derives
						 ! it provides derivate for each variable 
						 ! It should be changed so y = (pos, vel)
						 ! and then dydx = (dpos, dvel)
						 ! for each time step
						 ! for each particle
implicit none
EXTERNAL f
	r_temp = (/ r_x(nn),r_y(nn),r_z(nn) /)
	v_temp = (/ v_x(nn),v_y(nn),v_z(nn) /)
	call f(t(n),r_temp,v_temp))		
	k1(0,0) = h*dydx_x_temp
	k1(1,0) = h*dydx2_x_temp
	k1(0,1) = h*dydx_y_temp
	k1(1,1) = h*dydx2_y_temp
	k1(0,2) = h*dydx_z_temp
	k1(1,2) = h*dydx2_z_temp
		
	r_temp = (/ r_x(nn)+0.5*k1(0,0),r_y(nn)+0.5*k1(0,1),r_z(nn)+0.5*k1(0,2) /)
	v_temp = (/ v_x(nn)+0.5*k1(1,0),v_y(nn)+0.5*k1(1,1),v_z(nn)+0.5*k1(1,2) /)
	call f(t(n)+h/2,r_temp,v_temp))
	k2(0,0) = h*dydx_x_temp
	k2(1,0) = h*dydx2_x_temp
	k2(0,1) = h*dydx_y_temp
	k2(1,1) = h*dydx2_y_temp
	k2(0,2) = h*dydx_z_temp
	k2(1,2) = h*dydx2_z_temp
	
	r_temp = (/ r_x(nn)+0.5*k2(0,0),r_y(nn)+0.5*k2(0,1),r_z(nn)+0.5*k2(0,2) /)
	v_temp = (/ v_x(nn)+0.5*k2(1,0),v_y(nn)+0.5*k2(1,1),v_z(nn)+0.5*k2(1,2) /)
	call f(t(n)+h/2,r_temp,v_temp))
	k3(0,0) = h*dydx_x_temp
	k3(1,0) = h*dydx2_x_temp
	k3(0,1) = h*dydx_y_temp
	k3(1,1) = h*dydx2_y_temp
	k3(0,2) = h*dydx_z_temp
	k3(1,2) = h*dydx2_z_temp
	
	r_temp = (/ r_x(nn)+k3(0,0),r_y(nn)+k3(0,1),r_z(nn)+k3(0,2) /)
	v_temp = (/ v_x(nn)+k3(1,0),v_y(nn)+k3(1,1),v_z(nn)+k3(1,2) /)
	call f(t(n)+h,r_temp,v_temp))
	k4(0,0) = h*dydx_x_temp
	k4(1,0) = h*dydx2_x_temp
	k4(0,1) = h*dydx_y_temp
	k4(1,1) = h*dydx2_y_temp
	k4(0,2) = h*dydx_z_temp
	k4(1,2) = h*dydx2_z_temp
END

subroutine init_random_generator()
implicit none
double precision, parameter :: l0(3) = [10.0e-6,10.0e-6,10.0e-6];
    brng=VSL_BRNG_MCG31
    method=VSL_RNG_METHOD_GAUSSIAN_ICDF

    CALL RANDOM_SEED()
    call random_number(rand_num)
    rand_seed = int(rand_num(1)*1e4)
    errcode=vslnewstream( stream, brng, rand_seed )
    
    errcode=vdrnggaussian( method, stream, n_ions,r(0,:,0), 0.0d0, l0(1))
    errcode=vdrnggaussian( method, stream, n_ions,r(0,:,1), 0.0d0, l0(3))
    errcode=vdrnggaussian( method, stream, n_ions,r(0,:,2), 0.0d0, l0(2))

    v_z  = 0.0d0; a1_z = 0.0d0; a2_z = 0.0d0;
    v_r  = 0.0d0; a1_r = 0.0d0; a2_r = 0.0d0;
    v_y  = 0.0d0; a1_y = 0.0d0; a2_y = 0.0d0;
    
    return
end subroutine init_random_generator

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
