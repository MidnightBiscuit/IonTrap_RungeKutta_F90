!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!~ Nom ......... : Runge_Kutta_Nions.F90
!~ Role ........ : Dynamique moléculaire de N ions
!~ 				   piégés dans un champ électrique
!~ 				   avec frottements (dx/dt).
!~ 				   Integration avec méthode Runge-Kutta ordre 4
!~ Auteurs ..... : Jofre Pedregosa-Gutierrez et Adrien Poindron
!~ Version ..... : 16 février 2022
!~
!~
!~ Compilation :
!~ ifort -xhost -O3 -align -save -qopenmp -mkl -fma Runge_Kutta_Nions.F90 -o a.out
!~ ./a.out #  > output.txt 2> error.txt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

include 'mkl_vsl.f90'
program main
USE MKL_VSL_TYPE
USE MKL_VSL
use omp_lib
implicit none

include 'variable_definition.f90'

!~ double precision, parameter :: T_threshold = 100* 3/m_kb_x_inv_n_ions ! Kelvin 


TYPE (VSL_STREAM_STATE) :: stream
integer brng,method
integer(kind=4) errcode

!~ allocate(k_save(i_simu,4,2))
! Initialise ions
print*, 'dt', h
print*, 'n_ions', n_ions

call create_files_names()
call save_info()
call initialize_cos_sin()
call init_random_generator()
t(1) = 0
v(1,:,:) = v_jj
a(1,:,:) = 0 ! qe*2*V_rf /(mass*r_0**2)*cos(Omega*0)*r(0) - gam/mass*v(0)
r(1,:,:) = r_jj

print*, ''
print*, '***********************************************************'
print*, '                 Initialisation of ions'
print*, '***********************************************************'
time00 = omp_get_wtime()
call initialize_ions()
print*, 'Initialisation time', omp_get_wtime() - time00
call save_ion_photo('r_endinit')
print*, 'jj final init', jj

print*, ''
print*, '***********************************************************'
print*, '                    Simulation of ions'
print*, '***********************************************************'
! proper simulation with computation
! of energies, temperature, size of cloud ...
call fly()
call save_ion_photo('r_endfly')

call save_data_to_file('xyz_1ion')
call save_temperature_to_file('T')

!~ deallocate(k_save)
print*, 'Total time', omp_get_wtime() - time00
print*, 'jj final simu', jj
contains

! https://github.com/MidnightBiscuit/Runge_Kutta/blob/main/Runge_Kutta_ions.ipynb
! https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od

SUBROUTINE update_rk4(f) ! f is the function that derives
						 ! it provides derivate for each variable 
						 ! It should be changed so y = (pos, vel)
						 ! and then dydx = (dpos, dvel)
						 ! for each time step
						 ! for each particle
implicit none
EXTERNAL f

	r_jj_temp = r_jj
	v_jj_temp = v_jj	
	call f(t(jj))		
	k_rk4(1,1,:,:) = h*dydx_temp
	k_rk4(1,2,:,:) = h*dydx2_temp

	r_jj_temp = r_jj + 0.5*k_rk4(1,1,:,:)
	v_jj_temp = v_jj + 0.5*k_rk4(1,2,:,:)
	call f(t(jj)+h/2)
	k_rk4(2,1,:,:) = h*dydx_temp
	k_rk4(2,2,:,:) = h*dydx2_temp
	
	r_jj_temp = r_jj + 0.5*k_rk4(2,1,:,:)
	v_jj_temp = v_jj + 0.5*k_rk4(2,2,:,:)
	call f(t(jj)+h/2)
	k_rk4(3,1,:,:) = h*dydx_temp
	k_rk4(3,2,:,:) = h*dydx2_temp
	
	r_jj_temp = r_jj + k_rk4(3,1,:,:)
	v_jj_temp = v_jj + k_rk4(3,2,:,:)
	call f(t(jj)+h)
	k_rk4(4,1,:,:) = h*dydx_temp
	k_rk4(4,2,:,:) = h*dydx2_temp
		
END

! Subroutine that provides the derivate for the physical problem
SUBROUTINE derivs(time_between_steps)
implicit none
double precision time_between_steps
double precision, dimension(2) :: cte_aux
double precision               :: r2inv,cos_jj ! , sin_jj,Ux, Uy,Uz

! Computation of Coulomb force

    !$omp parallel default(none) &
    !$omp private(im, r2inv,rji,nn,jj) &
    !$omp shared (a_Coul_x,a_Coul_y,a_Coul_z) & ! Ec_aux
    !$omp firstprivate(ia,ib,cte_aux,r_jj_temp)
        im = omp_get_thread_num()+1
do nn = ia(im), ib(im)
	do j_ions = 1, n_ions
		if (nn.ne.j_ions) then
			rji(1)  = r_jj_temp(j_ions,1) - r_jj_temp(nn,1)
			rji(2)  = r_jj_temp(j_ions,2) - r_jj_temp(nn,2)
			rji(3)  = r_jj_temp(j_ions,3) - r_jj_temp(nn,3)
			r2inv   = 1.d0/dsqrt(rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3))
			r2inv   = r2inv * r2inv * r2inv * alpha
			a_Coul_x(nn) = rji(1)*r2inv
			a_Coul_y(nn) = rji(2)*r2inv
			a_Coul_z(nn) = rji(3)*r2inv
		endif
	enddo
enddo
    !$omp end parallel

! Computation of force induced by static potentials
!~     cos_jj = cosi(mod(jj,n_dt)+1)
    cte_aux(1) = a_cte_Quad_RF_LC(1) + a_cte_Quad_RF_LC(3)*cos(Omega*time_between_steps)
    cte_aux(2) = a_cte_Quad_RF_LC(2) - a_cte_Quad_RF_LC(3)*cos(Omega*time_between_steps)
    ! DC on radial
    a_Ust_x =          cte_aux(1) * r_jj_temp(:,1)
    a_Ust_y =          cte_aux(2) * r_jj_temp(:,2)
    ! DC on axial
    a_Udc_z = a_cte_Quad_RF_LC(4) * r_jj_temp(:,3)

! Actual derivation using the set of two first-order equations
! according to the Runge-Kutta method
	dydx_temp(:,1)  = v_jj_temp(:,1)
	dydx_temp(:,2)  = v_jj_temp(:,2)
	dydx_temp(:,3)  = v_jj_temp(:,3)
	dydx2_temp(:,1) = cte_aux(1)*r_jj_temp(:,1) + a_Coul_x - gam/mass*v_jj_temp(:,1)
	dydx2_temp(:,2) = cte_aux(2)*r_jj_temp(:,2) + a_Coul_y - gam/mass*v_jj_temp(:,2)
	dydx2_temp(:,3) = a_Udc_z + a_Coul_z - gam/mass*v_jj_temp(:,3)
		
END

! Subroutine to compute the sum of slopes
! given by RK4
SUBROUTINE k_rk4_sum(r_val,d_val)
integer r_val, d_val
k_sum = (k_rk4(1,d_val,:,r_val)+2*k_rk4(2,d_val,:,r_val)+2*k_rk4(3,d_val,:,r_val)+k_rk4(4,d_val,:,r_val))/6
END

SUBROUTINE fly()

! For each time step compute the dynamics
! Using RK4 algorithm
! And equation of motion for trapped ions in linear Paul trap
! With RF and DC electric fields
! With friction force
do jj = i_init , i_simu-1    ! jj is the time step
	! call the RK4 subroutine
	! to update x,v and a
	! at next time step
	call update_rk4(derivs)
	call Measure_of_Temperature()
	call Check_save_temp()

	! r_jj, v_jj are the pos and velocity at current step
	! r,v,a are pos, vel and acc stored for the whole simu
	do xyz_val = 1,3
	call k_rk4_sum(xyz_val,1)
	r(jj+1,:,xyz_val) = r_jj(:,xyz_val) + k_sum
	call k_rk4_sum(xyz_val,2)
	v(jj+1,:,xyz_val) = v_jj(:,xyz_val) + k_sum
	enddo
	! then advance to next step
	t(jj+1) = t(jj) + h
	t_act = t_act + h
	! reset the jj variables to the next time step
	r_jj(:,:) = r(jj+1,:,:)
	v_jj(:,:) = v(jj+1,:,:)
	! compute acceleration at next step
	r_jj_temp = r_jj
	v_jj_temp = v_jj	
	call derivs(t(jj+1))
	a(jj+1,:,:) = dydx2_temp
enddo
END

subroutine Measure_of_Temperature()
implicit none
integer :: i
    v_rf_avg(:,1) = v_rf_avg(:,1) + v_jj(:,1);
    v_rf_avg(:,2) = v_rf_avg(:,2) + v_jj(:,2);
    v_rf_avg(:,3) = v_rf_avg(:,3) + v_jj(:,3);

    do i = 1,3
        iRF(i) = iRF(i) + 1;
        if (iRF(i) == n_dt) then
            iRF(i)        = 0;
            T_CM (i)      = sum(v_rf_avg(:, i) )**2
            T_aux(i)      = sum(v_rf_avg(:, i  )**2)
            v_rf_avg(:,i) = 0.0d0
        endif
    enddo
    
    return
end subroutine

subroutine Check_save_temp()
implicit none
    if (n_save_temp == i_save_temp) then
        i_save_temp = 0;
        j_save_temp = j_save_temp + 1

        save_temperature(j_save_temp, 1 )    = t_act
        save_temperature(j_save_temp, 2 :4 ) = T_CM
        save_temperature(j_save_temp, 5 :7 ) = T_aux
    else
        i_save_temp = i_save_temp + 1;
    endif
endsubroutine


subroutine init_random_generator()
implicit none
double precision, parameter :: l0(3) = [r_0/10,r_0/10,r_0/10] ! [10.0e-6,10.0e-6,10.0e-6]
double precision, parameter :: v0(3) = [1,1,1] ! [10.0e-6,10.0e-6,10.0e-6]
    brng=VSL_BRNG_MCG31
    method=VSL_RNG_METHOD_GAUSSIAN_ICDF

    CALL RANDOM_SEED()
    call random_number(rand_num)
    rand_seed = int(rand_num(1)*1e4)
    errcode=vslnewstream( stream, brng, rand_seed )
    
    errcode=vdrnggaussian( method, stream, n_ions,r_jj(:,1), 0.0d0, l0(1))
    errcode=vdrnggaussian( method, stream, n_ions,r_jj(:,2), 0.0d0, l0(2))
    errcode=vdrnggaussian( method, stream, n_ions,r_jj(:,3), 0.0d0, l0(3))
    
    errcode=vdrnggaussian( method, stream, n_ions,v_jj(:,1), 0.0d0, v0(1))
    errcode=vdrnggaussian( method, stream, n_ions,v_jj(:,2), 0.0d0, v0(2))
    errcode=vdrnggaussian( method, stream, n_ions,v_jj(:,3), 0.0d0, v0(3))

    a_jj = 0.0d0
    
    iRF(1) = 0; 
    iRF(2) = floor(n_dt/2.0); 
    iRF(3) = 0;
    
    return
end subroutine init_random_generator

SUBROUTINE initialize_ions()

! For each time step compute the dynamics
! Using RK4 algorithm
! And equation of motion for trapped ions in linear Paul trap
! With RF and DC electric fields
! With friction force
do jj = 1 , i_init-1    ! jj is the time step

	! call the RK4 subroutine
	! to update x,v and a
	! at next time step
	call update_rk4(derivs)

	! r_jj, v_jj are the pos and velocity at current step
	! r,v,a are pos, vel and acc stored for the whole simu
	do xyz_val = 1,3
	call k_rk4_sum(xyz_val,1)
	r(jj+1,:,xyz_val) = r_jj(:,xyz_val) + k_sum
	call k_rk4_sum(xyz_val,2)
	v(jj+1,:,xyz_val) = v_jj(:,xyz_val) + k_sum
	enddo
	! then advance to next step
	t(jj+1) = t(jj) + h
	print*, t_act
	t_act = t_act + h
	! reset the jj variables to the next time step
	r_jj(:,:) = r(jj+1,:,:)
	v_jj(:,:) = v(jj+1,:,:)
	! compute acceleration at next step
	r_jj_temp = r_jj
	v_jj_temp = v_jj	
	call derivs(t(jj+1))
	a(jj+1,:,:) = dydx2_temp
enddo
END

subroutine initialize_cos_sin() !Done
implicit none
integer :: i
do i = 1, n_dt
    cosi(i) = dcos((i-1)*dt_n)
    sini(i) = dsin((i-1)*dt_n)
enddo
end

subroutine save_ion_photo(str_file_name)
implicit none
character(len=*) :: str_file_name

open(1, file=trim(adjustl(str_file_name))//'_'//trim(adjustl(str_file_aux))//'.dat', status='replace', access='sequential', action='write')
	do nn = 1, n_ions
        write(1,220) r(jj,nn,1),r(jj,nn,2),r(jj,nn,3)
    enddo

220     format(3(1X,e25.17e3))
close(1)
end

subroutine save_data_to_file(str_file_name)
implicit none
integer i
character(len=*) :: str_file_name

open(9, file = trim(str_file_name)//'_'//trim(str_file_aux)//'.dat', status='replace', access='sequential', action='write')
    do i = 1, i_simu-1
        write(9,221) t(i),&
                     r(i,2,1),r(i,2,2),r(i,2,3), &
                     v(i,2,1),v(i,2,2),v(i,2,3), &
                     a(i,2,1),a(i,2,2),a(i,2,3)
    enddo

221     format(10(1X,e25.17e3))
close(9)
end

subroutine save_temperature_to_file(str_file_name)
implicit none
integer i
character(len=*) :: str_file_name

    open(3, file = trim(str_file_name)//'_'//trim(str_file_aux)//'.dat', status='replace', access='sequential', action='write')
    print*, 'j_save_temp Temp',j_save_temp
        do i = 1, j_save_temp
            write(3,222)  save_temperature(i,1),&
                          m_kb_x_inv_n_ions2*save_temperature(i,2:4), &
                          m_kb_x_inv_n_ions *save_temperature(i,5:7)
        enddo
222    format( 7(1X,e27.19e3))
    close(3)

end


subroutine save_info()
implicit none
character(len=130)  :: str_file_info_tmp

    str_file_info_tmp = 'parameters.info'

	open(unit = 10, file=trim(adjustl(str_file_info_tmp)), status='replace', access='sequential', action='write')
        write(10,'(i16  , 3x, "%Number of LC ions")')                   n_ions
        write(10,'(i16  , 3x, "%Last index in the integration loop")')  jj
        write(10,'(e16.9, 3x, "%Last time")')                           t(i_simu)
        write(10,'(i16  , 3x, "%Mass   of the LC ions")')               int(mass/amu)
        write(10,'(i16  , 3x, "%Charge of the LC ions")')               int(charge/qe)

        write(10,'(e16.9  , 3x, "%V_st[V]")')       V_st
		write(10,'(e16.9  , 3x, "%U_dc[V]")')       Udc
        write(10,'(e16.9  , 3x, "%V_rf[V]")')       V_rf
        write(10,'(e16.9  , 3x, "%Omega/2pi[Hz]")') Omega/pi2

		write(10,'(e16.9  , 3x, "%Gamma")')         Gam
		write(10,'(i16  , 3x, "%i_init")')          i_init
        write(10,'(i16  , 3x, "%i_simu")')          i_simu
    close(10)
end subroutine save_info

subroutine create_files_names()
implicit none
character(len=10 )  :: str_N, str_N_aux
character(len=130)  :: str_header
character(len=50 )  :: str_F

!~ 	str_trap_aux = '_RK4BlumelStyle'
	str_header = ''

    write(str_N,"(I4.4)") n_ions
    str_N_aux = 'N'//trim(adjustl(str_N))

    ! Friction coefficient
    write(str_F,"(d9.2)") abs(gam)
    str_F = '_F'//trim(adjustl(str_F))

    ! Names used by the "from scrath" subroutine
    str_file_aux = trim(adjustl(str_header)) & ! trim(adjustl(str_trap_aux)) &
                // trim(adjustl(str_N_aux))    &
                // trim(adjustl(str_F))    &
                // trim(adjustl(str_extra))

!~ print*, str_file_aux
!~ stop
end subroutine

	
end program main
