! Math and Phys cte [SI units]:
double precision , parameter :: pi     = 3.141592653589793d0
double precision , parameter :: pi2    = 2.d0*pi
double precision , parameter :: inv_2pi= 1.0d0/pi2
double precision , parameter :: amu    = 1.66053886d-27   ![Kg]
double precision , parameter :: c      = 2.99792458d8     ! vitesse de la lumière
double precision , parameter :: hbar   = 1.055d-34
double precision , parameter :: hplanck= 6.62d-34
double precision , parameter :: qe     = 1.60217646d-19   ![C = sA] Elementary Charge
double precision , parameter :: kb     = 1.3806503d-23    ![m2 kg s-2] Boltzman cte
double precision , parameter :: ke     = 8.987551787d9    ![N m2 C-2] Coulomb cte =  1 / (4*pi*eps0)

! Particle parameters
double precision , parameter :: charge = qe*1
double precision , parameter :: mass   = amu*40
double precision , parameter :: q_m   = charge / mass
! Coulomb interaction cte (Not- Normalised):
double precision , parameter :: alpha = ke*charge*charge / mass

! OpenMP Parallel computing parameters

integer          , parameter :: ni   = 8 !number of threads
integer          , parameter :: n_ions    = 1024
integer          , parameter :: ia(ni)=[1,129,257,385,513,641,769,897]
integer          , parameter :: ib(ni)=[128,256,384,512,640,768,896,1024]

!~ integer          , parameter :: ni   = 4 !number of threads
!~ integer          , parameter :: n_ions    = 16
!~ integer          , parameter :: ia(ni)=[1,257,513,769]
!~ integer          , parameter :: ib(ni)=[256,512,768,1024]

!~ integer          , parameter :: n_ions    = 128
!~ integer          , parameter :: ia(ni)=[1,17,33,49,65,81,97,113]
!~ integer          , parameter :: ib(ni)=[16,32,48,64,80,96,112,128];

!~ integer          , parameter :: n_ions    = 64
!~ integer          , parameter :: ia(ni)=[1,9,17,25,33,41,49,57]
!~ integer          , parameter :: ib(ni)=[8,16,24,32,40,48,56,64];

!~ integer          , parameter :: ni   = 8 !number of threads
!~ integer          , parameter :: n_ions    = 16
!~ integer          , parameter :: ia(ni)=[1,3,5,7,9,11,13,15];
!~ integer          , parameter :: ib(ni)=[2,4,6,8,10,12,14,16];
integer :: im

! Simulations parameters
integer         , parameter :: n_dt = 100 ! number of time steps per RF period
double precision, parameter :: dt_n = pi2 / n_dt
double precision, parameter :: Omega = pi2 * 2.00d+06
double precision, parameter :: dt   = pi2 / (n_dt*Omega) ! simu time step
double precision, parameter :: dt1  = 0.5*dt
double precision, parameter :: dt2  = 0.5*dt*dt
integer         , parameter :: n_save_temp = 1*n_dt ! It has to be a multiple of n_dt
													! Number of points in one RF to save temp

integer         , parameter :: i_fly = 1e3 ! 2000
integer         , parameter :: i_simu = n_dt*2000 ! *50 length of simulation
integer         , parameter :: i_ss  = 1e1  ! 2000
double precision, parameter :: Temperature = 0.5d-3 ! Kelvin

!dec$ if defined(Eta0)
double precision, parameter :: eta = 1d-5   ![kg/s] friction  coefficient
!dec$ elseif defined(Eta1)
double precision, parameter :: eta = 1d-4   ![kg/s] friction  coefficient
!dec$ elseif defined(Eta2)
double precision, parameter :: eta = 1d-3   ![kg/s] friction  coefficient
!dec$ elseif defined(Eta3)
double precision, parameter :: eta = 1d-2   ![kg/s] friction  coefficient
!dec$ else
double precision, parameter :: eta = 1d-19   ![kg/s] friction  coefficient
!dec$ endif

! RK4 parameters
double precision, parameter :: gam = 1e-20 ! friction parameter
double precision, parameter :: h = dt

integer, parameter :: order = 2 ! order of ODE
double precision, dimension(4,order,n_ions,3) :: k_rk4
double precision :: k_sum(n_ions)
double precision, dimension(n_ions,3)  :: dydx_temp, dydx2_temp



character(len=150) :: str_file_to_load
character(len=100), parameter :: str_extra = '_nt1000'

! Trap dimensions
double precision, parameter :: r_0   = 2.865d-3/1.14511d0
double precision, parameter :: L     = 0.0140608827209d0
double precision, parameter :: d_0   = 4d-3 ! longueur du piege

! Trapping voltages
double precision, parameter :: Udc   = 1 ! axial static potential
double precision, parameter :: V_st  = 0.00d0 ! radial static potential
double precision, parameter :: V_rf  = 40 ! radial RF potential
double precision   , parameter :: wz_LC = pi2*100e3 !when 1V applied GiantMol  90806.9982303
double precision   , parameter :: wz2   = Udc * wz_LC**2
double precision, parameter :: a_cte_Quad_RF_LC(4) = (/ &
    -2*charge*V_st / (mass*r_0**2 ) + 0.5*wz2,    &
    +2*charge*V_st / (mass*r_0**2 ) + 0.5*wz2,    &
    +2*charge*V_rf / (mass*r_0**2 )             , &
     -wz2                                            /)



double precision, parameter :: softening = 1.00d-150

!~ ! Temperature cte:
double precision   , parameter :: m_kb_x_inv_n_ions  = mass/(kb*real(n_ions)*n_dt**2   )
double precision   , parameter :: m_kb_x_inv_n_ions2 = mass/(kb*real(n_ions)**2*n_dt**2)

!~ !dec$ if defined(Udc0)
!~ double precision   , parameter :: Udc   = 2.54    ! if other than 1V, then simply multiply by Udc
!~ !dec$ elseif defined(Udc1)
!~ double precision   , parameter :: Udc   = 5.71
!~ !dec$ elseif defined(Udc2)
!~ double precision   , parameter :: Udc   = 22.8
!~ !dec$ elseif defined(Udc3)
!~ double precision   , parameter :: Udc   = Udcstart*4
!~ !dec$ endif

!~ !dec$ if defined(Vrf0)
!~ double precision  , parameter :: V_rf  = 2.143d1
!~ character(len=100), parameter :: str_extra3 = ''
!~ !dec$ elseif defined(Vrf1)
!~ double precision  , parameter :: V_rf  = 3.215d1
!~ character(len=100), parameter :: str_extra3 = ''
!~ !dec$ elseif defined(Vrf2)
!~ double precision  , parameter :: V_rf  = 3.43d1
!~ character(len=100), parameter :: str_extra3 = ''
!~ !dec$ elseif defined(Vrf3)
!~ double precision  , parameter :: V_rf  = 6.3d1
!~ character(len=100), parameter :: str_extra3 = ''
!~ !dec$ endif

!~ double precision  , parameter :: V_rf  = 5.385d1 'q=0.5'
!~ double precision  , parameter :: V_rf  = 6.461d1 'q=0.6'
!~ double precision  , parameter :: V_rf  = 7.000d1 'q=0.65'


! Dynamics parameters
double precision :: t(i_simu)
double precision, dimension(i_simu,n_ions,3) :: r, v, a ! variables for ion dynamics
double precision, dimension(n_ions,3)        :: r_jj, v_jj, a_jj
double precision :: r_jj_temp(n_ions,3), v_jj_temp(n_ions,3)
double precision :: a_Coul_x(n_ions), a_Coul_y(n_ions), a_Coul_z(n_ions)
double precision :: a_Ust_x(n_ions), a_Ust_y(n_ions), a_Udc_z(n_ions)
double precision :: rji(3)



! S_aux_LC is source term from RF heating averaged over one RF period
double precision   , dimension(3)        :: S_aux_LC
double precision   , dimension(3)        :: T_CM__LC,T_aux_LC
double precision :: cosi(n_dt), sini(n_dt)
double precision  :: alphaL, betaL, time00, time_a, t_act
! Adrien 2020 07 08
double precision   , dimension(n_ions,3)   :: v_rf_avg, S_rf_avg
double precision   , dimension(3)          :: T_CM,T_aux
integer :: j_start  , size_Energy_array
integer :: j_save_temp = 0
integer :: i_save_temp = 1
integer :: j_end, iRF(3), j_save_E, idE_dt
integer :: j_save, j_save_next, j_aux

integer :: j,jj, nn, xyz_val ! loop variable
integer :: j_ions

!~ integer, dimension(ni)   :: ia , ib
character(len=150)       :: str_file_aux
double precision, allocatable, dimension(:,:) ::  save_temperature, save_trj, save_trj1, save_Energies

double precision        :: Ep, Ec, Ek, dE_dt, vdt, vdt1, vdt2, t_next_n_dt, S_rf_x, S_rf_y
logical :: flag_dt0
character(len=150)      :: str_file_Temp, str_file_xva, str_file_trj
!***********************************************************************
!~ double precision , parameter :: rSum_Adb = 188.57d-06 !(for N = 30ions) Inhomo


character(len=100), parameter :: str_extra_Lan   = '_5'
double precision  , parameter :: eta_mass = eta/mass
double precision  , parameter ::    w_p   = (dexp(-eta_mass*dt)+eta_mass*dt-1.0d0) /(eta_mass*dt*(1.0d0-dexp(-eta_mass*dt)))
double precision  , parameter ::    w_m   = 1-w_p
double precision  , parameter ::    cte1 = (kb*Temperature/mass) * (2*w_p*w_p * eta_mass * dt + w_p - w_m)
double precision  , parameter ::    cte2 = (kb*Temperature/mass) * (2*w_p*w_m * eta_mass * dt + w_m - w_p) ! b
double precision  , parameter ::    cte3 = (kb*Temperature/mass) * (2*w_m*w_m * eta_mass * dt + w_p - w_m) ! c
!dir$ if (integration_methode.eq.1) ! Langeving
double precision  , parameter ::    cte4 = w_p - 0.5d0 ! S vGB82
double precision  , parameter ::    dt_S  = dt*cte4
double precision  , parameter ::    dt_S2 = dt+dt_S
!dir$ elseif (integration_methode.eq.2) ! Langeving
double precision  , parameter ::    cte4 = 0           ! Impulse method
double precision  , parameter ::    dt_S = 0
double precision  , parameter ::    dt_S2 = dt
!dir$ endif
double precision  , parameter ::    e_eta_mass2 = dexp(-eta_mass*dt*0.5d0)
double precision  , parameter ::    e_eta_mass3 = (1.0d00-dexp(-eta_mass*dt))/(eta_mass*e_eta_mass2)
double precision  , parameter ::    cte13   = cte1 + cte3 !ac
!***********************************************************************

double precision, allocatable, dimension(:,:,:) ::  save_trjN
real    :: rand_num(1)
integer :: rand_seed
