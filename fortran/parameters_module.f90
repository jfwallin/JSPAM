!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
module parameters_module
!     
!
!     -----Description: sets of the run parameters for the r3 body runs
!
!
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------


    real (kind=8), parameter :: mass_gm = 1.98892d44 ! gm
    real (kind=8), parameter :: mass_solar = 1.98892d33 ! gm
    real (kind=8), parameter :: distance = 4.6285203749999994d22 ! cm
    real (kind=8), parameter :: time_s = 2.733342473337471d15 ! s
    real (kind=8), parameter :: vel_unit = distance / time_s  ! cm/s
    real (kind=8), parameter :: pc = 3.08568025d18 !cm
    real (kind=8), parameter :: kpc = pc * 1000.0d0
    real (kind=8), parameter :: year = 365.25d0 * 24.0d0*3600.0d0 
    real (kind=8), parameter :: km = 1d5 ! cm
    real (kind=8), parameter :: vel_km_sec = vel_unit / km
    real (kind=8), parameter :: a_mss = distance / (time_s * time_s) / 100.0d0
    real (kind=8), parameter :: a0_mks = 1.2d-10
    real (kind=8), parameter :: a0 = a0_mks / a_mss
    real (kind=8), parameter :: pi = 3.141592653589793d0
    
    real (kind=8), parameter :: hbase = 0.001d0
    integer(kind=4) :: potential_type = 0

    integer, parameter :: ndim = 3
    real (kind=8) :: mass1 
    real (kind=8) :: epsilon1
    real (kind=8) :: rin1
    real (kind=8) :: rout1
    real (kind=8), dimension(3):: rscale1
    real (kind=8) :: theta1
    real (kind=8) :: phi1
    integer (kind=4) :: opt1
    real (kind=8) :: heat1
  
    real (kind=8) :: mass2 
    real (kind=8) :: epsilon2
    real (kind=8):: rin2
    real (kind=8) :: rout2
    real (kind=8), dimension(3):: rscale2
    real (kind=8) :: theta2
    real (kind=8) :: phi2
    integer (kind=4) :: opt2
    real (kind=8) :: heat2
    real (kind=8) :: tcurrent

    real (kind=8), dimension(:,:), allocatable :: x0, xout
    integer (kind=4) :: n, n1, n2

    real (kind=8) :: time,tstart,tend
    real (kind=8) :: inclination_degree
    real (kind=8) :: omega_degree
    real (kind=8) :: rmin
    real (kind=8) :: velocity_factor
    real (kind=8) :: mratio
    real (kind=8) :: secondary_size
    real (kind=8), dimension(6):: sec_vec
    logical :: use_sec_vec, tIsSet

    real (kind=8) :: h
    integer (kind=4) :: nstep
    integer (kind=4) :: nout
  
    integer (kind=4) :: iout
    integer (kind=4) :: unit
    integer (kind=4) :: istep

!
!     ----------------------------------------------------------------
!
CONTAINS
!
!     ----------------------------------------------------------------
!
 



!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  subroutine standard_galaxy1(mass1, epsilon1, rin1, rout1, rscale1, theta1, &
    phi1, opt1, heat1)
!     
!
!     -----Description: initializes the test parameters for the standard run
!
!
!          on input:  
!
!          on output:
!          
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------
!
  real (kind=8) :: mass1 
  real (kind=8) :: epsilon1
  real (kind=8) :: rin1
  real (kind=8) :: rout1
  real (kind=8), dimension(3):: rscale1
  real (kind=8) :: theta1
  real (kind=8) :: phi1
  integer (kind=4) :: opt1
  real (kind=8) :: heat1

!
!     ----------------------------------------------------------------
!

  ! disk profile #1
  mass1 = 1.0d0
  epsilon1  = 0.3d0
  rin1 = 0.05d0
  rout1 = 1.0d0
  rscale1 = 3.0d0
  theta1 = 0.0d0
  phi1 = 0.0d0
  opt1 = 1

  heat1 = 0.0d0
!  seed1 = 0.0d0

end subroutine standard_galaxy1


!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  subroutine standard_galaxy2(mass2, epsilon2, rin2, rout2, rscale2, theta2, &
    phi2, opt2, heat2)
!     
!
!     -----Description: 
!
!
!          on input:  
!
!          on output:
!          
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------
!
  real (kind=8) :: mass2 
  real (kind=8) :: epsilon2
  real (kind=8):: rin2
  real (kind=8) :: rout2
  real (kind=8), dimension(3):: rscale2
  real (kind=8) :: theta2
  real (kind=8) :: phi2
  integer (kind=4) :: opt2
  real (kind=8) :: heat2
  
!
!     ----------------------------------------------------------------
!
  ! disk profile #2
  mass2 = 1.0d0
  epsilon2  = 0.3d0
  rin2 = 0.05d0
  rout2 = 1.0d0
  rscale2 = 3.0d0
  theta2 = 0.0d0
  phi2 = 0.0d0
  opt2 = 1

  heat2 = 0.0d0
!  seed2 = 0.0d0

  return

  end subroutine standard_galaxy2


!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  subroutine test_collision(n, n1, n2, time, inclination_degree, omega_degree, &
    rmin, velocity_factor, h, nstep, nout) 
!     
!
!     -----Description: 
!
!
!          on input:  
!
!          on output:
!          
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------
!
  integer (kind=4) :: n, n1, n2
  real (kind=8) :: time  
  real (kind=8) :: inclination_degree
  real (kind=8) :: omega_degree
  real (kind=8) :: rmin
  real (kind=8) :: velocity_factor
  
  real (kind=8) :: h
  integer (kind=4) :: nstep
  integer (kind=4) :: nout

!  real (kind=8) :: lnl
!
!     ----------------------------------------------------------------
!

  ! collision parameters
  inclination_degree = 90.0d0
  omega_degree = 0.0d0
  rmin = 1.0d0
  velocity_factor = 1.0d0
  time = -3.0

  ! time step parameters
  h = hbase
  nout = 5
  nstep = 500

  ! particle numbers for the total and each disk
  n1 = 1000
  n2 = 1000
  n = n1 + n2


  return

end subroutine test_collision




!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  subroutine unrotate_frame(x_in, y_in, z_in, stheta, ctheta, sphi, cphi, xn, yn, zn)
!     
!
!     -----Description: rotates the frame using input angles
!          phi and theta
!
!
!          on input:  
!
!          on output:
!          
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------
!
    real (kind=8), intent(in) :: x_in, y_in, z_in
    real (kind=8), intent(in) :: stheta, ctheta, sphi, cphi
    real (kind=8), intent(out) :: xn, yn, zn

    real (kind=8) :: x1, y1, z1
    real (kind=8) :: x2, y2, z2
!
!     ----------------------------------------------------------------
!

    x1 = x_in
    y1 = y_in
    z1 = z_in

    x2  =   x1 * ctheta +  z1 * stheta
    y2  =   y1
    z2  =  -x1 * stheta +  z1 * ctheta

    xn  =   x2  * cphi -  y2 * sphi
    yn  =   x2  * sphi +  y2 * cphi
    zn  =   z2
    
    return

  end subroutine unrotate_frame




!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  subroutine rotate_frame(x_in, y_in, z_in, stheta, ctheta, sphi, cphi, xn, yn, zn)
!     
!
!     -----Description: rotates the frame using input angles
!          phi and theta
!
!
!          on input:  
!
!          on output:
!          
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------
!
    real (kind=8), intent(in) :: x_in, y_in, z_in
    real (kind=8), intent(in) :: stheta, ctheta, sphi, cphi
    real (kind=8), intent(out) :: xn, yn, zn

    real (kind=8) :: x1, y1, z1
    real (kind=8) :: x2, y2, z2
!
!     ----------------------------------------------------------------
!
    x1 = x_in
    y1 = y_in
    z1 = z_in
    
    x2  =   x1  * cphi+  y1 * sphi
    y2  =  -x1  * sphi+  y1 * cphi
    z2  =   z1
    
    xn  =   x2 * ctheta -  z2 * stheta
    yn  =   y2
    zn  =   x2 * stheta +  z2 * ctheta

    return

  end subroutine rotate_frame

















end module parameters_module


