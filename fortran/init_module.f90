!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
module init_module


use parameters_module
use setup_module
use io_module
use integrator
use df_module
!     
!
!     -----Description: sets up the rk
!
!
!
!     ----------------------------------------------------------------
!
    implicit none
!
!     -----Variable declarations
!          ---------------------


  character (len=10) :: fname
  integer (kind=4) :: iostat
  logical :: header_on
  real (kind=8) :: xmin, xmax, ymin, ymax

  real (kind=8), dimension(3) :: r_diff
  real (kind=8) :: angle1, angle2, angle3, angle4
  real (kind=8) :: separation
  logical :: solution_found
  real (kind=8) :: separation_target

  integer (kind=4), dimension(:,:), allocatable :: target_array
  real (kind=8), dimension(:,:), allocatable :: normalized_target
  integer (kind=4) :: xsize, ysize, maxsize
  integer (kind=4) :: nunit

  real (kind=8):: target_scale, target_separation, target_pa
  real (kind=8), dimension(ndim, ndim) :: matrix_a
  real (kind=8), dimension(ndim, ndim) :: matrix_b
  real (kind=8), dimension(:,:), allocatable :: projected
  real (kind=8) :: x_offset, y_offset
  real (kind=8) :: primary_a, primary_b, primary_angle
  real (kind=8) :: secondary_a, secondary_b, secondary_angle

  integer (kind=4) :: xsize_original, ysize_original
  real (kind=8) :: xscale, yscale, final_scale
  real (kind=8) :: xoff, yoff
  integer (kind=4), dimension(:,:), allocatable :: simulated_image_array
  real (kind=8), dimension(:,:), allocatable :: image_array
  real (kind=8) :: array_flux


  real (kind=8), dimension(:,:), allocatable :: residual_array
  integer (kind=4), dimension(:,:), allocatable :: residual_image_array
  logical :: do_output_now
  integer (kind=4) :: nangle4
  integer (kind=4):: iangle4
  real (kind=8) :: dangle4


  integer (kind=4) :: irun

  real (kind=8), dimension(3) :: pos, vel

  character (len=100) :: target_dat, target_pair, pnm_sim, pnm_target, pnm_resim, states_file

  real (kind=8) :: t1, p1, t2, p2, rx, ry, vz, rx_scale
  real (kind=8) :: rz, v_x, v_y
  real (kind=8) :: tzcross
  real (kind=8), dimension(6) :: minloc, zcrossloc

  real (kind=8) :: xs_reprocess, ys_reprocess, xc_reprocess, yc_reprocess
  real (kind=8), dimension(6) :: original_rv


!
!     ----------------------------------------------------------------
!
CONTAINS
!
!     ----------------------------------------------------------------
!


!---------------------------------------------------------------------

subroutine default_parameters

  potential_type=0

!
!     -----load the generic parameters for the system
!
  call STANDARD_GALAXY1(mass1, epsilon1, rin1, rout1, rscale1, theta1, phi1, opt1, heat1 )
  call STANDARD_GALAXY2(mass2, epsilon2, rin2, rout2, rscale2, theta2, phi2, opt2, heat2)
  call TEST_COLLISION(n, n1, n2, time, inclination_degree, omega_degree, &
    rmin, velocity_factor, h, nstep, nout) 

  call CUSTOM_COLLISION

  ! allocate the space for the particles - n+1 here ONLY
  allocate(x0(n+1,6), stat=iostat)
  allocate(xout(n+1,6), stat=iostat)
  allocate(projected(n+1,3), stat=iostat)

  return

end subroutine default_parameters
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine print_run

  implicit none

  call PRINT_PROFILE(1, rin1, rout1, rscale1, mass1, eps1, &
     theta1, phi1, opt1, heat1, n1)
  call PRINT_PROFILE(2, rin2, rout2, rscale2, mass2, eps2, &
     theta2, phi2, opt2, heat2, n2)
  call PRINT_COLLISION(n, time, inclination_degree, omega_degree, &
    rmin, velocity_factor, h, nstep, nout) 

  return

end subroutine print_run

!---------------------------------------------------------------------

subroutine create_collision

  implicit none
  real (kind=8), dimension(7) :: rv4min
  real (kind=8), dimension(4) :: tminVals
  real (kind=8) :: tmpT

  call INIT_DISTRIBUTION

  ! create the disks
!  call SET_DIFFQ2_PARAMETERS(phi1, theta1, phi2, theta2, rscale1, rscale2, rout1, rout2)

  call PROFILE(rin1, rout1, rscale1, 1, n1, mass1, eps1, &
       theta1, phi1, opt1, heat1, x0)

  call PROFILE(rin2, rout2, rscale2, n1+1, n, mass2, eps2, &
       theta2, phi2, opt2, heat2, x0)

  ! determine if we need to calculate tStart
  if( .NOT. tIsSet ) then
	rv4min = (/sec_vec(1), sec_vec(2), sec_vec(3), -sec_vec(4), -sec_vec(5), -sec_vec(6), 0.0d0/)
	tminVals = getTStart(rv4min, mass1, mass2, sqrt(eps1), sqrt(eps2) , h,-30.0d0 ,10.0d0*(rout1), rout1, rout2)
		
	tmpT = tminVals(1)
	if ( tmpT < -12.0 ) then
	  tmpT = -5
	endif
		
	if ( abs(tmpT) < h) then
	  tmpT = -5
	endif
		tstart = tmpT
                time = tstart
     tIsSet = .true.
  endif

  ! set the perturber galaxy position
  if( use_sec_vec ) then
    call PERTURBER_POSITION_VEC(sec_vec, mass1, mass2, eps1, eps2, h, n, n1, time, x0, original_rv)
  else
    call PERTURBER_POSITION(inclination_degree, omega_degree, rmin, &
         velocity_factor, &
         mass1, mass2, eps1, eps2, h, n, n1, time, x0, original_rv)
  endif





  return

end subroutine create_collision

!---------------------------------------------------------------------

subroutine create_images

  iout = iout + 1
  write(fname,'(a2,i3.3)') 'a.',iout
  open(unit, file=fname)
  call OUTPUT_PARTICLES(unit, x0, mass1, mass2, &
       eps1, eps2, &
       n, n1, n2, &
       time, header_on)
  close(unit)

  return

end subroutine create_images

!---------------------------------------------------------------------

subroutine take_a_step
  ! take a step
  

  h = hbase
  call WRAP_RK4(x0, h, xout)
  x0 = xout
  time = time + h

  return

end subroutine take_a_step
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  subroutine custom_collision()

    implicit none
    integer::narg
    character *300 buffer
    character *50 shortbuf

    tIsSet = .false.
! If command line arguments were passed, set them here
  narg = IARGC()
  if(narg .gt. 0) then
    print*,'custom collision ---------'
    call GETARG(1,buffer)
    shortbuf = TRIM(buffer)

    if(shortbuf .eq. '-f') then
! grab the filename
        call GETARG(2,buffer)
        shortbuf = TRIM(buffer)
        open(unit, file=shortbuf)
        CALL read_parameter_file(unit)
        close(unit)
    else
! this is a string specifying the state
        call parse_state_info_string(buffer)
        potential_type=1
        h = hbase
        tstart = -5
        tend = 0
        time = -5
        if(narg .gt. 1) then
          call GETARG(2,buffer)
          shortbuf = TRIM(buffer)
          read(shortbuf,*)tstart
          time = tstart
          tIsSet = .true.
        end if
    end if
  else
  
    phi1   = 5.0d0
    theta1 = 5.0d0
    rscale1 = 1.0d0
    rout1   = 1.0d0
    mass1   = 1.0d0
    epsilon1 = 0.3d0
    eps1 = epsilon1*epsilon1
    n1 = 1000
    heat1 = 0.0
    opt1 = 1

    phi2   = 0.0d0
    theta2 = 0.0d0
    rscale2 = 0.30d0
    rout2   = .50d0
    mass2   = .50d0
    epsilon2 = 0.3d0
    eps2 = epsilon2*epsilon2
    n2 = 500
    heat2 = 0.0
    opt2 = 1

    inclination_degree =20.0d0
    omega_degree = 0.0d0
    rmin = 0.90d0
    velocity_factor = 0.90d0
  
    h = hbase
    time = -5
    tstart = time
    tIsSet = .true.
  end if

  n = n1 + n2
  eps1 = epsilon1*epsilon1
  eps2 = epsilon2*epsilon2

  return
end subroutine custom_collision


!---------------------------------------------------------------------
!
   subroutine rotation_vector(theta, phi, xr, yr, zr)

     real (kind=8), intent(in) :: theta, phi
     real (kind=8), intent(out) :: xr, yr, zr

     real (kind=8) :: stheta, ctheta, sphi, cphi
     real (kind=8) :: x1, y1, z1

     stheta = sin(theta * pi / 180.0d0)
     ctheta = cos(theta * pi / 180.0d0)
     sphi = sin(phi * pi / 180.0d0)
     cphi = cos(phi * pi / 180.0d0)

     x1 = 0.0d0
     y1 = 0.0d0
     z1 = 1.0d0

     call ROTATE_FRAME(x1, y1, z1, stheta, ctheta, &
          sphi, cphi, xr, yr, zr)        
     return

   end subroutine rotation_vector

!---------------------------------------------------------------------
!
   subroutine cross_product(x1,y1,z1, x2,y2,z2, xc,yc,zc, mc)

     real (kind=8), intent(in) :: x1, y1, z1
     real (kind=8), intent(in) :: x2, y2, z2
     real (kind=8), intent(out) :: xc, yc, zc, mc
     
     xc =  (y1*z2 - z1*y2)
     yc = -(x1*z2 - z1*x2)
     zc =  (x1*y2 - y1*x2)
     
     mc = sqrt(xc*xc + yc*yc + zc*zc)
     return
   end subroutine cross_product

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
   subroutine rotate_position(x1, y1, z1, theta, phi, xr, yr, zr)

     real (kind=8), intent(in) :: x1, y1, z1
     real (kind=8), intent(in) :: theta, phi
     real (kind=8), intent(out) :: xr, yr, zr
     real (kind=8) :: stheta, ctheta, sphi, cphi


     stheta = sin(theta * pi / 180.0d0)
     ctheta = cos(theta * pi / 180.0d0)
     sphi = sin(phi * pi / 180.0d0)
     cphi = cos(phi * pi / 180.0d0)

     call ROTATE_FRAME(x1, y1, z1, stheta, ctheta, &
          sphi, cphi, xr, yr, zr)        

     return

   end subroutine rotate_position

  
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!


!---------------------------------------------------------------------
!---------------------------------------------------------------------



 

end module init_module

