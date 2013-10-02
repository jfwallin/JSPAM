program basic_run


  use init_module
  use parameters_module
  use setup_module
  use io_module
  use integrator


  implicit none

  character (len=20):: target_name
  real (kind=8) :: t0, time_interval
  integer (kind=4) :: nstep_local

!------------------------------------------------------
!
!
!

! set the disk parameters
  call RANDOM_SEED()


! set the target parameters
  call DEFAULT_PARAMETERS

  call CREATE_COLLISION

!
!     -----loop over the system for the output
!

! initialize rk routine
  call INIT_RKVAR(x0, mass1, mass2, epsilon1, epsilon2, theta1, phi1, &
       theta2, phi2, rscale1, rscale2, rout1, rout2, n, n1, n2)


  t0 = tstart

  nstep = int( (tend - t0) / h) + 2
  nstep_local = nstep 

  time_interval = (tend - t0) * 2

  nunit = 50
  call OCTAVE_PARAMETERS_OUT(mass1, theta1, phi1, rout1, mass2, &
       theta2, phi2, rout2, original_rv(1:3), original_rv(4:6), time_interval, x0(n,:), n, nunit)


! main integration loop
  iout = 0
  do istep = 1, nstep_local
    call TAKE_A_STEP
    if (mod(istep, 50) == 5) then
      print*,istep
!      call CREATE_IMAGES
    endif
  enddo

!      call CREATE_IMAGES
  write(fname,'(a5)') 'a.101'
  open(unit, file=fname)
  call OUTPUT_PARTICLES(unit, x0, mass1, mass2, &
       eps1, eps2, &
       n, n1, n2, &
       time, header_on)
  close(unit)

! this creates a simple script for animating the output with gnuplot
! gnuplot 
! i = 1
! j = 2
! load 'gscript
  if (.not. header_on) then
    call CREATE_GNUPLOT_SCRIPT(x0, iout)
  endif

! clean up memory
  deallocate(x0)
  deallocate(xout)
  call DEALLOCATE_RKVAR


! enddo


end program basic_run


