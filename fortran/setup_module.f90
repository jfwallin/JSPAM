module setup_module

  use parameters_module
  use io_module
  use df_module

  real (kind=8) :: phi_i1, phi_i2, theta_i1, theta_i2
  real (kind=8), dimension(3) :: rscale_i1, rscale_i2
  real (kind=8) :: rrout1, rrout2


contains

!---------------------------------------------------
! Use this method so that neither the caller of rk41
! nor implementation need to know which potential
! is being used
!---------------------------------------------------
subroutine  wrap_rk41(xx0, h, mass1, mass2, eps1, eps2, xxe)
    implicit none
    real(kind=8), dimension(6), intent(in):: xx0
    real(kind=8), intent(in) :: h
    real(kind=8), intent(in) :: mass1, mass2
    real(kind=8), intent(in) :: eps1, eps2
    real(kind=8), dimension(6), intent(out):: xxe

    if ( potential_type .eq. 0) then
      call rk41(xx0, h, mass1, mass2, eps1, eps2, xxe, diffq_spm)
    else if ( potential_type .eq. 1) then
      call rk41(xx0, h, mass1, mass2, eps1, eps2, xxe, diffq_nbi)
    else if ( potential_type .eq. 2) then
      call rk41(xx0, h, mass1, mass2, eps1, eps2, xxe, diffq_mond)
    endif

    return

end subroutine wrap_rk41


!---------------------------------------------------
  subroutine perturber_position(inclination_degree, omega_degree, &
       rmin, velocity_factor, &
       mass1, mass2, eps1, eps2, h, n, n1, t0, x0, original_rv)

    implicit none
    real (kind=8), intent(in) :: inclination_degree
    real (kind=8), intent(in) :: omega_degree
    real (kind=8), intent(in) :: rmin
    real (kind=8), intent(in) :: velocity_factor
    real (kind=8), intent(in) :: mass1
    real (kind=8), intent(in) :: mass2
    real (kind=8), intent(in) :: eps1
    real (kind=8), intent(in) :: eps2
    real (kind=8), intent(in) :: h
    real (kind=8), intent(in) :: t0
    real (kind=8), dimension(6), intent(out) :: original_rv
    integer (kind=4), intent(inout) :: n
    integer (kind=4), intent(in) :: n1
    real (kind=8), dimension(:,:), intent(inout) :: x0

    real (kind=8) :: en, v1
    real (kind=8), dimension(6) :: xx0  !!!, xxe
!!!    integer(kind=4) :: i
    real (kind=8) :: omega, incl  !!!, dist1
!!!    real (kind=8) :: tcurrent
    real (kind=8) :: epsilon1, epsilon2

    epsilon1 = sqrt(eps1)
    epsilon2 = sqrt(eps2)

! change inclination and omega into radians
    incl = inclination_degree * pi/180.0d0
    omega = omega_degree * pi/180.0d0

! energy from mass1
    if ( eps1 .gt. 0.0d0) then
!      en = mass1/eps1 * (pi/2.0d0 - atan(rmin/eps1 ) ) 
      en = mass1/epsilon1 * (pi/2.0d0 - atan(rmin/epsilon1 ) ) 
    else
      en = mass1 / rmin
    endif

! energy from mass2
    if ( eps2 .gt. 0.) then
!      en = en + mass2/eps1 * (pi/2.0d0 - atan(rmin/eps2 ) )
      en = en + mass2/epsilon2 * (pi/2.0d0 - atan(rmin/epsilon2 ) )
    else
      en = en + mass2 / rmin
    endif

! calculate escape velocity and velocity at rmin
    v1 = sqrt(2.0d0 * en)

    v1 = sqrt(2.0d0)*circular_velocity(mass1+mass2,rmin,rrout1,eps1)

! adjust velocity for MOND
    v1= -v1 * velocity_factor


!	setup the transformaton based on the matrix in
!	fundimentals of astrodynamics p. 82 by
!	bates, mueller, and white (1971)
!
    xx0(1) = cos(omega) * rmin
    xx0(2) = sin(omega) * cos(incl) * rmin
    xx0(3) = sin(omega) * sin(incl) * rmin

    xx0(4) = -sin(omega) * v1
    xx0(5) =  cos(omega) * cos(incl) * v1
    xx0(6) =  cos(omega) * sin(incl) * v1

! update sec_vec
    sec_vec = xx0
    sec_vec(4) = -sec_vec(4)
    sec_vec(5) = -sec_vec(5)
    sec_vec(6) = -sec_vec(6)

    call perturber_position_vec(xx0, mass1,mass2,eps1,eps2, h, n , n1, t0, x0, original_rv)

end subroutine perturber_position

  subroutine perturber_position_vec(xx0, mass1, mass2, eps1, eps2, h, n, n1, t0, x0, original_rv)

    implicit none
    real (kind=8), intent(in) :: mass1
    real (kind=8), intent(in) :: mass2
    real (kind=8), intent(in) :: eps1
    real (kind=8), intent(in) :: eps2
    real (kind=8), intent(in) :: h
    real (kind=8), intent(in) :: t0
    integer (kind=4), intent(inout) :: n
    integer (kind=4), intent(in) :: n1
    real (kind=8), dimension(:,:), intent(inout) :: x0
    real (kind=8), dimension(6),intent(inout) :: xx0
    real (kind=8), dimension(6),intent(inout) :: original_rv


    real (kind=8), dimension(6) :: xxe
    integer(kind=4) :: i
    real (kind=8) :: dist1
    real (kind=8) :: tcurrent
    real (kind=8) :: epsilon1, epsilon2


!   copy the original input vector
    original_rv = xx0


    epsilon1 = sqrt(eps1)
    epsilon2 = sqrt(eps2)

!   reverse the velocity for backward integration 
    xx0(4) = -xx0(4)
    xx0(5) = -xx0(5)
    xx0(6) = -xx0(6)

!
!       now move position back to t0 from t=0.0
!
    tcurrent = 0
    do while (t0 < tcurrent) 
      call wrap_rk41(xx0, h, mass1, mass2, eps1, eps2, xxe)
      dist1 = sqrt(xx0(1)**2+ xx0(2)**2 + xx0(3)**2 )
      xx0 = xxe
      tcurrent  = tcurrent - h
    enddo

!
!   reverse the velocity for forward integration 
    xx0(4) = -xx0(4)
    xx0(5) = -xx0(5)
    xx0(6) = -xx0(6)

!	now adjust the test particles from the 
!	second disk to the proper velocity and positions
!
    if (n > n1) then
      do i = n1+1, n
        x0(i,:) = x0(i,:) + xx0(:)
      enddo
    endif

!
! 	include the perturbing galaxy
!
    n = n + 1
    x0(n,:) = xx0



    return
  end subroutine perturber_position_vec


!---------------------------------------------------
  subroutine reset_perturber_position(pos, vel, &
       mass1, mass2, eps1, eps2, h, n, n1, t0, x0, tend, &
       minloc, zcrossloc, tzcross)

    implicit none

    real (kind=8), dimension(3) :: pos, vel
    real (kind=8), intent(in) :: mass1
    real (kind=8), intent(in) :: mass2
    real (kind=8), intent(in) :: eps1
    real (kind=8), intent(in) :: eps2
    real (kind=8), intent(in) :: h
    real (kind=8), intent(in) :: t0
    integer (kind=4), intent(inout) :: n
    integer (kind=4), intent(in) :: n1
    real (kind=8), dimension(:,:), intent(inout) :: x0
    real (kind=8), intent(out) :: tend, tzcross
    real (kind=8), dimension(6), intent(out) :: minloc, zcrossloc

    real (kind=8), dimension(6) :: xx0, xxe
    integer(kind=4) :: i, istep

    real (kind=8) :: tcurrent, dist, dist_old, zdist, zdist_old
    logical (kind=4) :: zmin_flag, min_flag
    integer (kind=4) :: izmin


    real (kind=8) :: rtime
  
!       set the positions and velocity of the companion

    xx0(1:3) =   pos(1:3)
    xx0(4:6) =  -vel(1:3)

!
!       now move position back to t0 from t=0.0
!
    tcurrent = 0.0d0
    tend = 0.0d0
    istep = 0
    rtime = 0.0d0

    dist_old = 1.0d10
    zdist_old = 1.0d10
    min_flag = .true.
    zmin_flag = .true.


!    xx0(4:6) = -xx0(4:6)

    do while (t0 < tcurrent  ) 
      call wrap_rk41(xx0, h, mass1, mass2, eps1, eps2, xxe)

      write(18,*) xx0

      dist = sqrt(xx0(1)**2+ xx0(2)**2 + xx0(3)**2 )
      xx0 = xxe


!      if the distance is larger than the last step, update the
!      clock 
      if (dist > dist_old) then
        tcurrent  = tcurrent - h 

!      record the minimum location                                
        if (min_flag) then
          minloc = xx0          
        endif
        min_flag = .false.

      else
!     if the distance is larger thant he last step, update the 
!     ending time of the simulation and the closest point
        tend = tend + h 
        dist_old = dist
      endif


!     if the distance from the z plan is larger than the last 
!     step, set the crossing location and time
      zdist = xx0(3)
      if (abs(zdist) > abs(zdist_old) .and. zmin_flag) then
        zmin_flag = .false.
        zcrossloc = xx0
        izmin = istep
      endif
      zdist_old = zdist

      rtime = rtime + h
      istep = istep + 1
    enddo
    tzcross = tend - h * izmin

    print*,'t0, tcurrent ', t0, tcurrent
    print*,'tend, istep ', tend, istep
    print*,'rtime ', rtime


!
!       set the time to t0
    xx0(4) = -xx0(4)
    xx0(5) = -xx0(5)
    xx0(6) = -xx0(6)


!	now move adjust the test particles from the 
!	second disk to the proper velocity and positions
!
    if (n > n1) then
      do i = n1+1, n
        x0(i,:) = x0(i,:) + xx0(:)
      enddo
    endif
!
! 	include the perturbing galaxy
!
    n = n + 1
    x0(n,:) = xx0

!!!    pscale = 1.1


    return
  end subroutine reset_perturber_position


!---------------------------------------------------------
  subroutine rk41(xx0, h, mass1, mass2, eps1, eps2, xxe, diffq1)

    implicit none
    real(kind=8), dimension(6), intent(in):: xx0
    real(kind=8), intent(in) :: h
    real(kind=8), intent(in) :: mass1, mass2
    real(kind=8), intent(in) :: eps1, eps2
    real(kind=8), dimension(6), intent(out):: xxe

    real(kind=8), dimension(6):: x, f

!---------------------------------------------------
    interface
      subroutine diffq1(x, mass1, mass2, eps1, eps2, f)
        real(kind=8), dimension(6), intent(in) :: x
        real(kind=8), intent(in) :: mass1, mass2
        real(kind=8), intent(in) :: eps1, eps2
        real(kind=8), dimension(6), intent(out):: f
      end subroutine
    end interface

    x = xx0
    call diffq1(x, mass1, mass2, eps1, eps2, f)

    xxe = xx0 + h * f / 6.0d0
    x   = xx0 + h * f / 2.0d0
    call diffq1(x, mass1, mass2, eps1, eps2, f)

    xxe = xxe + h * f / 3.0d0
    x   = xx0 + h * f / 2.0d0
    call diffq1(x, mass1, mass2, eps1, eps2, f)

    xxe = xxe + h * f / 3.0d0
    x   = xx0 + h*f
    call diffq1(x, mass1, mass2, eps1, eps2, f)

    xxe = xxe + h * f /6.0d0

    return
  end subroutine rk41






!
!
!-----------------------------------------------------------
  subroutine diffq_spm(x, mass1, mass2, eps1, eps2, f)

    implicit none
    real(kind=8), dimension(6), intent(in) :: x
    real(kind=8), intent(in) :: mass1, mass2
    real(kind=8), intent(in) :: eps1, eps2
    real(kind=8), dimension(6), intent(out):: f

    real(kind=8) :: r21, r1, a1



    r21=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    r1=sqrt(r21)

    a1 = -mass1 / (r21 + eps1) - mass2 / (r21 + eps2)

    f(1) = x(4)
    f(2) = x(5)
    f(3) = x(6)
    f(4) = a1 * x(1)/r1
    f(5) = a1 * x(2)/r1
    f(6) = a1 * x(3)/r1


    return
  end subroutine diffq_spm





!
!
!-----------------------------------------------------------
  subroutine diffq_nbi(x, mass1, mass2, eps1, eps2, f)

    implicit none
    real(kind=8), dimension(6), intent(in) :: x
    real(kind=8), intent(in) :: mass1, mass2
    real(kind=8), intent(in) :: eps1, eps2
    real(kind=8), dimension(6), intent(out):: f

    real(kind=8) :: r21, r1, a1, a2, at

    real (kind=8) :: c1, c2, c3, v21, v1, xvalue
    real (kind=8) :: sqrtpi

    integer (kind=4) :: ival, ival2
    real (kind=8) :: df_force1, df_force2
    real (kind=8) :: df_sigma, df_rho
    real (kind=8) :: ee1, ee2

!   fix to eliminate a compilation warning message for unused variables
    ee1 = eps1
    ee2 = eps2

    sqrtpi = sqrt(pi)

    r21 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    r1  = sqrt(r21)


! get the index for the calculations
    ival  =  df_index(r1, rrout1   )
    ival2 =  df_index(r1, rrout2   )

! get the forces, sigma and rho, and rescale them

    df_force1 = acceleration_particle(ival) * rs_internal * rs_internal
    df_force2 = acceleration_particle(ival2)* rs_internal * rs_internal

    df_sigma  = new_vr2(ival) * rs_internal * rs_internal
    df_rho    = new_rho(ival) * ( rs_internal * rs_internal * rs_internal )



! interpolated forces 
    a1 = -mass1 * df_force1
    a2 = -mass2 * df_force2
    at = a1 + a2

    write(24,*) r1, ival, df_force1, a1, rrout1


! df
    v21 = x(4)*x(4)+x(5)*x(5)+x(6)*x(6)
    v1  = sqrt(v21)


    xvalue = v1 / df_sigma
    c1 = erf(xvalue) - 2.0d0 * xvalue / SqrtPI * exp(-xvalue*xvalue)
    write(21,*) r1, c1, xvalue, v1, df_sigma, df_rho


! df formula with G=1
    c2 = -4.0d0 * pi * mass2 * lnl / v21
    c3 = c1 * c2 * df_rho    


    f(1) = x(4)
    f(2) = x(5)
    f(3) = x(6)
    f(4) = at * x(1)/r1  - c3 * x(4)/ v1 
    f(5) = at * x(2)/r1  - c3 * x(5)/ v1
    f(6) = at * x(3)/r1  - c3 * x(6)/ v1




    return
  end subroutine diffq_nbi


!
!
!-----------------------------------------------------------
  subroutine diffq_mond(x, mass1, mass2, eps1, eps2, f)

    implicit none
    real(kind=8), dimension(6), intent(in) :: x
    real(kind=8), intent(in) :: mass1, mass2
    real(kind=8), intent(in) :: eps1, eps2
    real(kind=8), dimension(6), intent(out):: f

    real(kind=8) :: r21, r1, a1, tmp, a2



    r21=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    r1=sqrt(r21)

    a1 = -mass1 / (r21 + eps1)
    a2 = -mass2 / (r21 + eps2)
  
    tmp = 2*A0/a1;
    a1 = a1/sqrt(2.0d0) * sqrt(1.0d0 + sqrt(1.0d0 + tmp*tmp));
  
    tmp = 2*A0/a2;
    a2 = a2/sqrt(2.0d0) * sqrt(1.0d0 + sqrt(1.0d0 + tmp*tmp));

    a1 = a1 + a2

    f(1) = x(4)
    f(2) = x(5)
    f(3) = x(6)
    f(4) = a1 * x(1)/r1
    f(5) = a1 * x(2)/r1
    f(6) = a1 * x(3)/r1


    return
  end subroutine diffq_mond




!------------------------------------------------------
!
!
!
  subroutine profile(rin, rout, rscale, nstart, ntot, mass, eps, &
       theta, phi, opt, heat, x0)

!
!
!       variables-
!       opt - option for the distribution
!       rin - inner radius
!       rout - outer radius
!       rscale - scale of brightness drop
!	nstart - start number for placement of particles
!       ntot - number of particles to be placed
!       heat - heat parameter
!       m - mass of galaxy
!       sl - softening length
!       nring - number of rings
!       npart - number of particle per ring (opt)
!       x0 - position of center of mass
!

    implicit none
    real (kind=8), intent(in) :: rin
    real (kind=8), intent(in) :: rout
    real (kind=8), dimension(3), intent(in) :: rscale
    integer (kind=4), intent(in) ::  nstart
    integer (kind=4), intent(in) ::  ntot 
    real (kind=8), intent(in) :: mass 
    real (kind=8), intent(in) :: eps
    real (kind=8), intent(in) :: theta
    real (kind=8), intent(in) :: phi
    integer(kind=4), intent(in) :: opt
    real (kind=8), intent(in) :: heat
!    real (kind=8), intent(in) :: seed
!    real (kind=8), intent(in) :: t0
    real (kind=8), intent(out), dimension(:,:) :: x0

    real (kind=8) :: stheta,ctheta,sphi,cphi
    real (kind=8) :: x3,y3,z3,xv3,yv3,zv3,x2,y2,z2,xv2,yv2,zv2
    real (kind=8) :: x,y,z,xv,yv,zv

    integer (kind=4) :: i, j, n, nprof

    real (kind=8) :: rnorm
    real (kind=8), dimension(:), allocatable :: rp, r, angle, v, p_ring, cp_ring
    real (kind=8) :: st, ct, dr, ran, r1, r2, ptot
    integer (kind=4), dimension(:), allocatable :: n_ring
    integer (kind=4) :: nring, dnring, is, ie, iring, tct


    n = size(x0,1)
    allocate(r(n))
    allocate(angle(n))
    allocate(v(n))

    stheta = sin(theta*pi/180.0d0)
    ctheta = cos(theta*pi/180.0d0)
    sphi   = sin(phi*pi/180.0d0)
    cphi   = cos(phi*pi/180.0d0)

!
!       set up the probablity distribution for the disk



    nprof = 1000
    nring = nprof / 10

    dnring = nprof/nring
    allocate(rp(nprof))
    allocate(n_ring(nprof))
    allocate(p_ring(nprof))
    allocate(cp_ring(nprof))


! set the differential sum of the probability function into a vector
    rnorm = 0.0d0
    dr = (rout - rin)/float(nprof)
    do i = 1, nprof
      r1 = float(i)*dr + rin
      rp(i) =  distrb(r1, opt, rscale ) * r1 * dr * 2.0d0 * pi
      rnorm = rnorm + rp(i)
    enddo

! normalize the vector
!
    rp = rp / rnorm

!  take the fine bins and put them into the selection bins
!
    tct = 0
    do iring =  1, nring
      is = (iring -1) * dnring + 1
      ie = (iring) * dnring 

      ptot = 0.0d0
      do i = is, ie
        ptot = ptot + rp(i) 
      enddo
      p_ring(iring) = ptot
    enddo

! formulative cumulative distribution function
!
    cp_ring(1) = p_ring(1)
    do iring = 2, nring
      cp_ring(iring) = cp_ring(iring -1) + p_ring(iring)

    enddo


! find the number of particles in each bin
!
    n_ring = 0
    do i = nstart, ntot     
! find the radial position bin
      ran = randm()
      j = 1
      do while (ran > cp_ring(j) .and. j < nring) 
        j = j + 1
      enddo
      n_ring(j) = n_ring(j) + 1
    enddo

    tct = 0
    i = nstart
    do iring =  1, nring
      is = (iring -1) * dnring + 1
      ie = (iring) * dnring 

      r1 = float(is)*dr + rin
      r2 = float(ie)*dr + rin

      do j = 1, n_ring(iring)
        ran = randm()
        r(i) = r1 + ran * (r2 - r1)
        i = i + 1
      enddo
    enddo


!   set the angular positions and orbital velocities
!
    do i=nstart,ntot
      angle(i) = 2.0d0 * pi * randm()

      v(i) = circular_velocity(mass,r(i),rout,eps)

    enddo
!    stop



! set position and velocity based on the distribution parameters
    do i= nstart,ntot

      st  =  sin(angle(i))
      ct  =  cos(angle(i))

      x   =  ct*r(i)
      y   =  st*r(i)
      z   =  0.0

      xv  = -v(i)*st
      yv  =  v(i)*ct
      zv  =  0.0

      x2  =   x * ctheta +  z * stheta
      y2  =   y
      z2  =  -x * stheta +  z * ctheta
      xv2 =  xv * ctheta + zv * stheta
      yv2 =  yv
      zv2 = -xv * stheta + zv * ctheta

      x3  =  x2  * cphi -  y2 * sphi
      y3  =  x2  * sphi +  y2 * cphi
      z3  =  z2
      xv3 =  xv2 * cphi - yv2 * sphi
      yv3 =  xv2 * sphi + yv2 * cphi
      zv3 =  zv2


!!    write(16,*) r(i)

      x0(i,1) = x3 
      x0(i,2) = y3 
      x0(i,3) = z3 
      x0(i,4) = xv3  + randm()*heat
      x0(i,5) = yv3  + randm()*heat
      x0(i,6) = zv3  + randm()*heat
!    x0(6,i) = t0

    enddo

    deallocate(r)
    deallocate(angle)
    deallocate(v)
    deallocate(rp)
    deallocate(n_ring)
    deallocate(p_ring)
    deallocate(cp_ring)

    return
  end subroutine profile
!
!
!-----------------------------------------------------------------
!
!
  function distrb(r1,opt,rscale)
    real (kind=8) :: distrb
    real(kind=8) :: r1, rscale(3)
    integer opt


    if (opt.eq.1) then
      distrb = 1./r1
    elseif (opt.eq.2) then
      distrb = exp(-r1/rscale(1))
    elseif (opt.eq.3) then
      distrb = exp(-r1*r1*rscale(1) - rscale(2)*r1 - rscale(3) )
    endif

    return
  end function distrb

!
!
!-----------------------------------------------------------------
!
!
  function randm()
!seed)
    real (kind=8) :: randm
!    real (kind=8) seed

    call random_number(randm )
    return
  end function randm

!
!
!-----------------------------------------------------------------
!
!
  function circular_velocity(mass, r, rout, eps)
    real (kind=8) :: circular_velocity 
    real (kind=8) :: r, eps, rout, mass
    real (kind=8) :: ftotal,tmp
    integer (kind=4) :: ival

    if(potential_type .eq. 0) then
      ftotal = mass / ( r*r + eps )
    else if(potential_type .eq. 1) then
      ival = df_index(r, rout)
      ftotal = mass * acceleration_particle(ival) * rs_internal * rs_internal
    else if(potential_type .eq. 2) then
      ftotal = mass / ( r*r + eps )
      tmp = 2*A0/ftotal;
      ftotal = ftotal/sqrt(2.0d0) * sqrt(1.0d0 + sqrt(1.0d0 + tmp*tmp));
    endif

    circular_velocity = sqrt(ftotal * r)

    return

  end function circular_velocity

!---------------------------------------------------
  subroutine set_perturber_position(pos, vel, t0, x0, n1, n )

    implicit none

    real (kind=8), dimension(3) :: pos, vel
    real (kind=8), intent(in) :: t0
    integer (kind=4), intent(inout) :: n
    integer (kind=4), intent(in) :: n1
    real (kind=8), dimension(:,:), intent(inout) :: x0

    real (kind=8), dimension(6) :: xx0  !, xxe
    integer(kind=4) :: i
    real (kind=8) :: tcurrent
    

    xx0(1:3) = pos
    xx0(4:6) = vel
    tcurrent = t0


!	now move adjust the test particles from the 
!	second disk to the proper velocity and positions
!
    if (n > n1) then
      do i = n1+1, n
        x0(i,:) = x0(i,:) + xx0(:)
      enddo
    endif
!
! 	include the perturbing galaxy
!
    n = n + 1
    x0(n,:) = xx0


    return
  end subroutine set_perturber_position


  function rvToCoe(r, v, mu)
    implicit none

    real (kind=8), intent(in), dimension(3) :: r, v
    real (kind=8), intent(in) :: mu
    real (kind=8), dimension(7) :: rvToCoe
    
    real (kind=8), dimension(3) :: h, n, v1, v2, ev
    real (kind=8), dimension(3) :: K = (/0.0,0.0,1.0/)
    real (kind=8) :: muInv, rmag, vmag, hmag, nmag, tmp1, tmp2, p, ecc, cosi, cosO, cosw, cosv, cosu, arg
	
    muInv = 1.0/mu

    rmag = mag(r)
    vmag = mag(v)
    
    h = cross(r,v)
    hmag = mag(h)
    
    n = cross(K,h)
    nmag = mag(n)
    
    tmp1 = vmag*vmag - mu/rmag
    tmp2 = dot(r,v)
    
    v1 = scaleVec(tmp1,r)
    v2 = scaleVec(tmp2,v)
    
    ev = sub(v1,v2)
    ev = scaleVec(muInv,ev)
    
    p = hmag*hmag*muInv
    ecc = mag(ev)
    cosi = h(3)/hmag
    cosO = n(1)/nmag
    cosw = dot(n,ev)/(nmag*ecc)
    cosv = dot(ev,r)/(ecc*rmag)
    cosu = dot(n,r)/(nmag*rmag)

    rvToCoe(1) = p
    rvToCoe(2) = ecc
    rvToCoe(3) = acos(cosi)
    
    tmp1 = acos(cosO)
    
    if (n(1) < 0) then
      tmp1 = 2.0*pi-tmp1
    endif
    
    rvToCoe(4) = tmp1
        
    tmp1 = acos(cosw)
    
    if (ev(2) < 0 ) then
      tmp1 = 2.0*pi-tmp1
    endif
    
    rvToCoe(5) = tmp1
    
    tmp1 = acos(cosv)
    
    if (dot(r,v) < 0) then
      tmp1 = 2.0*pi-tmp1
    endif
    
    rvToCoe(6) = tmp1
        
    if(cosu > 1.0 .or. cosu < -1.0) then
      arg = -1
      rvToCoe(7) = sqrt(arg) !Double.NaN IN CODE RELEASE VERSION
    else
      tmp1 = acos(cosu)
      if(r(2)>0) then
        tmp1 = 2.0*pi-tmp1 ! NOT IN CODE RELEASE VERSION
      endif 
      rvToCoe(7) = tmp1
    endif
    
    return
  end function rvToCoe


    function dot(v1, v2)
      real (kind=8), intent(in), dimension(3) :: v1, v2
      integer(kind=4) :: i, length
      real (kind=8) :: cp, dot
    	
      length = size(v1)
      cp = 0.0;
    	
      do i=1,length
        cp = cp + v1(i)*v2(i)
      enddo
    	
      dot = cp
      return
    end function dot
   

    function mag(v)
      real (kind=8), intent(in), dimension(:) :: v
      real (kind=8) :: mag
      mag = sqrt(dot(v,v))
      return
    end function mag
    
	
    function cross(v1, v2)
      real (kind=8), intent(in), dimension(3) :: v1, v2
      real (kind=8), dimension(3) :: cross
		
      cross(1) = v1(2)*v2(3)-v1(3)*v2(2)
      cross(2) = v1(3)*v2(1)-v1(1)*v2(3)
      cross(3) = v1(1)*v2(2)-v1(2)*v2(1)
    	
      return
    end function cross
    
    function scaleVec(sc, v1)
      real (kind=8), intent(in), dimension(:) :: v1
      real (kind=8), intent(in) :: sc
      integer(kind=4) :: i, length
      real (kind=8), allocatable, dimension(:) :: scaleVec
      length = size(v1)
     
      allocate(scaleVec(length))
      scaleVec(1:length) = sc*v1(1:length)
        
      return
    end function scaleVec
    
    function sub(v1,v2)
      real (kind=8), intent(in), dimension(:) :: v1, v2
      integer(kind=4) :: i, length
      real (kind=8), allocatable, dimension(:) :: sub
      length = size(v1)
      allocate(sub(length))
      sub = v1 - v2
    end function sub
    
    
    
	! * Find the time of rmin, assuming earlier than now, given
	! * the r and v values.  Returns r and v at time of rmin
	! * by replacing r and v.  r and v are given as
	! * {rx,ry,rz,vx,vy,vz}.
	
    function getTStart(rv, mass1, mass2, eps1, eps2, h, tmin, mind, rout1, rout2)
      real (kind=8), intent(inout), dimension(7) :: rv
      real (kind=8), intent(in) :: mass1, mass2, eps1, eps2, h, tmin, mind, rout1, rout2

      real (kind=8) :: t, distOld, distNew, mu, ecc, a, period, apocenter, a2, tApp, distNearApp
      real (kind=8) :: minDist, minVel
      real (kind=8), dimension(3) :: r, v
      real (kind=8), dimension(4) :: outStuff, getTStart
      real (kind=8), dimension(7) :: coe, xxe, rvmin
      logical :: isEllipse
		
      integer (kind=4) :: i, length

      mu = mass1+mass2
      t = 0d0
		
      r = (/rv(1), rv(2), rv(3)/)
      v = (/-rv(4),-rv(5), -rv(6)/)
      coe = rvToCoe(r,v,mu)
      ecc = coe(2)
      a = coe(1)/(1.0-ecc*ecc)
      period = 0.0d0
      apocenter = a*(1+ecc)
      a2 = apocenter*apocenter
      tApp = 0.0d0
    	
      isEllipse = .false.
    	    	
      if (ecc < 1.0) then
        isEllipse = .true.
        period = 2d0*pi/sqrt(mu)*(a**1.5d0)
        period = period * 1.0d0
      endif
    	
      rvmin = rv ! should assign rvmin the values of rv
      
      distNew = rv(1)*rv(1)+ rv(2)*rv(2) + rv(3)*rv(3)
      distOld = 2.0*distNew
    	
      distNearApp = -1e30
   
      ! keep looping as long as distance is decreasing
      do while(tmin < t) 
        r = (/rv(1), rv(2), rv(3)/)
        v = (/rv(4), rv(5), rv(6)/)
        coe = rvToCoe(r,v,mu)
        xxe(7)=t +h
        call wrap_rk41(rv, h, mass1, mass2, eps1, eps2, xxe)
    	   
        distNew = xxe(1)*xxe(1)+ xxe(2)*xxe(2) + xxe(3)*xxe(3)
    	    
        ! if it's ellipse and it's near apocenter, take this time
        if ( isEllipse .and. (abs(distNew-a2)/a2 < 0.05d0) ) then
          if(distNew > distNearApp) then
            distNearApp = distNew
            tApp = t
          endif
        endif
    	    
        if (distNew < distOld) then
          distOld = distNew
          rvmin(1:7) = xxe(1:7)
          rvmin(7) = rvmin(7) - h
        endif

        rv(1:7) = xxe(1:7)
        rv(7) = xxe(7) - h * 2.0d0
        t = t - h
      enddo

      rv(1:7) = rvmin(1:7)
    	
      minDist = sqrt(rv(1)*rv(1) + rv(2)*rv(2) + rv(3)*rv(3))
      minVel =  sqrt(rv(4)*rv(4) + rv(5)*rv(5) + rv(6)*rv(6))
    	
      t = rv(7)
    	
      if(isEllipse .and. tApp < 0.0d0) then
        t = tApp
      else
        t = t - mind/minVel;
      endif
    	
      outStuff = (/t,minDist,minVel,rv(7)/)
    	
      getTStart = outStuff
    	
      return
    end function getTStart
    
end module setup_module


