module integrator

  use parameters_module
  use df_module
 
  real (kind=8), dimension(:,:), allocatable, target:: xe
  real (kind=8), dimension(:,:), allocatable :: x, f
  real (kind=8), dimension(:), allocatable :: r22, r21, r2n
  real (kind=8), dimension(:), allocatable :: r1, r2, rn
  real (kind=8), dimension(:), allocatable :: a1, a2, a3
  real (kind=8), dimension(:), allocatable :: mond_tmp

  real (kind=8) :: m1, m2, m3
  real (kind=8) :: eps1, eps2
  real (kind=8) :: theta_i1, phi_i1, theta_i2, phi_i2
  real (kind=8), dimension(3):: rscale_i1, rscale_i2
  real (kind=8) :: rrout1, rrout2

  real (kind=8), dimension(:), allocatable :: df_force11, df_force22, df_forcen, c3n
  integer (kind=4), dimension(:), allocatable :: ival11, ival22, ivaln

  integer (kind=4) :: pn, pn1, pn2

contains

!-----------------------------------------------------------------
!
subroutine init_rkvar(x0, mass1, mass2, epsilon1, epsilon2, t1, p1, t2, p2, rs1, rs2, ro1, ro2, nn, nn1, nn2)
  implicit none
  
  real (kind=8), dimension(:,:), intent(in) :: x0
  real (kind=8), intent(in) :: mass1, mass2, epsilon1, epsilon2, p1, t1, p2, t2
  real (kind=8), intent(in), dimension(3):: rs1, rs2
  real (kind=8), intent(in) :: ro1, ro2
  integer (kind=4), intent(in) :: nn, nn1, nn2
  integer (kind=4) :: n


  pn = nn
  pn1 = nn1
  pn2 = nn2

  n = size(x0, 1)

  allocate(x(n,6))
  allocate(f(n,6))
  allocate(xe(n,6))

  allocate(r22(n))
  allocate(r21(n))
  allocate(r2n(n))
  allocate(r1(n))
  allocate(r2(n))
  allocate(rn(n))
  allocate(a1(n))
  allocate(a2(n))
  allocate(a3(n))
  allocate(mond_tmp(n))


  m1 = mass1
  m2 = mass2
  m3 = mass2
  eps1 = epsilon1*epsilon1
  eps2 = epsilon2*epsilon2
  


!
!     -----allocate space for helper arrays
!

  allocate(ival11(n))
  allocate(ival22(n))
  allocate(ivaln(n))
  allocate(df_force11(n))
  allocate(df_force22(n))
  allocate(df_forcen(n))
  allocate(c3n(n))


  phi_i1 = p1
  theta_i1 = t1
  phi_i2 = p2
  theta_i2 =t2

  rscale_i1 = rs1
  rscale_i2 = rs2

  rrout1 = ro1
  rrout2 = ro2

  return
end subroutine init_rkvar

!-----------------------------------------------------------------
!
subroutine deallocate_rkvar

  deallocate(x)
  deallocate(f)
  deallocate(xe)

  deallocate(r22)
  deallocate(r21)
  deallocate(r2n)
  deallocate(r1)
  deallocate(r2)
  deallocate(rn)
  deallocate(a1)
  deallocate(a2)  
  deallocate(a3)
  deallocate(mond_tmp)


!
!     -----deallocate scratch space
!


  return
end subroutine deallocate_rkvar


!---------------------------------------------------
! Use this method so that neither the caller of rk4
! nor implementation need to know which potential
! is being used
!---------------------------------------------------
subroutine  wrap_rk4(x0, h, xout)
    implicit none
    real (kind=8), dimension(:,:), intent(in) :: x0
    real (kind=8) :: h
    real (kind=8), dimension(:,:), intent(inout) :: xout

    if ( potential_type .eq. 0) then
      call rk4(x0, h, xout, diffeq_spm)
    else if ( potential_type .eq. 1) then
      call rk4(x0, h, xout, diffeq_nbi)
    else if ( potential_type .eq. 2) then
      call rk4(x0, h, xout, diffeq_mond)
    endif

    return

end subroutine wrap_rk4

!-----------------------------------------------------------------
!
subroutine rk4(x0, h, xout, diffeq)
  
  implicit none

  real (kind=8), dimension(:,:), intent(in) :: x0
  real (kind=8) :: h
  real (kind=8), dimension(:,:), intent(inout) :: xout

  integer (kind=4) ::  n

! Define an interface so that we can integrate any diffeq
  interface
    subroutine diffeq(x)
      real (kind=8), dimension(:,:), intent(in) :: x  
    end subroutine 
  end interface

  n = size(x0,1)
  x = x0

  call diffeq(x)
  xe = x0 + h * f / 6.0d0
  x  = x0 + h * f / 2.0d0

  call diffeq(x)
  xe = xe + h * f / 3.0d0
  x  = x0 + h * f / 2.0d0

  call diffeq(x)
  xe = xe + h * f / 3.0d0
  x  = x0 + h * f

  call diffeq(x)
  xe = xe + h * f / 6.0d0

  xout = xe
  
return
end subroutine rk4


!
!
!----------------------------------------------------------------
subroutine diffeq_spm(x)

  implicit none

  real (kind=8), dimension(:,:), intent(in) :: x  
  real (kind=8), dimension(6) ::xn
  integer (kind=4) :: n

  n = size(x,1)
  xn = x(n,:)

  r22 = (x(:,1)-xn(1))**2  +(x(:,2)-xn(2))**2 + (x(:,3)-xn(3))**2
  r21 = x(:,1)**2  + x(:,2)**2  + x(:,3)**2
  r2n = xn(1)**2 + xn(2)**2 + xn(3)**2

  r2 = sqrt(r22)
  r1 = sqrt(r21)
  rn = sqrt(r2n)

  a1 = -m1 / (r21 + eps1)
  a2 = -m2 / (r22 + eps2)
  a3 = -m3 / (r2n + eps2)

  ! this is a correction to prevent NaN errors in the vectorized
  ! function evalution at the location of the second mass
  r2(n) = 1.0d0

  ! calculate the RHS of the diffeq

  f(:,1) = x(:,4)
  f(:,2) = x(:,5)
  f(:,3) = x(:,6)

  f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn
  f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn
  f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn


! i = 900
! write(20,*) x(i,1), x(i,2), x(i,3)
! write(20,*) a1(i)*x(i,1)/r1(i), a1(i)*x(i,2)/r1(i), a1(i)*x(i,3)/r1(i)
! write(20,*)  a2(i)*(x(i,1)-xn(1))/r1(i), a2(i)*(x(i,2)-xn(2))/r1(i), &
!           a2(i)*(x(i,3)-xn(3))/r1(i)
! write(20,*)  a3(i)*(xn(1))/r1(i), a3(i)*(xn(2))/r1(i), &
!           a3(i)*(xn(3))/r1(i) 
!  stop

  return
end subroutine diffeq_spm

!----------------------------------------------------------------




!----------------------------------------------------------------
subroutine diffeq_nbi(x)

  implicit none

  real (kind=8), dimension(:,:), intent(in) :: x  
  real (kind=8), dimension(6) ::xn
  integer (kind=4) :: n
  integer (kind=4) :: i

!  real (kind=8) :: df_forcen
  real (kind=8) :: df_sigma, df_rho
  real (kind=8) :: c1, c2, xvalue, v1, v21
  real (kind=8) :: sqrtpi
    
  sqrtpi = sqrt(pi)

  n = size(x,1)
  xn = x(n,:)

! distance between the main galaxy and the particle
  r21 = x(:,1)**2  + x(:,2)**2  + x(:,3)**2
  r1 = sqrt(r21)

! distance between the companion and the particle
  r22 = (x(:,1)-xn(1))**2  +(x(:,2)-xn(2))**2 + (x(:,3)-xn(3))**2
  r2 = sqrt(r22)


! distance between the two galaxies - the tidal force
  r2n = xn(1)**2 + xn(2)**2 + xn(3)**2
  rn = sqrt(r2n)
    
  
  do i = 1, n
    ival11(i) =  df_index(r1(i), rrout1)
  enddo

  do i = 1, n
    ival22(i) =  df_index(r2(i), rrout2)
  enddo

  do i = 1, n
    ivaln(i)  =  df_index(rn(i), rrout2)
  enddo


  df_force11 = acceleration_particle(ival11) * rs_internal * rs_internal
  df_force22 = acceleration_particle(ival22) * rs_internal * rs_internal
  df_forcen  = acceleration_particle(ivaln)  * rs_internal * rs_internal



! get the forces, sigma and rho, and rescale them
  df_sigma  =  new_vr2(ivaln(1)) * rs_internal * rs_internal
  df_rho    =  new_rho(ivaln(1)) * ( rs_internal * rs_internal * rs_internal )


! interpolated forces 
  a1 = -m1 * df_force11
  a2 = -m2 * df_force22
  a3 = -m3 * df_forcen

!  write(25,*) r1(1000), ival11(1000), df_force11(1000), a1(1000), rrout1
!  print*, df_index(5.3d0, rrout1), '  tests '
!  stop

! df
    v21 = xn(4)*xn(4)+xn(5)*xn(5)+xn(6)*xn(6)
    v1  = sqrt(v21)


    xvalue = v1 / df_sigma
    c1 = erf(xvalue) - 2.0d0 * xvalue / SqrtPI * exp(-xvalue*xvalue)

! df formula with G=1
    c2  = 4.0d0 * pi * m2 * lnl / v21
    c3n(1:n-1) = 0.0d0
    c3n(pn1+1:n)  = c1 * c2 * df_rho

!    c3 = 0.0d0

  ! this is a correction to prevent NaN errors in the vectorized
  ! function evalution at the location of the second mass
  r2(n) = 1.0d0


  ! calculate the RHS of the diffeq

  f(:,1) = x(:,4)
  f(:,2) = x(:,5)
  f(:,3) = x(:,6)

  f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn - c3n * xn(4)/ v1 
  f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn - c3n * xn(5)/ v1
  f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn - c3n * xn(6)/ v1




  return
end subroutine diffeq_nbi

!
!
!----------------------------------------------------------------
subroutine diffeq_mond(x)
  implicit none

  real (kind=8), dimension(:,:), intent(in) :: x  
  real (kind=8), dimension(6) ::xn
  integer (kind=4) :: n

  n = size(x,1)
  xn = x(n,:)

  r22 = (x(:,1)-xn(1))**2  +(x(:,2)-xn(2))**2 + (x(:,3)-xn(3))**2
  r21 = x(:,1)**2  + x(:,2)**2  + x(:,3)**2
  r2n = xn(1)**2 + xn(2)**2 + xn(3)**2

  r2 = sqrt(r22)
  r1 = sqrt(r21)
  rn = sqrt(r2n)

  a1 = -m1 / (r21 + eps1)
  a2 = -m2 / (r22 + eps2)
  a3 = -m3 / (r2n + eps2)

  ! this is a correction to prevent NaN errors in the vectorized
  ! function evalution at the location of the second mass
  r2(n) = 1.0d0

  ! scale the accelerations to reflect mond
  mond_tmp = 2*A0/a1;
  a1 = a1/sqrt(2.0d0) * sqrt(1.0d0 + sqrt(1.0d0 + mond_tmp*mond_tmp));

  mond_tmp = 2*A0/a2;
  a2 = a2/sqrt(2.0d0) * sqrt(1.0d0 + sqrt(1.0d0 + mond_tmp*mond_tmp));

  mond_tmp = 2*A0/a3;
  a3 = a3/sqrt(2.0d0) * sqrt(1.0d0 + sqrt(1.0d0 + mond_tmp*mond_tmp));


  n = size(x,1)
  xn = x(n,:)

  ! calculate the RHS of the diffeq

  f(:,1) = x(:,4)
  f(:,2) = x(:,5)
  f(:,3) = x(:,6)

  f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn
  f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn
  f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn


  return
end subroutine diffeq_mond

end module integrator
