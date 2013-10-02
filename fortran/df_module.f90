!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
module df_module
!     
!
!     -----Description: dynamical friction module
!
!
!
!     ----------------------------------------------------------------
!
  use parameters_module
  implicit none
!
!     -----Variable declarations
!          ---------------------
!
    integer (kind=4), parameter :: nnn = 10000
    real (kind=8), dimension(nnn) :: rad    
    real (kind=8), dimension(nnn) :: rho_halo, mass_halo
    real (kind=8), dimension(nnn) :: rho_disk, mass_disk
    real (kind=8), dimension(nnn) :: rho_bulge, mass_bulge
    real (kind=8), dimension(nnn) :: rho_total, mass_total
    real (kind=8), dimension(nnn) :: masses, radius, density
    real (kind=8), dimension(nnn) :: vr2, vr, new_vr2, new_vr
    real (kind=8), dimension(nnn) :: acceleration, acceleration_particle
    real (kind=8), dimension(nnn) :: new_mass, new_rho, phi

    real (kind=8), parameter :: rs_internal = 10.0d0

    real (kind=8) :: pscale

    real (kind=8) :: lnl
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
  subroutine init_distribution
!     
!
!     -----Description: initializes the distribution 
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


!
!     ----------------------------------------------------------------
!
    real (kind=8) :: rmax
    real (kind=8) :: mold, dmold, mtot
    real (kind=8) :: rscale
    real (kind=8) :: dx, x
    real (kind=8) :: alphahalo, qhalo, gammahalo, mhalo, rchalo, rhalo, epsilon_halo
    real (kind=8) :: zdisk, hdisk, zdiskmax
    real (kind=8) :: hbulge, mbulge
    real (kind=8) :: rho_tmp
    real (kind=8) :: G, factor
    real (kind=8) :: r, m, sqrtpi
    real (kind=8) :: p1, rd, rho_local
    real (kind=8) :: p, rr, dr, rh, dp, mnew, dm
    real (kind=8) :: acc_merge, rad_merge, acc

    integer (kind=4) :: j, nmax, k, nmerge, ntotal, jj



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3

!!!!!
! set the constant for dynamical friction
    lnl = 0.00d0
! default for Merger Zoo
    lnl = 0.001d0

!!!!!
! set up the parameters for the halo
    mhalo = 5.8d0  
    rhalo = 10.0d0
    rchalo = 10.0d0
    gammahalo = 1.0d0
    epsilon_halo = 0.4d0 
    SqrtPI = sqrt(pi)

!!!!!
! derive additional constants
    qhalo = gammahalo / rchalo
    alphahalo = 1.0d0 / ( 1.0d0 - SqrtPI * qhalo * exp(qhalo**2) * (1.0d0 - erf(qhalo)) )

!!!!!
! set the integration limits and zero integration constants
    rmax = 20
    nmax = 2000
    dr = rmax / (nmax)
    mold = 0

    rscale = 5
!    ntotal = nmax * rscale
    ntotal = nnn

!!!!!
! set the limits for integration, and zero integration constants
    k = nmax / 2
    dx = 1.0  / k
    x = 0.0d0
    dmold = 0.0d0
    mtot = 0.0d0
    rad = 0.0d0
    m = 0.0d0
    G = 1
    

!!!!!
! set the fundamental disk parameters
    zdisk = 0.2
    zdiskmax = 3.5
    hdisk = 1.0
    

!!!!!
! set the fundamental bulge parameters
    hbulge = 0.2
    mbulge = 0.3333
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
!!!!! set up the radius array
    do j = 1, nmax
      x = x + dx
      rad(j)= x * rchalo
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
!!!!!
 
   dr = rad(2) - rad(1)
   dx = dr / rchalo 
   
   do j = 1, nmax

! set the position
      r = rad(j)
      x = r / rchalo

! calculate the local rho based 
      rho_tmp = alphahalo / (2*SqrtPI**3 ) * (exp(-x**2) / (x**2 + qhalo**2))

! renormalize for the new halo size
      rho_tmp = rho_tmp / ( rchalo * rchalo * rchalo) 

! calculate mass in local shell, and update total mass
!      dm = rho_tmp * 4 * pi * x * x *dx
      dm = rho_tmp * 4 * pi * r * r *dr
      mtot = mtot + dm

! store values in an array
      rho_halo(j) = rho_tmp * mhalo 
      mass_halo(j) = mtot * mhalo
    end do

!!!!!
! now calculate the potential
    do j = 1, nmax
      r = rad(j)
      m = mass_halo(j)
      p1 = -G * m / r
      phi(j) = p1

    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
! disk model

!!!!!
! loop over the distribution
    do j = 1,nmax

! set the radius
      rd = rad(j)
  
! find the local density in the disk
      rho_local  = exp(-rd/hdisk)/ (8*pi*hdisk**2.0d0) 
      rho_disk(j) = rho_local
      
! find the mass in the spherical shell
      mnew = 4 * pi * rho_local * rd *rd * dr
      
      mass_disk(j) = mnew + mold
      mold = mass_disk(j)
    end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
! bulge model

!!!!!
! loop over the distribution
    mold = 0.0
    do j = 1,nmax


! set the radius
      rd = rad(j)
  
! find the local density in the disk
      rho_local  = exp(-rd**2/hbulge**2)
      rho_bulge(j) = rho_local

! find the mass in the spherical shell
      mnew = 4 * pi * rho_local * rd *rd * dr
      
      mass_bulge(j) = mnew + mold
      mold = mass_bulge(j)
    end do

! renormalize distribution
    factor = mbulge / mass_bulge(nmax)
    do j = 1,nmax
      mass_bulge(j) = mass_bulge(j) * factor
      rho_bulge(j)  = rho_bulge(j)  * factor
    end do

  
    dr = rad(2) - rad(1)      
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
    j = 1
    mass_total(j)=  (mass_halo(j) + mass_disk(j) + mass_bulge(j)) 
    r = rad(j)
    rho_total(j) = mass_total(j) /  (4.0d0/3.0d0 * pi * r * r * r)
    dr = rad(2) - rad(1)

    do j = 2,nmax
      r = rad(j)
      mass_total(j)=  (mass_halo(j) + mass_disk(j) + mass_bulge(j)) 

      dm = mass_total(j) - mass_total(j-1)
      rho_total(j) = dm / (4 * pi * r * r * dr)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the velocity dispersion v_r**2

    masses = mass_total
    radius = rad
    density = rho_total


    do j = 1,nmax

      p = 0.0d0
      rr = radius(j)
      dr = radius(nmax) / nmax
      do jj = j,nmax
      
        m  = masses(jj)
        rh = density(jj)
        rr = rr + dr
        
        dp = rh * G * m / rr**2 * dr
        p = p + dp
      end do
      
      vr2(j) = 1/density(j) * p
      vr(j) = sqrt(vr2(j))
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the velocity dispersion v_r**2
    masses = mass_total
    radius = rad
    density = rho_total
    
    
    do j = 1,nmax

      p = 0.0d0
      rr = radius(j)
      dr = radius(nmax) / nmax
      do jj = j,nmax
      
      m  = masses(jj)
      rh = density(jj)
      rr = rr + dr
      
      dp = rh * G * m / rr**2 * dr
      p = p + dp
    end do
  
    vr2(j) = 1/density(j) * p
    vr(j) = sqrt(vr2(j))
  enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the accelerations felt by the particles and center of mass
    masses = mass_total
    radius = rad
    density = rho_total


    do j = 1,nmax
      rr = radius(j)
      m  = masses(j)
      acceleration(j) = G * m / rr**2
    end do


    acceleration_particle = acceleration
    nmerge = 50
    acc_merge = acceleration(nmerge)
    rad_merge = rad(nmerge)
    
    do j = 1, nmerge
      rr = radius(j)
      m  = masses(j)

! smoothed acceleration
      acc = G * m / (rr**2 + .1* (rad_merge -rr)) 
      acceleration_particle(j) = acc
      
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rederive the masses from the new particle acceleration
    radius = rad
    dr = rad(2) - rad(1)

! find the accelerations felt by the particles and center of mass
    radius = rad
    
    do j = 1, nmax
      rr = radius(j)
      new_mass(j) = rr**2 * acceleration_particle(j) / G
      new_rho(j)  = new_mass(j) / (4 * pi * rr * rr * dr)
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the velocity dispersion v_r**2 using the new density and masses
    masses = new_mass
    radius = rad
    density = new_rho
    

    do j = 1, nmax

      p = 0.0d0
      rr = radius(j)
      dr = radius(nmax) / nmax
      do jj = j,nmax
      
      m  = masses(jj)
      rh = density(jj)
      rr = rr + dr
      
      dp = rh * G * m / rr**2 * dr
      p = p + dp
    end do
    
    new_vr2(j) = 1/density(j) * p
    new_vr(j) = sqrt(new_vr2(j))
    
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! extend the values to large rmax
do  j= nmax+1, ntotal
  mass_total(j) = mass_total(nmax)
  mass_halo(j) = mass_halo(nmax)
  mass_disk(j) = mass_disk(nmax)
  mass_bulge(j) = mass_bulge(nmax)
  new_mass(j) = new_mass(nmax)

!  rho_total(j) = 1e-3
!  new_rho(j) = new_rho(nmax)
  rho_total(j) = 0.0d0
  new_rho(j)   = 0.0d0

  vr(j)      = 1d-6
  vr2(j)     = 1d-6
  new_vr(j)  = 1d-6
  new_vr2(j) = 1d-6

  m = mass_total(nmax)
  rr = rad(nmax) + dr*(j - nmax)
  rad(j) = rr
  acc = G * m / rr**2  
  acceleration_particle(j) = acc
  acceleration(j) = acc

end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! normalize to the unit mass

do  j= 1, ntotal
  mass_total(j)  = mass_total(j) / 7.13333d0
  mass_halo(j)   = mass_halo(j)  / 7.13333d0
  mass_disk(j)   = mass_disk(j)  / 7.13333d0
  mass_bulge(j)  = mass_bulge(j) / 7.13333d0
  new_mass(j)    = new_mass(j)   / 7.13333d0

  rho_total(j)   = rho_total(j)  / 7.13333d0
  new_rho(j)     = new_rho(j)    / 7.13333d0
 
  vr(j)          = vr(j)      / 7.13333d0
  vr2(j)         = vr2(j)     / 7.13333d0
  new_vr(j)      = new_vr(j)  / 7.13333d0
  new_vr2(j)     = new_vr2(j) / 7.13333d0

  rad(j)         = rad(j) 

  acceleration_particle(j) = acceleration_particle(j) / 7.13333d0
  acceleration(j)          = acceleration(j)  / 7.13333d0

!!  write(11,*) rad(j), new_rho(j), new_mass(j),  new_vr(j)
end do



pscale = 1.0d0



!% ! tabulate the right hand side of the dynamical friction formula
!% xmax = 10
!% x = 0.0d0
!% dx = xmax / j
!% do j = 1,j
!%   x = x + dx
!%   rhs(j) = erf(x) - 2.0*x / sqrt(pi) * exp(-x**2)
!%   xx(j) = x
!% end do




  end subroutine init_distribution



!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
  integer (kind=4) function df_index(rin, rs)
!     
!
!     -----Description: interpolates the force for a given particle
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
    real (kind=8), intent(in) :: rin
    real (kind=8), intent(in) :: rs
    integer (kind=4) :: ival
    real (kind=8), parameter :: rmax_scale = 100.0d0

!
!     ----------------------------------------------------------------
!

    
!    ival = min(int(  (rin * rs_internal/ rmax_scale) * nnn + 1), nnn) 
    ival = min(int(  (rin * pscale * rs_internal/ rmax_scale) * nnn + 1), nnn) 
    df_index = ival

!!    print*,rin, pscale, rs_internal, rmax_scale, nnn, 'dddd'

    return

  end function df_index






  




end module df_module


