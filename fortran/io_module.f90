 module io_module


  use parameters_module

  implicit none

contains


!----------------------------------------------------------------
!
  subroutine input_particles(unit, n, x0)

    implicit none
    integer (kind=4), intent(in) :: unit
    integer (kind=4), intent(in) :: n
    real (kind=8), dimension(:,:), intent(inout) :: x0
    integer (kind=4) :: i

    do i = 1, n
      read(unit,'(6e16.8)')  x0(i,:)
    enddo

    return
  end subroutine input_particles

!----------------------------------------------------------------
!
  subroutine output_particles(unit, x0, mass1, mass2, eps1, eps2, n, n1, n2, time, header_on)

    implicit none
    integer (kind=4), intent(in) :: unit
    real (kind=8), dimension(:,:), intent(inout) :: x0
    integer (kind=4), intent(in) :: n, n1, n2
    real (kind=8), intent(in) :: mass1, mass2, eps1, eps2, time
    logical, intent(in) :: header_on
    integer (kind=4) :: i
    real (kind=8), dimension(6) :: cm_sys


   if( header_on) then
     write(unit, '(e16.8)') time
     write(unit, '(2e16.8)') mass1, mass2
     write(unit, '(2e16.8)') eps1, eps2
     write(unit, '(3i8)') n, n1, n2
   endif


   cm_sys = mass2 * x0(n,:) / (mass1 + mass2)
   cm_sys = 0.0d0

   do i = 1, n
     write(unit,'(6e16.8)')  x0(i,:) - cm_sys
   enddo




   return
  end subroutine output_particles
!----------------------------------------------------------------
!
  subroutine create_gnuplot_script(x0, iout)

    real (kind=8), dimension(:,:), intent(inout) :: x0
    integer (kind=4), intent(in) :: iout

    real (kind=8) :: xmin, xmax, ymin, ymax
    real (kind=8) :: amax
    integer (kind=4) :: unit
    integer (kind=4) :: i
    
    xmin = minval(x0(:,1))
    xmax = maxval(x0(:,1))
    ymin = minval(x0(:,2))
    ymax = maxval(x0(:,2))

    
    amax = max(-xmin, xmax)
    amax = max( amax, -ymin)
    amax = max( amax, ymax)

    unit = 9
    open(unit, file="gscript")
    write(unit,'(a11,f15.6,a1,f15.6,a1)') 'set xrange[,', -amax, ':', amax,']'
    write(unit,'(a11,f15.6,a1,f15.6,a1)') 'set yrange[,', -amax, ':', amax,']'
!    write(unit,'(a16)') 'set terminal png'

    do i = 1, iout
!      write(unit,'(a,i3.3,a)') "set output'a.",i,".png'"
      write(unit,'(a9,i3.3,a)') "plot 'a.",i,"' using 1:2"
    enddo
    close(unit)


    return
  end subroutine

!------------------------------------------------------
!
!
!
  subroutine print_profile(galaxy, rin, rout, rscale, mass, eps, &
       theta, phi, opt, heat, n)

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
    integer (kind=4), intent(in) :: galaxy
    real (kind=8), intent(in) :: rin
    real (kind=8), intent(in) :: rout
    real (kind=8), dimension(3), intent(in) :: rscale
    real (kind=8), intent(in) :: mass 
    real (kind=8), intent(in) :: eps
    real (kind=8), intent(in) :: theta
    real (kind=8), intent(in) :: phi
    integer(kind=4), intent(in) :: opt
    real (kind=8), intent(in) :: heat
    integer (kind=4) :: n

    write(*,*)
    write(*,'(a)') '----------------------------------'
    write(*,'(a,i3)') 'GALAXY =',galaxy
    write(*,'(a,f12.6)') 'mass        = ', mass
    write(*,'(a,f12.6)') 'epsilon     = ', eps
    write(*,'(a,f12.6)') 'rin         = ', rin 
    write(*,'(a,f12.6)') 'rout        = ', rout 
    write(*,'(a,f12.6)') 'rscale      = ', rscale(1)
    write(*,'(a,f12.6)') 'rscale      = ', rscale(2)
    write(*,'(a,f12.6)') 'rscale      = ', rscale(3)

    write(*,'(a,f12.6)') 'theta       = ', theta  
    write(*,'(a,f12.6)') 'phi         = ', phi  
    write(*,'(a,i12)')   'opt         = ', opt 
    write(*,'(a,f12.6)') 'heat        = ', heat  
    write(*,'(a,i12)')   'particles   = ', n  
    write(*,'(a)') '----------------------------------'

    return
  end subroutine print_profile


!------------------------------------------------------
!
!
!
  subroutine print_collision(n, time, inclination_degree, omega_degree, &
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
    integer (kind=4) :: n
    real (kind=8) :: time  
    real (kind=8) :: inclination_degree
    real (kind=8) :: omega_degree
    real (kind=8) :: rmin
    real (kind=8) :: velocity_factor

    real (kind=8) :: h
    integer (kind=4) :: nstep
    integer (kind=4) :: nout

    write(*,*)
    write(*,'(a)') '----------------------------------'
    write(*,'(a)') 'COLLISION PARAMETERS'
    write(*,'(a,i12)')    'n           = ', n
    write(*,'(a,f12.5)')  'time        = ', time
    write(*,'(a,f12.5)')  'inclination = ', inclination_degree
    write(*,'(a,f12.5)')  'omega       = ', omega_degree
    write(*,'(a,f12.5)')  'rmin        = ', rmin
    write(*,'(a,f12.5)')  'velocity    = ', velocity_factor
    write(*,'(a,f12.5)')  'h           = ', h
    write(*,'(a,i12)')    'nstep       = ', nstep
    write(*,'(a,i12)')    'nout        = ', nout 
    write(*,'(a)') '----------------------------------'
    write(*,*)

    return
  end subroutine print_collision

!
!-----------------------------------------------------------------
!

!----------------------------------------------------------------
!
  subroutine octave_parameters_out(mass1, theta1, phi1, rout1, mass2, &
       theta2, phi2, rout2, pos, vel, tend, x00, n, nunit)

    real (kind=8), intent(in) :: mass1, theta1, phi1, rout1
    real (kind=8), intent(in) :: mass2, theta2, phi2, rout2
    real (kind=8), dimension(3) , intent(in) :: pos, vel
    real (kind=8), intent(in) :: tend
    real (kind=8), dimension(6) , intent(in) :: x00
    integer (kind=4), intent(in) :: n, nunit
    
    write(nunit,'(a9,f12.6,a1)') '$mass1 = ',mass1,';'
    write(nunit,'(a9,f12.6,a1)') '$t1    = ',theta1,';'
    write(nunit,'(a9,f12.6,a1)') '$p1    = ', phi1,';'
    write(nunit,'(a9,f12.6,a1)') '$rout1 = ', rout1,';'
    
    write(nunit,'(a9,f12.6,a1)') '$mass2 = ',mass2,';'
    write(nunit,'(a9,f12.6,a1)') '$t2    = ', theta2,';'
    write(nunit,'(a9,f12.6,a1)') '$p2    = ', phi2,';'
    write(nunit,'(a9,f12.6,a1)') '$rout2 = ', rout2,';'
    
    write(nunit,'(a9,f12.6,a1)') '$xf    = ', pos(1),';'
    write(nunit,'(a9,f12.6,a1)') '$yf    = ', pos(2),';'
    write(nunit,'(a9,f12.6,a1)') '$zf    = ', pos(3),';'
    write(nunit,'(a9,f12.6,a1)') '$vxf   = ', vel(1),';'
    write(nunit,'(a9,f12.6,a1)') '$vyf   = ', vel(2),';'
    write(nunit,'(a9,f12.6,a1)') '$vzf   = ', vel(3),';'
    
    write(nunit,'(a9,f12.6,a1)') '$x     = ', x00(1),';'
    write(nunit,'(a9,f12.6,a1)') '$y     = ', x00(2),';'
    write(nunit,'(a9,f12.6,a1)') '$z     = ', x00(3),';'
    write(nunit,'(a9,f12.6,a1)') '$vx    = ', x00(4),';'
    write(nunit,'(a9,f12.6,a1)') '$vy    = ', x00(5),';'
    write(nunit,'(a9,f12.6,a1)') '$vz    = ', x00(6),';'
    write(nunit,'(a9,f12.6,a1)') '$t     = ', tend,';'
    
    return
  end subroutine octave_parameters_out


!
!-----------------------------------------------------------------
!
  subroutine read_parameter_file(unit)
    implicit none
    integer (kind=4), intent(in) :: unit
    character *40 label
    character *80 buffer 
    real (kind=8) :: val 
    integer (kind=4) :: stat
    
    READ(unit,*,IOSTAT=stat) buffer
    DO WHILE(stat .eq. 0) 
        call split_str(buffer,label,val)

!        print*,'buffer = ',buffer
!        print*,'label = ',label
!        print*,'val = ',val

!       determine parameter value and set it
        if(label .eq. 'potential_type') then
            potential_type = val
        else if(label .eq. 'mass1') then
            mass1 = val
        else if(label .eq. 'mass2') then
            mass2 = val
        else if(label .eq. 'epsilon1') then
            epsilon1 = val
        else if(label .eq. 'epsilon2') then
            epsilon2 = val
        else if(label .eq. 'rin1') then
            rin1 = val
        else if(label .eq. 'rin2') then
            rin2 = val
        else if(label .eq. 'rout1') then
            rout1 = val
        else if(label .eq. 'rout2') then
            rout2 = val
        else if(label .eq. 'theta1') then
            theta1 = val
        else if(label .eq. 'theta2') then
            theta2 = val
        else if(label .eq. 'phi1') then
            phi1 = val
        else if(label .eq. 'phi2') then
            phi2 = val
        else if(label .eq. 'opt1') then
            opt1 = val
        else if(label .eq. 'opt2') then
            opt2 = val
        else if(label .eq. 'heat1') then
            heat1 = val
        else if(label .eq. 'heat2') then
            heat2 = val
        else if(label .eq. 'n1') then
            n1 = val
        else if(label .eq. 'n2') then
            n2 = val
        else if(label .eq. 'inclination_degree') then
            inclination_degree = val
        else if(label .eq. 'omega_degree') then
            omega_degree = val
        else if(label .eq. 'rmin') then
            rmin = val
        else if(label .eq. 'velocity_factor') then
            velocity_factor = val
        else if(label .eq. 'tstart') then
            time = val
            tstart = val
            tIsSet = .true.
        else if(label .eq. 'tend') then
            tend = val
        else if(label .eq. 'h') then
            h = val
        else if(label .eq. 'rx') then
            use_sec_vec = .true. 
            sec_vec(1) = val
        else if(label .eq. 'ry') then
            use_sec_vec = .true. 
            sec_vec(2) = val
        else if(label .eq. 'rz') then
            use_sec_vec = .true. 
            sec_vec(3) = val
        else if(label .eq. 'vx') then
            use_sec_vec = .true. 
            sec_vec(4) = val
        else if(label .eq. 'vy') then
            use_sec_vec = .true. 
            sec_vec(5) = val
        else if(label .eq. 'vz') then
            use_sec_vec = .true. 
            sec_vec(6) = val
        else if(label .eq. 'rscale1') then
            rscale1 = val
        else if(label .eq. 'rscale2') then
            rscale2 = val
        else if(label .eq. 'rscale11') then
            rscale1(1) = val
        else if(label .eq. 'rscale12') then
            rscale1(2) = val
        else if(label .eq. 'rscale13') then
            rscale1(3) = val
        else if(label .eq. 'rscale21') then
            rscale2(1) = val
        else if(label .eq. 'rscale22') then
            rscale2(2) = val
        else if(label .eq. 'rscale23') then
            rscale2(3) = val
        else
            print*,'skipping line ',buffer
        end if

        READ(unit,*,IOSTAT=stat) buffer
    END DO
  end subroutine read_parameter_file


!----------------------------------------------------------------
!
  subroutine split_str(str,label,val)
    character *40, intent(out) :: label
    character *80, intent(in) :: str 
    real (kind=8), intent(out) :: val 
    character strt
    integer (kind=4) :: ind, slen

    strt = str(1:1)

    if(strt .eq. '!' .OR. strt .eq. '#' .OR. strt .eq. '/') then
        label = '!'
        val = 0
        return
    end if

    ind = INDEX(str,'=',.false.)

    if(ind .eq. 0) then
        label = '!'
        val = 0
    end if

    label=STR(1:ind-1)

    slen = LEN_TRIM(str)
    read(str(ind+1:slen),*)val

  end subroutine split_str

!----------------------------------------------------------------
!
  subroutine write_parameter_file(unit)
    implicit none
    integer (kind=4), intent(in) :: unit

    write(unit,*),'potential_type=',potential_type
    write(unit,*),'mass1=',mass1
    write(unit,*),'mass2=',mass2
    write(unit,*),'epsilon1=',epsilon1
    write(unit,*),'epsilon2=',epsilon2
    write(unit,*),'rin1=',rin1
    write(unit,*),'rin2=',rin2
    write(unit,*),'rout1=',rout1
    write(unit,*),'rout2=',rout2
    write(unit,*),'theta1=',theta1
    write(unit,*),'theta2=',theta2
    write(unit,*),'phi1=',phi1
    write(unit,*),'phi2=',phi2
    write(unit,*),'opt1=',opt1
    write(unit,*),'opt2=',opt2
    write(unit,*),'heat1=',heat1
    write(unit,*),'heat2=',heat2
    write(unit,*),'n1=',n1
    write(unit,*),'n2=',n2
    write(unit,*),'inclination_degree=',inclination_degree
    write(unit,*),'omega_degree=',omega_degree
    write(unit,*),'rmin=',rmin
    write(unit,*),'velocity_factor=',velocity_factor
    write(unit,*),'tstart=',tstart
    write(unit,*),'tend=',tend
    write(unit,*),'h=',h
    write(unit,*),'rx=',sec_vec(1)
    write(unit,*),'ry=',sec_vec(2)
    write(unit,*),'rz=',sec_vec(3)
    write(unit,*),'vx=',sec_vec(4)
    write(unit,*),'vy=',sec_vec(5)
    write(unit,*),'vz=',sec_vec(6)
    write(unit,*),'rscale11=',rscale1(1)
    write(unit,*),'rscale12=',rscale1(2)
    write(unit,*),'rscale13=',rscale1(3)
    write(unit,*),'rscale21=',rscale2(1)
    write(unit,*),'rscale22=',rscale2(2)
    write(unit,*),'rscale23=',rscale2(3)
  end subroutine write_parameter_file

!----------------------------------------------------------------
!
  subroutine parse_state_info_string(buffer)
    implicit none
    character *300, intent(in) ::  buffer
    real (kind=8), dimension(22) :: infos
    integer :: ind, ind2, len, i

    len = LEN_TRIM(buffer)
    ind = 0
    ind2 = 0
    i=1
    infos = 0

! There are 22 parameters to parse
    do while(i.lt.23)
        ind = ind2+1
        ind2 = INDEX(buffer(ind:len), ',', .false.)

        if(ind2 .gt. 0) then
            ind2 = ind2 + ind-2
            read(buffer(ind:ind2),*)infos(i)
            ind2 = ind2 + 1
        else
            exit
        end if
        i= i+1
    end do    

! set the parameters
            potential_type = 0
    sec_vec(1) = infos(1)
    sec_vec(2) = infos(2)
    sec_vec(3) = infos(3)
    sec_vec(4) = infos(4)
    sec_vec(5) = infos(5)
    sec_vec(6) = infos(6)
    mass1 = infos(7)
    mass2 = infos(8)
    rout1 = infos(9)
    rout2 = infos(10)
    phi1 = infos(11)
    phi2 = infos(12)
    theta1 = infos(13)
    theta2 = infos(14)
    epsilon1 = infos(15)
    epsilon2 = infos(16)
    rscale1(1) = infos(17)
    rscale1(2) = infos(18)
    rscale1(3) = infos(19)
    rscale2(1) = infos(20)
    rscale2(2) = infos(21)
    rscale2(3) = infos(22)
    use_sec_vec = .true. 

  end subroutine parse_state_info_string
  

!----------------------------------------------------------------
!
!----------------------------------------------------------------
!
!----------------------------------------------------------------
!
!----------------------------------------------------------------
!
!----------------------------------------------------------------
!
end module io_module
