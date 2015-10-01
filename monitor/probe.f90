module probe
  use precision
  use geometry
  use partition
  implicit none
  
  ! Number of probes and variables to probe
  integer :: nprobefiles
  integer :: nprobes
  integer :: nvars
  character(len=str_medium), dimension(:), allocatable :: varnames
  character(len=str_medium), dimension(:), pointer :: probe_filename
  integer, dimension(:), pointer :: lengths ! number of probes in each file

  ! Probe locations
  real(WP), dimension(:), pointer :: probe_x,probe_y,probe_z
  logical,  dimension(:), pointer :: probe_local
  
  ! Interpolation
  real(WP), dimension(:), pointer :: interp_x,interp_xm
  real(WP), dimension(:), pointer :: interp_y,interp_ym
  real(WP), dimension(:), pointer :: interp_z,interp_zm
  integer,  dimension(:), pointer :: index_x,index_xm
  integer,  dimension(:), pointer :: index_y,index_ym
  integer,  dimension(:), pointer :: index_z,index_zm
  
  ! Values to monitor
  real(WP), dimension(:), pointer :: mval



end module probe


! ===================== !
! Initialize the probes !
! ===================== !
subroutine probe_init
  use probe
  use parser
  use data
  implicit none
  integer :: nargs,n,i,j,k,isc,m,ipfile
  logical :: isdef, istrue
  real(WP), dimension(:,:), pointer :: list
  integer ::  max_probe_files = 10 
  real(WP), dimension(:), pointer :: xvals 

  ! Any probes ?
  call parser_is_defined('Probe locations',isdef)
  if (isdef) then
     call parser_getsize('Probe locations',nargs)
     if (mod(nargs,3).ne.0) then
        call die('probe_init: incorrect number of coordinates for probes locations')
     else
        nprobes = nargs/3
     end if
  else
     nprobes = 0
     return
  end if

  ! Additional variable to probe
  call parser_is_defined('Other probe variables',isdef)
  if (isdef) then
     call parser_getsize('Other probe variables',nvars)
     allocate(varnames(nvars))
     call parser_read('Other probe variables',varnames)
  else
     nvars = 0
  end if

  ! Allocate arrays
  allocate(list(3,nprobes))
  allocate(probe_x(nprobes))
  allocate(probe_y(nprobes))
  allocate(probe_z(nprobes))
  allocate(probe_local(nprobes))
  allocate(interp_x (nprobes))
  allocate(interp_xm(nprobes))
  allocate(interp_y (nprobes))
  allocate(interp_ym(nprobes))
  allocate(interp_z (nprobes))
  allocate(interp_zm(nprobes))
  allocate(index_x (nprobes))
  allocate(index_xm(nprobes))
  allocate(index_y (nprobes))
  allocate(index_ym(nprobes))
  allocate(index_z (nprobes))
  allocate(index_zm(nprobes))
  allocate(lengths(max_probe_files))
  allocate(xvals(max_probe_files))

  ! Read the locations
  call parser_read('Probe locations',list)
  probe_x = list(1,:)
  probe_y = list(2,:)
  probe_z = list(3,:)
   
  ! Separate probe files ? 
  nprobefiles = 1
  lengths(1) = nprobes
  xvals(1) = 0
  call parser_is_defined('Separate probe files',isdef)
  if (isdef) then
     call parser_read('Separate probe files',istrue)
     if (istrue) then
        lengths(1) = 0
        xvals(1) = list(1,1)
        do n=1,nprobes
           if (list(1,n).ne.xvals(nprobefiles)) then
              nprobefiles = nprobefiles + 1
              if (nprobefiles.gt.max_probe_files) then
                 call die('probe_init: too many probe files (x-locations)')
              end if
              xvals(nprobefiles) = list(1,n)
              lengths(nprobefiles) = 1
           else
              lengths(nprobefiles) = lengths(nprobefiles) + 1
           end if 
        end do
     end if
  end if

  ! Get interpolation
  do n=1,nprobes
     if (probe_x(n).lt.x(imin) .or. probe_x(n).ge.x(imax+1)) &
          call die('probe_init: probe x location not in domain')
     if (probe_y(n).lt.y(jmin) .or. probe_y(n).ge.y(jmax+1)) &
          call die('probe_init: probe y location not in domain')
     if (probe_z(n).lt.z(kmin) .or. probe_z(n).ge.z(kmax+1)) &
          call die('probe_init: probe z location not in domain')
     
     if ( probe_x(n).lt.x(imin_) .or. probe_x(n).ge.x(imax_+1) .or. &
          probe_y(n).lt.y(jmin_) .or. probe_y(n).ge.y(jmax_+1) .or. &
          probe_z(n).lt.z(kmin_) .or. probe_z(n).ge.z(kmax_+1)) then
        probe_local(n) = .false.
     else
        probe_local(n) = .true.

        ! Interpolation in x
        i = imin_
        do while(x(i+1).le.probe_x(n)) 
           i = i+1
        end do
        index_x(n)  = i
        interp_x(n) = (probe_x(n)-x(i))/(x(i+1)-x(i))
        i = imino_
        do while(xm(i+1).le.probe_x(n)) 
           i = i+1
        end do
        index_xm(n)  = i
        interp_xm(n) = (probe_x(n)-xm(i))/(xm(i+1)-xm(i))
        
        ! Interpolation in y
        j = jmin_
        do while(y(j+1).le.probe_y(n)) 
           j = j+1
        end do
        index_y(n)  = j
        interp_y(n) = (probe_y(n)-y(j))/(y(j+1)-y(j))
        j = jmino_
        do while(ym(j+1).le.probe_y(n)) 
           j = j+1
        end do
        index_ym(n)  = j
        interp_ym(n) = (probe_y(n)-ym(j))/(ym(j+1)-ym(j))
        
        ! Interpolation in z
        k = kmin_
        do while(z(k+1).le.probe_z(n)) 
           k = k+1
        end do
        index_z(n)  = k
        interp_z(n) = (probe_z(n)-z(k))/(z(k+1)-z(k))
        k = kmino_
        do while(zm(k+1).le.probe_z(n)) 
           k = k+1
        end do
        index_zm(n)  = k
        interp_zm(n) = (probe_z(n)-zm(k))/(zm(k+1)-zm(k))
     end if
  end do
  
  ! Create files to monitor
  allocate(probe_filename(nprobefiles))
  do ipfile=1,nprobefiles
     write(probe_filename(ipfile),'(a9,e12.6e2)') 'probes-x-',xvals(ipfile)
     if ( (ipfile.eq.1) .and. (.not.(istrue.and.isdef)) ) probe_filename(1)='probes'
     call monitor_create_file_step(probe_filename(ipfile),(4+nscalar+nvars)*lengths(ipfile))
     do n=1,lengths(ipfile)
        call monitor_set_header(1+(n-1)*(4+nscalar+nvars),'U','d')
        call monitor_set_header(2+(n-1)*(4+nscalar+nvars),'V','d')
        call monitor_set_header(3+(n-1)*(4+nscalar+nvars),'W','d')
        call monitor_set_header(4+(n-1)*(4+nscalar+nvars),'P','d')
        do isc=1,nscalar
           call monitor_set_header(4+isc+(n-1)*(4+nscalar+nvars),trim(SC_name(isc)),'d')
        end do
        do m=1,nvars
           call monitor_set_header(4+nscalar+m+(n-1)*(4+nscalar+nvars),trim(varnames(m)),'d')
        end do
     end do
  end do
  allocate(mval((4+nscalar+nvars)*nprobes))

  return
end subroutine probe_init


! ================== !
! Monitor the probes !
! ================== !
subroutine probe_monitor
  use probe
  use data
  implicit none
  integer  :: n,isc,i,j,k,m,ipfile,nprev
  real(WP) :: wx1,wy1,wz1,wx2,wy2,wz2,tmp
  real(WP), dimension(:,:,:,:), allocatable :: varsdata 
  
  ! Nothing to do if no probes
  if (nprobes.eq.0) return

  ! Get chemtable data for relevant vairables 
  if (nvars.ne.0) then
     allocate(varsdata(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nvars))
     do n=1,nvars
        call chemtable_lookup(trim(varnames(n)), varsdata(:,:,:,n))
     end do
  end if

  ! Compute local values at the probes
  mval = 0.0_WP
  do n=1,nprobes
     if (probe_local(n)) then
        ! Velocity - U
        i = index_x (n); wx2 = interp_x (n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        mval(1+(n-1)*(4+nscalar+nvars)) = &
             + wx1 * ( wy1 * (wz1*U(i  ,j  ,k)+wz2*U(i  ,j  ,k+1)) &
                     + wy2 * (wz1*U(i  ,j+1,k)+wz2*U(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*U(i+1,j  ,k)+wz2*U(i+1,j  ,k+1)) &
                     + wy2 * (wz1*U(i+1,j+1,k)+wz2*U(i+1,j+1,k+1)) )
        ! Velocity - V
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_y (n); wy2 = interp_y (n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        mval(2+(n-1)*(4+nscalar+nvars)) = &
             + wx1 * ( wy1 * (wz1*V(i  ,j  ,k)+wz2*V(i  ,j  ,k+1)) &
                     + wy2 * (wz1*V(i  ,j+1,k)+wz2*V(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*V(i+1,j  ,k)+wz2*V(i+1,j  ,k+1)) &
                     + wy2 * (wz1*V(i+1,j+1,k)+wz2*V(i+1,j+1,k+1)) )
        ! Velocity - W
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_z (n); wz2 = interp_z (n); wz1 = 1.0_WP-wz2
        mval(3+(n-1)*(4+nscalar+nvars)) = &
             + wx1 * ( wy1 * (wz1*W(i  ,j  ,k)+wz2*W(i  ,j  ,k+1)) &
                     + wy2 * (wz1*W(i  ,j+1,k)+wz2*W(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*W(i+1,j  ,k)+wz2*W(i+1,j  ,k+1)) &
                     + wy2 * (wz1*W(i+1,j+1,k)+wz2*W(i+1,j+1,k+1)) )
        ! Pressure
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        mval(4+(n-1)*(4+nscalar+nvars)) = &
             + wx1 * ( wy1 * (wz1*P(i  ,j  ,k)+wz2*P(i  ,j  ,k+1)) &
                     + wy2 * (wz1*P(i  ,j+1,k)+wz2*P(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*P(i+1,j  ,k)+wz2*P(i+1,j  ,k+1)) &
                     + wy2 * (wz1*P(i+1,j+1,k)+wz2*P(i+1,j+1,k+1)) )
        ! Scalars
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        do isc=1,nscalar
           mval(4+isc+(n-1)*(4+nscalar+nvars)) = &
                + wx1 * ( wy1 * (wz1*SC(i  ,j  ,k,isc)+wz2*SC(i  ,j  ,k+1,isc)) &
                        + wy2 * (wz1*SC(i  ,j+1,k,isc)+wz2*SC(i  ,j+1,k+1,isc)) )&
                + wx2 * ( wy1 * (wz1*SC(i+1,j  ,k,isc)+wz2*SC(i+1,j  ,k+1,isc)) &
                        + wy2 * (wz1*SC(i+1,j+1,k,isc)+wz2*SC(i+1,j+1,k+1,isc)) )
        end do
        ! Other Variables
        do m=1,nvars
           mval(4+nscalar+m+(n-1)*(4+nscalar+nvars)) = &
                + wx1 * ( wy1 * (wz1*varsdata(i  ,j  ,k,m)+wz2*varsdata(i  ,j  ,k+1,m)) &
                        + wy2 * (wz1*varsdata(i  ,j+1,k,m)+wz2*varsdata(i  ,j+1,k+1,m)) )&
                + wx2 * ( wy1 * (wz1*varsdata(i+1,j  ,k,m)+wz2*varsdata(i+1,j  ,k+1,m)) &
                        + wy2 * (wz1*varsdata(i+1,j+1,k,m)+wz2*varsdata(i+1,j+1,k+1,m)) )
        end do

     end if
  end do
  
  ! Get the global value
  do n=1,(4+nscalar+nvars)*nprobes
     call parallel_sum(mval(n),tmp)
     mval(n) = tmp
  end do

  ! Transfer to monitor
  nprev = 0
  do ipfile = 1,nprobefiles
     call monitor_select_file(probe_filename(ipfile))
     call monitor_set_array_values( mval(nprev+1 : nprev+(4+nvars+nscalar)*lengths(ipfile) ) )
     nprev = nprev + (4 + nvars + nscalar)*lengths(ipfile)
  end do

  return
end subroutine probe_monitor

