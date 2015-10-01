module implicit
  use precision
  use string
  use velocity
  implicit none
  
  ! Use of implicit solver
  logical :: implicit_x
  logical :: implicit_y
  logical :: implicit_z
  logical :: implicit_any
  logical :: implicit_auto
  
  ! Specific use of implicit solver
  logical :: imp_conv_x,imp_conv_y,imp_conv_z
  logical :: imp_visc_x,imp_visc_y,imp_visc_z
  integer :: ndx,ndy,ndz
  
  ! Work variables
  real(WP), dimension(:), pointer :: Rcyl
  real(WP), dimension(:), pointer :: Rcyl2
  
  real(WP) :: inverse_burning_timescale
  
  ! Values to monitor
  real(WP) :: CFLc_x,CFLc_y,CFLc_z,CFLc_axis
  real(WP) :: CFLv_x,CFLv_y,CFLv_z,CFLv_axis
  real(WP) :: CFLd_x,CFLd_y,CFLd_z,CFLd_axis
  real(WP) :: CFLa_x,CFLa_y,CFLa_z,CFLa_axis
  
contains
  
  ! Check which terms should be implicit
  subroutine implicit_check
    use parallel
    implicit none
    
    if (.not.implicit_auto) return
    
    ! Check if implicit is possible -> at best:
    ! - pentadiagonal in parallel
    ! - polydiagonal in serial
    imp_conv_x = .false.;imp_visc_x = .false.
    if (implicit_x) then
       if ((vel_conv_order.eq.2 .or. npx.eq.1) .and. (CFLc_x.ge.0.5_WP)) imp_conv_x = .true.
       if ((vel_visc_order.eq.2 .or. npx.eq.1) .and. (CFLv_x.ge.0.5_WP)) imp_visc_x = .true.
    end if
    imp_conv_y = .false.;imp_visc_y = .false.
    if (implicit_y) then
       if ((vel_conv_order.eq.2 .or. npy.eq.1) .and. (CFLc_y.ge.0.5_WP)) imp_conv_y = .true.
       if ((vel_visc_order.eq.2 .or. npy.eq.1) .and. (CFLv_y.ge.0.5_WP)) imp_visc_y = .true.
    end if
    imp_conv_z = .false.;imp_visc_z = .false.
    if (implicit_z) then
       if ((vel_conv_order.eq.2 .or. npz.eq.1) .and. (CFLc_z.ge.0.5_WP)) imp_conv_z = .true.
       if ((vel_visc_order.eq.2 .or. npz.eq.1) .and. (CFLv_z.ge.0.5_WP)) imp_visc_z = .true.
    end if
    
    ! Obtain number of diagonals - scheme dependent
    if (imp_conv_x) ndx = 2*vel_conv_order-1
    if (imp_visc_x) ndx = max(ndx,2*vel_visc_order-1)
    if (imp_conv_y) ndy = 2*vel_conv_order-1
    if (imp_visc_y) ndy = max(ndy,2*vel_visc_order-1)
    if (imp_conv_z) ndz = 2*vel_conv_order-1
    if (imp_visc_z) ndz = max(ndz,2*vel_visc_order-1)
    
    return
  end subroutine implicit_check
  
end module implicit


! ============================== !
! Initialize the implicit module !
! ============================== !
subroutine implicit_init
  use implicit
  use parser
  use parallel
  use data
  implicit none
  
  character(len=str_medium) :: implicit_dir
  
  implicit_x = .false.
  implicit_y = .false.
  implicit_z = .false.
  
  ! Read in the implicit directions
  call parser_read('Implicit directions',implicit_dir)
  if (index(implicit_dir,'x').ne.0) implicit_x = .true.
  if (index(implicit_dir,'y').ne.0) implicit_y = .true.
  if (index(implicit_dir,'z').ne.0) implicit_z = .true.
  if (trim(implicit_dir).eq.'auto') then
     implicit_x = .true.
     implicit_y = .true.
     implicit_z = .true.
     implicit_auto = .true.
  end if
  
  ! Test if everything has been correctly detected
  implicit_any = implicit_x .or. implicit_y .or. implicit_z
  if (.not.(implicit_any .or. trim(implicit_dir).eq.'none')) call die('implicit_init: unknown directions')
  
  ! Account for 2D and 1D problems
  if (nx.eq.1) implicit_x = .false.
  if (ny.eq.1) implicit_y = .false.
  if (nz.eq.1) implicit_z = .false.
  implicit_any = implicit_x .or. implicit_y .or. implicit_z
  
  ! Check if implicit is possible -> at best:
  ! - pentadiagonal in parallel
  ! - polydiagonal in serial
  imp_conv_x = .false.;imp_visc_x = .false.
  if (implicit_x) then
     if (vel_conv_order.eq.2 .or. npx.eq.1) imp_conv_x = .true.
     if (vel_visc_order.eq.2 .or. npx.eq.1) imp_visc_x = .true.
  end if
  imp_conv_y = .false.;imp_visc_y = .false.
  if (implicit_y) then
     if (vel_conv_order.eq.2 .or. npy.eq.1) imp_conv_y = .true.
     if (vel_visc_order.eq.2 .or. npy.eq.1) imp_visc_y = .true.
  end if
  imp_conv_z = .false.;imp_visc_z = .false.
  if (implicit_z) then
     if (vel_conv_order.eq.2 .or. npz.eq.1) imp_conv_z = .true.
     if (vel_visc_order.eq.2 .or. npz.eq.1) imp_visc_z = .true.
  end if
  
  ! Obtain number of diagonals - scheme dependent
  if (imp_conv_x) ndx = 2*vel_conv_order-1
  if (imp_visc_x) ndx = max(ndx,2*vel_visc_order-1)
  if (imp_conv_y) ndy = 2*vel_conv_order-1
  if (imp_visc_y) ndy = max(ndy,2*vel_visc_order-1)
  if (imp_conv_z) ndz = 2*vel_conv_order-1
  if (imp_visc_z) ndz = max(ndz,2*vel_visc_order-1)
  
  ! Implicit boundary conditions in cylindrical coordinates
  if (icyl.eq.1) allocate(Rcyl (imin_:imax_))
  if (icyl.eq.1) allocate(Rcyl2(imin_:imax_))
  
  ! Create new file to monitor at each iterations
  if (nscalar.eq.0) then
     if (.not.compressible) then
        call monitor_create_file_step('timestep',10)
     else
        call monitor_create_file_step('timestep',14)
     end if
  else
     if (.not.compressible) then
        call monitor_create_file_step('timestep',14)
     else
        call monitor_create_file_step('timestep',18)
     end if
  end if
  call monitor_set_header(1, 'dt','r')
  call monitor_set_header(2, 'CFL_max','r')
  call monitor_set_header(3, 'CFLc_x','r')
  call monitor_set_header(4, 'CFlc_y','r')
  call monitor_set_header(5, 'CFLc_z','r')
  call monitor_set_header(6, 'CFLc_axis','r')
  call monitor_set_header(7, 'CFLv_x','r')
  call monitor_set_header(8, 'CFLv_y','r')
  call monitor_set_header(9, 'CFLv_z','r')
  call monitor_set_header(10,'CFLv_axis','r')
  if (nscalar.ne.0) then
     call monitor_set_header(11, 'CFLd_x','r')
     call monitor_set_header(12, 'CFLd_y','r')
     call monitor_set_header(13, 'CFLd_z','r')
     call monitor_set_header(14, 'CFLd_axis','r')
  end if
  if (compressible) then
     if (nscalar.eq.0) then
        call monitor_set_header(11, 'CFLa_x','r')
        call monitor_set_header(12, 'CFLa_y','r')
        call monitor_set_header(13, 'CFLa_z','r')
        call monitor_set_header(14, 'CLFa_axis','r')
     else
        call monitor_set_header(15, 'CFLa_x','r')
        call monitor_set_header(16, 'CFLa_y','r')
        call monitor_set_header(17, 'CFLa_z','r')
        call monitor_set_header(18, 'CLFa_axis','r')
     end if
  end if
  
  return
end subroutine implicit_init


! ================================================ !
! Compute the maximum CFL from velocity components !
! ================================================ !
subroutine velocity_CFL
  use implicit
  use data
  use parallel
  use memory
  implicit none
  
  integer  :: i,j,k
  real(WP) :: max_CFLc_x,max_CFLc_y,max_CFLc_z,max_CFLc_axis
  real(WP) :: max_CFLv_x,max_CFLv_y,max_CFLv_z,max_CFLv_axis
  real(WP) :: max_CFLd_x,max_CFLd_y,max_CFLd_z,max_CFLd_axis
  real(WP) :: max_CFLa_x,max_CFLa_y,max_CFLa_z,max_CFLa_axis
  real(WP) :: max_CFL_lvlset
  real(WP) :: a
  
  ! Set the CFLs to zero
  max_CFLc_x = 0.0_WP
  max_CFLc_y = 0.0_WP
  max_CFLc_z = 0.0_WP
  max_CFLc_axis = 0.0_WP
  max_CFLv_x = 0.0_WP
  max_CFLv_y = 0.0_WP
  max_CFLv_z = 0.0_WP
  max_CFLv_axis = 0.0_WP
  max_CFLd_x = 0.0_WP
  max_CFLd_y = 0.0_WP
  max_CFLd_z = 0.0_WP
  max_CFLd_axis = 0.0_WP
  max_CFLa_x = 0.0_WP
  max_CFLa_y = 0.0_WP
  max_CFLa_z = 0.0_WP
  max_CFLa_axis = 0.0_WP
  
  ! Calculate the CFLs
  if (icyl.EQ.1) then
     ! Cylindrical
     !$OMP PARALLEL DO PRIVATE(CFLc_x,CFLc_y,CFLc_z,CFLv_x,CFLv_y,CFLv_z) &
     !$OMP REDUCTION(max:max_CFLc_x,max_CFLc_y,max_CFLc_z,max_CFLv_x,max_CFLv_y,max_CFLv_z)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              CFLc_x = abs( U(i,j,k) ) * dxmi(i-1)
              CFLc_y = abs( V(i,j,k) ) * dymi(j-1)
              CFLc_z = abs( W(i,j,k) ) * dzi*ymi(j)
              CFLv_x = abs( 4.0_WP*VISC(i,j,k) ) * dxi(i)**2       / RHOmid(i,j,k)
              CFLv_y = abs( 4.0_WP*VISC(i,j,k) ) * dyi(j)**2       / RHOmid(i,j,k)
              CFLv_z = abs( 4.0_WP*VISC(i,j,k) ) * (dzi*ymi(j))**2 / RHOmid(i,j,k)
              max_CFLc_x = max(max_CFLc_x,CFLc_x)
              max_CFLc_y = max(max_CFLc_y,CFLc_y)
              max_CFLc_z = max(max_CFLc_z,CFLc_z)
              max_CFLv_x = max(max_CFLv_x,CFLv_x)
              max_CFLv_y = max(max_CFLv_y,CFLv_y)
              max_CFLv_z = max(max_CFLv_z,CFLv_z)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     ! Centerline CFL
     if (jproc.eq.1) then
        !$OMP PARALLEL DO PRIVATE(j,CFLc_axis,CFLv_axis) &
        !$OMP REDUCTION(max:max_CFLc_axis,max_CFLv_axis)
        do k=kmin_,kmax_
           j = jmin_
           do i=imin_,imax_
              CFLc_axis = abs( V(i,j,k) ) * dymi(j-1)
              CFLv_axis = abs( 4.0_WP*VISC(i,j,k) ) * dyi(j)**2 / RHOmid(i,j,k)
              max_CFLc_axis = max(max_CFLc_axis,CFLc_axis)
              max_CFLv_axis = max(max_CFLv_axis,CFLv_axis)
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     ! For scalars
     if (nscalar.ne.0) then
        !$OMP PARALLEL DO PRIVATE(CFLd_x,CFLd_y,CFLd_z) &
        !$OMP REDUCTION(max:max_CFLd_x,max_CFLd_y,max_CFLd_z)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 CFLd_x = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * dxi(i)**2       / RHOmid(i,j,k)
                 CFLd_y = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * dyi(j)**2       / RHOmid(i,j,k)
                 CFLd_z = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * (dzi*ymi(j))**2 / RHOmid(i,j,k)
                 max_CFLd_x = max(max_CFLd_x,CFLd_x)
                 max_CFLd_y = max(max_CFLd_y,CFLd_y)
                 max_CFLd_z = max(max_CFLd_z,CFLd_z)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        if (jproc.eq.1) then
           !$OMP PARALLEL DO PRIVATE(j,CFLd_axis) REDUCTION(max:max_CFLd_axis)
           do k=kmin_,kmax_
              j = jmin_
              do i=imin_,imax_
                 CFLd_axis = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * dyi(j)**2 / RHOmid(i,j,k)
                 max_CFLd_axis = max(max_CFLd_axis,CFLd_axis)
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     end if
     ! For acoustics
     if (compressible) then
        call combustion_soundspeed(tmp1)
        !$OMP PARALLEL DO PRIVATE(CFLa_x,CFLa_y,CFLa_z,a) &
        !$OMP REDUCTION(max:max_CFLa_x,max_CFLa_y,max_CFLa_z)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 CFLa_x = tmp1(i,j,k) * dxmi(i-1)
                 CFLa_y = tmp1(i,j,k) * dymi(j-1)
                 CFLa_z = tmp1(i,j,k) * dzi*ymi(j)
                 max_CFLa_x = max(max_CFLa_x,CFLa_x)
                 max_CFLa_y = max(max_CFLa_y,CFLa_y)
                 max_CFLa_z = max(max_CFLa_z,CFLa_z)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        if (jproc.eq.1) then
           !$OMP PARALLEL DO PRIVATE(j,CFLa_axis,a) REDUCTION(max:max_CFLa_axis)
           do k=kmin_,kmax_
              j = jmin_
              do i=imin_,imax_
                 CFLa_axis = tmp1(i,j,k) * dymi(j-1)
                 max_CFLa_axis = max(max_CFLa_axis,CFLa_axis)
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     end if
  else
     ! Cartesian
     !$OMP PARALLEL DO PRIVATE(CFLc_x,CFLc_y,CFLc_z,CFLv_x,CFLv_y,CFLv_z) &
     !$OMP REDUCTION(max:max_CFLc_x,max_CFLc_y,max_CFLc_z,max_CFLv_x,max_CFLv_y,max_CFLv_z)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              CFLc_x = abs( U(i,j,k) ) * dxmi(i-1)
              CFLc_y = abs( V(i,j,k) ) * dymi(j-1)
              CFLc_z = abs( W(i,j,k) ) * dzi
              CFLv_x = abs( 4.0_WP*VISC(i,j,k) ) * dxi(i)**2 / RHOmid(i,j,k)
              CFLv_y = abs( 4.0_WP*VISC(i,j,k) ) * dyi(j)**2 / RHOmid(i,j,k)
              CFLv_z = abs( 4.0_WP*VISC(i,j,k) ) * dzi**2    / RHOmid(i,j,k)
              max_CFLc_x = max(max_CFLc_x,CFLc_x)
              max_CFLc_y = max(max_CFLc_y,CFLc_y)
              max_CFLc_z = max(max_CFLc_z,CFLc_z)
              max_CFLv_x = max(max_CFLv_x,CFLv_x)
              max_CFLv_y = max(max_CFLv_y,CFLv_y)
              max_CFLv_z = max(max_CFLv_z,CFLv_z)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     ! For scalars
     if (nscalar.ne.0) then
        !$OMP PARALLEL DO PRIVATE(CFLd_x,CFLd_y,CFLd_z) &
        !$OMP REDUCTION(max:max_CFLd_x,max_CFLd_y,max_CFLd_z)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 CFLd_x = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * dxi(i)**2 / RHOmid(i,j,k)
                 CFLd_y = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * dyi(j)**2 / RHOmid(i,j,k)
                 CFLd_z = abs( 4.0_WP*maxval(DIFF(i,j,k,:)) ) * dzi**2    / RHOmid(i,j,k)
                 max_CFLd_x = max(max_CFLd_x,CFLd_x)
                 max_CFLd_y = max(max_CFLd_y,CFLd_y)
                 max_CFLd_z = max(max_CFLd_z,CFLd_z)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     ! For acoustics
     if (compressible) then
        call combustion_soundspeed(tmp1)
        !$OMP PARALLEL DO PRIVATE(CFLa_x,CFLa_y,CFLa_z,a) &
        !$OMP REDUCTION(max:max_CFLa_x,max_CFLa_y,max_CFLa_z)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 CFLa_x = tmp1(i,j,k) * dxmi(i-1)
                 CFLa_y = tmp1(i,j,k) * dymi(j-1)
                 CFLa_z = tmp1(i,j,k) * dzi
                 max_CFLa_x = max(max_CFLa_x,CFLa_x)
                 max_CFLa_y = max(max_CFLa_y,CFLa_y)
                 max_CFLa_z = max(max_CFLa_z,CFLa_z)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  end if
  
  ! Get the global maxs
  max_CFLc_x = max_CFLc_x * dt_uvw
  max_CFLc_y = max_CFLc_y * dt_uvw
  max_CFLc_z = max_CFLc_z * dt_uvw
  max_CFLc_axis = max_CFLc_axis * dt_uvw
  max_CFLv_x = max_CFLv_x * dt_uvw
  max_CFLv_y = max_CFLv_y * dt_uvw
  max_CFLv_z = max_CFLv_z * dt_uvw
  max_CFLv_axis = max_CFLv_axis * dt_uvw
  max_CFLd_x = max_CFLd_x * dt_uvw
  max_CFLd_y = max_CFLd_y * dt_uvw
  max_CFLd_z = max_CFLd_z * dt_uvw
  max_CFLd_axis = max_CFLd_axis * dt_uvw
  max_CFLa_x = max_CFLa_x * dt_uvw
  max_CFLa_y = max_CFLa_y * dt_uvw
  max_CFLa_z = max_CFLa_z * dt_uvw
  max_CFLa_axis = max_CFLa_axis * dt_uvw
  call parallel_max(max_CFLc_x,CFLc_x)
  call parallel_max(max_CFLc_y,CFLc_y)
  call parallel_max(max_CFLc_z,CFLc_z)
  call parallel_max(max_CFLc_axis,CFLc_axis)
  call parallel_max(max_CFLv_x,CFLv_x)
  call parallel_max(max_CFLv_y,CFLv_y)
  call parallel_max(max_CFLv_z,CFLv_z)
  call parallel_max(max_CFLv_axis,CFLv_axis)
  call parallel_max(max_CFLd_x,CFLd_x)
  call parallel_max(max_CFLd_y,CFLd_y)
  call parallel_max(max_CFLd_z,CFLd_z)
  call parallel_max(max_CFLd_axis,CFLd_axis)
  call parallel_max(max_CFLa_x,CFLa_x)
  call parallel_max(max_CFLa_y,CFLa_y)
  call parallel_max(max_CFLa_z,CFLa_z)
  call parallel_max(max_CFLa_axis,CFLa_axis)
  
  ! Set the main CFL based on the implicit directions
  CFL = max(CFLc_x,CFLc_y)
  if (.not.imp_visc_x .and. nx.gt.1) CFL = max(CFL,CFLv_x,CFLd_x)
  
  if (.not.imp_conv_y .and. ny.gt.1) CFL = max(CFL,CFLc_y)
  if (.not.imp_visc_y .and. ny.gt.1) CFL = max(CFL,CFLv_y,CFLd_y)
  
  if (.not.imp_conv_z .and. nz.gt.1) CFL = max(CFL,CFLc_z)
  if (.not.imp_visc_z .and. nz.gt.1) CFL = max(CFL,CFLv_z,CFLd_z)
  
  if (icyl.eq.1 .and. nz.gt.1 .and. isect.eq.0) then
     CFL = max(CFL,CFLc_axis)
     CFL = max(CFL,CFLv_axis)
     CFL = max(CFL,CFLd_axis)
  end if
  if (icyl.eq.1 .and. nz.gt.1) then
     CFL = max(CFL,0.1_WP*CFLc_z)
  end if

  ! Set the main CFL with acoustic constraint: Always explicit
  if (compressible) then
     if (nx.gt.1) CFL = max(CFL,CFLa_x)
     if (ny.gt.1) CFL = max(CFL,CFLa_y)
     if (nz.gt.1) CFL = max(CFL,CFLa_z)
     if (icyl.eq.1 .and. nz.gt.1 .and. isect.eq.0) then
        CFL = max(CFL,CFLa_axis)
     end if
  end if
  
  ! Consider the level set source term
  if (use_lvlset) then
     max_CFL_lvlset = inverse_burning_timescale * dt_uvw
     CFL = max(CFL,max_CFL_lvlset)
  else
     max_CFL_lvlset = 0.0_WP
  end if
  
  ! Check if we can do better for the implicit solver
  call implicit_check
  
  return
end subroutine velocity_CFL


! ================================== !
! Treatment of the centerline CFL    !
! Modification of the modeled VISC   !
! and DIFF such that viscous effects !
! can never be cfl-limiting          !
! ================================== !
subroutine velocity_CFL_centerline
  use implicit
  use data
  use parallel
  implicit none
  
  integer :: i,j,k,isc
  real(WP) :: pred_CFLv_axis
  real(WP), parameter :: CFLv_max = 0.5_WP
  
  ! Only if cylindrical and not DNS
  if (icyl.eq.0 .or. .not.use_sgs) return
  
  ! For each centerline point
  if (jproc.eq.1) then
     !$OMP PARALLEL DO PRIVATE(pred_CFLv_axis)
     do k=kmin_,kmax_
        j=jmin_
        do i=imin_,imax_
           ! Predict diffusive centerline CFL
           pred_CFLv_axis = dt_uvw * abs( 4.0_WP*VISC(i,j,k) ) * dyi(j)**2 / RHOmid(i,j,k)
           ! Adjust VISC
           if (pred_CFLv_axis.gt.CFLv_max) then
              VISC(i,j,k) = VISC(i,j,k)*CFLv_max/pred_CFLv_axis
           end if
           ! Same for scalars
           do isc=1,nscalar
              ! Predict diffusive centerline CFL
              pred_CFLv_axis = dt_uvw * abs( 4.0_WP*DIFF(i,j,k,isc) ) * dyi(j)**2 / RHOmid(i,j,k)
              ! Adjust DIFF
              if (pred_CFLv_axis.gt.CFLv_max) then
                 DIFF(i,j,k,isc) = DIFF(i,j,k,isc)*CFLv_max/pred_CFLv_axis
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Update borders
  call boundary_update_border(VISC,'+','ym')
  do isc=1,nscalar
     call boundary_update_border(DIFF(:,:,:,isc),'+','ym')
  end do
  
  return
end subroutine velocity_CFL_centerline


! ================================== !
! Inverse the Residual of u          !
! By using approximate factorization !
! ================================== !
subroutine velocity_inverse_u
  use implicit
  use data
  use metric_generic
  use metric_velocity_conv
  use metric_velocity_visc
  use memory
  use time_info
  implicit none
  
  integer  :: i,j,k,st,n
  real(WP) :: RHOi,VISCi
  real(WP) :: dt2
  
  dt2 = dt_uvw/2.0_WP
  
  ! If purely explicit return
  if (.not.(imp_conv_x .or. imp_visc_x .or. &
            imp_conv_y .or. imp_visc_y .or. &
            imp_conv_z .or. imp_visc_z)) return
  
  ! X-direction
  if (imp_conv_x .or. imp_visc_x) then

     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)
     
     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
              Ax(j,k,i,:) = 0.0_WP
              Ax(j,k,i,0) = RHOi
              Rx(j,k,i)   = RHOi * ResU(i,j,k)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_x) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=max(jmin_,jmin+icyl),jmax_ ! No convective implicit correction at the axis => tight coupling x-z
              do i=imin_,imax_
                 do st=-stc2,+stc1
                    n = interp_xx(i,j,st)
                    Ax(j,k,i,st-stc1:st+stc2) = Ax(j,k,i,st-stc1:st+stc2) + &
                         0.5_WP*dt2*divc_xx(i,j,st)*(rhoU(i+st+n+1,j,k)+rhoU(i+st-n,j,k))*interp_Ju_xm(i+st,j,:)
                    Ax(j,k,i,st+n+1) = Ax(j,k,i,st+n+1) + 0.5_WP*dt2*divc_xx(i,j,st)*rhoUi(i+st,j,k)
                    Ax(j,k,i,st-n)   = Ax(j,k,i,st-n)   + 0.5_WP*dt2*divc_xx(i,j,st)*rhoUi(i+st,j,k)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_x) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv2,stv1
                    do n=-stv1,stv2
                       Ax(j,k,i,st+n) = Ax(j,k,i,st+n) - &
                            dt2*divv_xx(i,j,st)*2.0_WP*VISC(i+st,j,k)*(grad_u_x(i+st,j,n)-1.0_WP/3.0_WP*divv_u(i+st,j,n))
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_x(ndx)
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Rx(j,k,i) = ResU(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Y-direction
  if (imp_conv_y .or. imp_visc_y) then
     
     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)

     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
              Ay(i,k,j,:) = 0.0_WP
              Ay(i,k,j,0) = RHOi
              Ry(i,k,j)   = RHOi * Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_y) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc1,+stc2
                    n = interp_xy(i,j,st)
                    Ay(i,k,j,st+n-1) = Ay(i,k,j,st+n-1) + 0.5_WP*dt2*divc_xy(i,j,st)*rhoVi(i,j+st,k)
                    Ay(i,k,j,st-n)   = Ay(i,k,j,st-n)   + 0.5_WP*dt2*divc_xy(i,j,st)*rhoVi(i,j+st,k)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_y) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv1,stv2
                    VISCi = sum(interp_sc_xy(i,j+st,:,:)*VISC(i-st2:i+st1,j+st-st2:j+st+st1,k))
                    do n=-stv2,stv1
                       Ay(i,k,j,st+n) = Ay(i,k,j,st+n) - dt2*divv_xy(i,j,st)*VISCi*grad_u_y(i,j+st,n)
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Boundary conditions
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP DO
        do k=kmin_,kmax_
           do n=0,(ndy-1)/2-1
              do st=-(ndy-1)/2,-n-1
                 Ay(:,k,jmin+n,-n) = Ay(:,k,jmin+n,-n) + Ay(:,k,jmin+n,st)
                 Ay(:,k,jmin+n,st) = 0.0_WP
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_y(ndy)
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Ry(i,k,j) = Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Z-direction
  if (imp_conv_z .or. imp_visc_z) then
     
     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)

     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
              Az(i,j,k,:) = 0.0_WP
              Az(i,j,k,0) = RHOi
              Rz(i,j,k)   = RHOi * Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_z) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc1,+stc2
                    n = interp_xz(i,j,st)
                    Az(i,j,k,st+n-1) = Az(i,j,k,st+n-1) + 0.5_WP*dt2*divc_xz(i,j,st)*rhoWi(i,j,k+st)
                    Az(i,j,k,st-n)   = Az(i,j,k,st-n)   + 0.5_WP*dt2*divc_xz(i,j,st)*rhoWi(i,j,k+st)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_z) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv1,stv2
                    VISCi = sum(interp_sc_xz(i,j,:,:)*VISC(i-st2:i+st1,j,k+st-st2:k+st+st1))
                    do n=-stv2,stv1
                       Az(i,j,k,st+n) = Az(i,j,k,st+n) - dt2*divv_xz(i,j,st)*VISCi*grad_u_z(i,j,n)
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_z(ndz)
     
     ! Get back the residual
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ResU(i,j,k) = Rz(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ResU(i,j,k) = Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine velocity_inverse_u


! ================================== !
! Inverse the Residual of v          !
! By using approximate factorization !
! ================================== !
subroutine velocity_inverse_v
  use implicit
  use data
  use metric_generic
  use metric_velocity_conv
  use metric_velocity_visc
  use memory
  use time_info
  implicit none
  
  integer  :: i,j,k,st,n
  real(WP) :: RHOi,VISCi
  real(WP) :: Acos,Asin
  real(WP) :: dt2
  
  dt2 = dt_uvw/2.0_WP
  
  ! If purely explicit return
  if (.not.(imp_conv_x .or. imp_visc_x .or. &
            imp_conv_y .or. imp_visc_y .or. &
            imp_conv_z .or. imp_visc_z)) return
  
  ! X-direction
  if (imp_conv_x .or. imp_visc_x) then

     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)
     
     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
              Ax(j,k,i,:) = 0.0_WP
              Ax(j,k,i,0) = RHOi
              Rx(j,k,i)   = RHOi * ResV(i,j,k)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_x) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc1,+stc2
                    n = interp_yx(i,j,st)
                    Ax(j,k,i,st+n-1) = Ax(j,k,i,st+n-1) + 0.5_WP*dt2*divc_yx(i,j,st)*rhoUi(i+st,j,k)
                    Ax(j,k,i,st-n)   = Ax(j,k,i,st-n)   + 0.5_WP*dt2*divc_yx(i,j,st)*rhoUi(i+st,j,k)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_x) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv1,stv2
                    VISCi = sum(interp_sc_xy(i+st,j,:,:)*VISC(i+st-st2:i+st+st1,j-st2:j+st1,k))
                    do n=-stv2,stv1
                       Ax(j,k,i,st+n) = Ax(j,k,i,st+n) - dt2*divv_yx(i,j,st)*VISCi*grad_v_x(i+st,j,n)
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_x(ndx)
     
     ! Axis treatment
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP PARALLEL DO PRIVATE(Acos,Asin)
        do i=imin_,imax_
           Acos = 2.0_WP*sum(Rx(jmin,:,i)*cos(zm(kmin:kmax)))/real(nz,WP)
           Asin = 2.0_WP*sum(Rx(jmin,:,i)*sin(zm(kmin:kmax)))/real(nz,WP)
           Rx(jmin,:,i) = Acos*cos(zm(kmin:kmax)) + Asin*sin(zm(kmin:kmax))
        end do
        !$OMP END PARALLEL DO
     end if
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Rx(j,k,i) = ResV(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Y-direction
  if (imp_conv_y .or. imp_visc_y) then
     
     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)

     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
              Ay(i,k,j,:) = 0.0_WP
              Ay(i,k,j,0) = RHOi
              Ry(i,k,j)   = RHOi * Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_y) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc2,+stc1
                    n = interp_yy(i,j,st)
                    Ay(i,k,j,st-stc1:st+stc2) = Ay(i,k,j,st-stc1:st+stc2) + &
                         0.5_WP*dt2*divc_yy(i,j,st)*(rhoV(i,j+st+n+1,k)+rhoV(i,j+st-n,k))*interp_Jv_ym(i,j+st,:)
                    Ay(i,k,j,st+n+1) = Ay(i,k,j,st+n+1) + 0.5_WP*dt2*divc_yy(i,j,st)*rhoVi(i,j+st,k)
                    Ay(i,k,j,st-n)   = Ay(i,k,j,st-n)   + 0.5_WP*dt2*divc_yy(i,j,st)*rhoVi(i,j+st,k)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_y) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv2,stv1
                    do n=-stv1,stv2
                       Ay(i,k,j,st+n) = Ay(i,k,j,st+n) - dt2*2.0_WP*VISC(i,j+st,k)*( &
                            +divv_yy(i,j,st)*(grad_v_y(i,j+st,n)-1.0_WP/3.0_WP*divv_v(i,j+st,n)) &
                            -yi(j)*interpv_cyl_F_y(i,j,st)*( ymi(j+st)*interpv_cyl_v_ym(i,j+st,n)-1.0_WP/3.0_WP*divv_v(i,j+st,n)) )
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Boundary conditions
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP DO
        do k=kmin_,kmax_
           do n=0,(ndy-1)/2-1
              do st=-(ndy-1)/2,-n-1
                 Ay(:,k,jmin+n,-n) = Ay(:,k,jmin+n,-n) + Ay(:,k,jmin+n,st)
                 Ay(:,k,jmin+n,st) = 0.0_WP
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_y(ndy)
     
     ! Axis treatment
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP PARALLEL DO PRIVATE(Acos,Asin)
        do i=imin_,imax_
           Acos = 2.0_WP*sum(Ry(i,:,jmin)*cos(zm(kmin:kmax)))/real(nz,WP)
           Asin = 2.0_WP*sum(Ry(i,:,jmin)*sin(zm(kmin:kmax)))/real(nz,WP)
           Ry(i,:,jmin) = Acos*cos(zm(kmin:kmax)) + Asin*sin(zm(kmin:kmax))
        end do
        !$OMP END PARALLEL DO
     end if
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Ry(i,k,j) = Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Z-direction
  if (imp_conv_z .or. imp_visc_z) then
     
     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n,Acos,Asin)

     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
              Az(i,j,k,:) = 0.0_WP
              Az(i,j,k,0) = RHOi
              Rz(i,j,k)   = RHOi * Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_z) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc1,+stc2
                    n = interp_yz(i,j,st)
                    Az(i,j,k,st+n-1) = Az(i,j,k,st+n-1) + 0.5_WP*dt2*divc_yz(i,j,st)*rhoWi(i,j,k+st)
                    Az(i,j,k,st-n)   = Az(i,j,k,st-n)   + 0.5_WP*dt2*divc_yz(i,j,st)*rhoWi(i,j,k+st)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_z) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv1,stv2
                    VISCi = sum(interp_sc_yz(i,j,:,:)*VISC(i,j-st2:j+st1,k+st-st2:k+st+st1))
                    do n=-stv2,stv1
                       Az(i,j,k,st+n) = Az(i,j,k,st+n) - dt2*divv_yz(i,j,st)*VISCi*grad_v_z(i,j,n)
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Morinishi's axis treatment
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP DO
        do i=imin_,imax_
           Acos = 2.0_WP*sum(ResV(i,jmin,:)*cos(zm(kmin:kmax)))/real(nz,WP)
           Asin = 2.0_WP*sum(ResV(i,jmin,:)*sin(zm(kmin:kmax)))/real(nz,WP)
           ResV(i,jmin,:) = Acos*cos(zm(kmin:kmax)) + Asin*sin(zm(kmin:kmax))
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_z(ndz)
     
     ! Get back the residual
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ResV(i,j,k) = Rz(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ResV(i,j,k) = Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine velocity_inverse_v


! ================================== !
! Inverse the Residual of w          !
! By using approximate factorization !
! ================================== !
subroutine velocity_inverse_w
  use implicit
  use data
  use metric_generic
  use metric_velocity_conv
  use metric_velocity_visc
  use memory
  use time_info
  implicit none
  
  integer  :: i,j,k,st,n
  real(WP) :: RHOi,VISCi
  real(WP) :: dt2
  
  dt2 = dt_uvw/2.0_WP
  
  ! If purely explicit return
  if (.not.(imp_conv_x .or. imp_visc_x .or. &
            imp_conv_y .or. imp_visc_y .or. &
            imp_conv_z .or. imp_visc_z)) return
  
  ! X-direction
  if (imp_conv_x .or. imp_visc_x) then

     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)
     
     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
              Ax(j,k,i,:) = 0.0_WP
              Ax(j,k,i,0) = RHOi
              Rx(j,k,i)   = RHOi * ResW(i,j,k)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_x) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc1,+stc2
                    n = interp_zx(i,j,st)
                    Ax(j,k,i,st+n-1) = Ax(j,k,i,st+n-1) + 0.5_WP*dt2*divc_zx(i,j,st)*rhoUi(i+st,j,k)
                    Ax(j,k,i,st-n)   = Ax(j,k,i,st-n)   + 0.5_WP*dt2*divc_zx(i,j,st)*rhoUi(i+st,j,k)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_x) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv1,stv2
                    VISCi = sum(interp_sc_xz(i+st,j,:,:)*VISC(i+st-st2:i+st+st1,j,k-st2:k+st1))
                    do n=-stv2,stv1
                       Ax(j,k,i,st+n) = Ax(j,k,i,st+n) - dt2*divv_zx(i,j,st)*VISCi*grad_w_x(i+st,j,n)
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_x(ndx)
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Rx(j,k,i) = ResW(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Y-direction
  if (imp_conv_y .or. imp_visc_y) then

     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)
     
     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
              Ay(i,k,j,:) = 0.0_WP
              Ay(i,k,j,0) = RHOi
              Ry(i,k,j)   = RHOi * Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_y) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc1,+stc2
                    n = interp_zy(i,j,st)
                    Ay(i,k,j,st+n-1) = Ay(i,k,j,st+n-1) + 0.5_WP*dt2*divc_zy(i,j,st)*rhoVi(i,j+st,k)
                    Ay(i,k,j,st-n)   = Ay(i,k,j,st-n)   + 0.5_WP*dt2*divc_zy(i,j,st)*rhoVi(i,j+st,k)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_y) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv1,stv2
                    VISCi = sum(interp_sc_yz(i,j+st,:,:)*VISC(i,j+st-st2:j+st+st1,k-st2:k+st1))
                    do n=-stv2,stv1
                       Ay(i,k,j,st+n) = Ay(i,k,j,st+n) - dt2*VISCi*( &
                            +divv_zy(i,j,st)*grad_w_y(i,j+st,n) &
                            +ymi(j)*interpv_cyl_F_ym(i,j,st)*(grad_w_y(i,j+st,n)-yi(j+st)*interpv_cyl_w_y(i,j+st,n)) )
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Boundary conditions
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP DO
        do k=kmin_,kmax_
           do n=0,(ndy-1)/2-1
              do st=-(ndy-1)/2,-n-1
                 Ay(:,k,jmin+n,-n) = Ay(:,k,jmin+n,-n) + Ay(:,k,jmin+n,st)
                 Ay(:,k,jmin+n,st) = 0.0_WP
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_y(ndy)
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Ry(i,k,j) = Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Z-direction
  if (imp_conv_z .or. imp_visc_z) then

     !$OMP PARALLEL PRIVATE(RHOi,VISCi,n)
     
     ! Set matrix diagonal and RHS
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHOi = sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
              Az(i,j,k,:) = 0.0_WP
              Az(i,j,k,0) = RHOi
              Rz(i,j,k)   = RHOi * Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Convective part
     if (imp_conv_z) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stc2,+stc1
                    n = interp_zz(i,j,st)
                    Az(i,j,k,st-stc1:st+stc2) = Az(i,j,k,st-stc1:st+stc2) + &
                         0.5_WP*dt2*divc_zz(i,j,st)*(rhoW(i,j,k+st+n+1)+rhoW(i,j,k+st-n))*interp_Jw_zm(i,j,:) + &
                         dt2*ymi(j)*interp_cyl_F_z(i,j,st)*sum(interp_cyl_v_ym(i,j,:)*rhoV(i,j-stc1:j+stc2,k+st))*interp_cyl_w_zm(i,j,:)
                    Az(i,j,k,st+n+1) = Az(i,j,k,st+n+1) + 0.5_WP*dt2*divc_zz(i,j,st)*rhoWi(i,j,k+st)
                    Az(i,j,k,st-n)   = Az(i,j,k,st-n)   + 0.5_WP*dt2*divc_zz(i,j,st)*rhoWi(i,j,k+st)
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if
     
     ! Viscous part
     if (imp_visc_z) then
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 do st=-stv2,stv1
                    do n=-stv1,stv2
                       Az(i,j,k,st+n) = Az(i,j,k,st+n) - &
                            dt2*divv_zz(i,j,st)*2.0_WP*VISC(i,j,k+st)*(grad_w_z(i,j,n)-1.0_WP/3.0_WP*divv_w(i,j,n))
                    end do
                 end do
              end do
           end do
        end do
        !$OMP END DO
     end if

     !$OMP END PARALLEL
     
     ! Solve
     call implicit_solve_z(ndz)
     
     ! Get back the residual
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ResW(i,j,k) = Rz(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     
  else
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ResW(i,j,k) = Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine velocity_inverse_w


! ============================ !
! Monitor the implicit solvers !
! ============================ !
subroutine implicit_monitor
  use implicit
  use time_info
  use data
  implicit none
  
  ! Compute the CFL
  call velocity_CFL
  
  ! Transfer the values to monitor
  call monitor_select_file('timestep')
  call monitor_set_single_value(1,dt)
  call monitor_set_single_value(2,CFL)
  call monitor_set_single_value(3,CFLc_x)
  call monitor_set_single_value(4,CFLc_y)
  call monitor_set_single_value(5,CFLc_z)
  call monitor_set_single_value(6,CFLc_axis)
  call monitor_set_single_value(7,CFLv_x)
  call monitor_set_single_value(8,CFLv_y)
  call monitor_set_single_value(9,CFLv_z)
  call monitor_set_single_value(10,CFLv_axis)
  if (nscalar.ne.0) then
     call monitor_set_single_value(11,CFLd_x)
     call monitor_set_single_value(12,CFLd_y)
     call monitor_set_single_value(13,CFLd_z)
     call monitor_set_single_value(14,CFLd_axis)     
  end if
  if (compressible) then
     if (nscalar.eq.0) then
        call monitor_set_single_value(11,CFLa_x)
        call monitor_set_single_value(12,CFLa_y)
        call monitor_set_single_value(13,CFLa_z)
        call monitor_set_single_value(14,CFLa_axis)
     else
        call monitor_set_single_value(15,CFLa_x)
        call monitor_set_single_value(16,CFLa_y)
        call monitor_set_single_value(17,CFLa_z)
        call monitor_set_single_value(18,CFLa_axis)
     end if
  end if
  
  return
end subroutine implicit_monitor


subroutine implicit_solve_x(nd)
  use precision
  use memory
  use partition
  implicit none
  
  ! Number of diagonals
  integer, intent(in) :: nd
  real(WP), dimension(:,:), allocatable :: A,B,C,D,E
  real(WP), dimension(:,:,:), allocatable :: P
  integer :: i,j,k,n
  
  ! Choose based on number of diagonals
  select case(nd)
  case(3)
     allocate(A(ny_*nz_,nx_))
     allocate(B(ny_*nz_,nx_))
     allocate(C(ny_*nz_,nx_))
     !$OMP PARALLEL DO
     do i=imin_,imax_
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              A((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i,-1)
              B((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i, 0)
              C((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i, 1)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     ! LAST TWO: ny*nz x nx
     call tridiagonal(A,B,C,Rx,nx_,ny_*nz_,'x')
     deallocate(A); deallocate(B); deallocate(C)
  case(5)
     allocate(A(ny_*nz_,nx_)); allocate(B(ny_*nz_,nx_)); allocate(C(ny_*nz_,nx_)); allocate(D(ny_*nz_,nx_)); allocate(E(ny_*nz_,nx_))
     !$OMP PARALLEL DO
     do i=imin_,imax_
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              A((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i,-2)
              B((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i,-1)
              C((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i, 0)
              D((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i, 1)
              E((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1) = Ax(j,k,i, 2)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     ! LAST TWO: ny*nz x nx
     call pentadiagonal(A,B,C,D,E,Rx,nx_,ny_*nz_,'x')
     deallocate(A); deallocate(B); deallocate(C); deallocate(D); deallocate(E)
  case default
     allocate(P(ny_*nz_,nx_,-(nd-1)/2:(nd-1)/2))
     !$OMP PARALLEL
     do n=-(nd-1)/2,(nd-1)/2
        !$OMP DO
        do i=imin_,imax_
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 P((k-kmin_)*ny_+(j-jmin_)+1,i-imin_+1,n) = Ax(j,k,i,n)
              end do
           end do
        end do
        !$OMP END DO
     end do
     !$OMP END PARALLEL
     !LAST ONE: ny*nz x nx x nd
     call polydiagonal((nd-1)/2,P,Rx,nx_,ny_*nz_,'x')
     deallocate(P)
  end select
  
  return
end subroutine implicit_solve_x


subroutine implicit_solve_y(nd)
  use precision
  use memory
  use partition
  use parallel
  implicit none
  
  ! Number of diagonals
  integer, intent(in) :: nd
  real(WP), dimension(:,:), allocatable :: A,B,C,D,E
  real(WP), dimension(:,:,:), allocatable :: P
  integer :: i,j,k,n
  
  ! Choose based on number of diagonals
  select case(nd)
  case(3)
     allocate(A(nx_*nz_,ny_)); allocate(B(nx_*nz_,ny_)); allocate(C(nx_*nz_,ny_))
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do k=kmin_,kmax_
           do i=imin_,imax_
              A((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j,-1)
              B((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j, 0)
              C((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j, 1)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call tridiagonal(A,B,C,Ry,ny_,nx_*nz_,'y')
     deallocate(A); deallocate(B); deallocate(C)
  case(5)
     allocate(A(nx_*nz_,ny_)); allocate(B(nx_*nz_,ny_)); allocate(C(nx_*nz_,ny_)); allocate(D(nx_*nz_,ny_)); allocate(E(nx_*nz_,ny_))
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do k=kmin_,kmax_
           do i=imin_,imax_
              A((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j,-2)
              B((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j,-1)
              C((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j, 0)
              D((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j, 1)
              E((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1) = Ay(i,k,j, 2)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call pentadiagonal(A,B,C,D,E,Ry,ny_,nx_*nz_,'y')
     deallocate(A); deallocate(B); deallocate(C); deallocate(D); deallocate(E)
  case default
     allocate(P(nx_*nz_,ny_,-(nd-1)/2:(nd-1)/2))
     !$OMP PARALLEL
     do n=-(nd-1)/2,(nd-1)/2
        !$OMP DO
        do j=jmin_,jmax_
           do k=kmin_,kmax_
              do i=imin_,imax_
                 P((k-kmin_)*nx_+(i-imin_)+1,j-jmin_+1,n) = Ay(i,k,j,n)
              end do
           end do
        end do
        !$OMP END DO
     end do
     !$OMP END PARALLEL
     call polydiagonal((nd-1)/2,P,Ry,ny_,nx_*nz_,'y')
     deallocate(P)
  end select
  
  return
end subroutine implicit_solve_y


subroutine implicit_solve_z(nd)
  use precision
  use memory
  use partition
  implicit none
  
  ! Number of diagonals
  integer, intent(in) :: nd
  real(WP), dimension(:,:), allocatable :: A,B,C,D,E
  real(WP), dimension(:,:,:), allocatable :: P
  integer :: i,j,k,n
  
  ! Choose based on number of diagonals
  select case(nd)
  case(3)
     allocate(A(nx_*ny_,nz_)); allocate(B(nx_*ny_,nz_)); allocate(C(nx_*ny_,nz_))
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              A((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k,-1)
              B((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k, 0)
              C((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k, 1)
           end do
        end do
     end do
     call tridiagonal(A,B,C,Rz,nz_,nx_*ny_,'z')
     deallocate(A); deallocate(B); deallocate(C)
  case(5)
     allocate(A(nx_*ny_,nz_)); allocate(B(nx_*ny_,nz_)); allocate(C(nx_*ny_,nz_)); allocate(D(nx_*ny_,nz_)); allocate(E(nx_*ny_,nz_))
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              A((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k,-2)
              B((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k,-1)
              C((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k, 0)
              D((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k, 1)
              E((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1) = Az(i,j,k, 2)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call pentadiagonal(A,B,C,D,E,Rz,nz_,nx_*ny_,'z')
     deallocate(A); deallocate(B); deallocate(C); deallocate(D); deallocate(E)
  case default
     allocate(P(nx_*ny_,nz_,-(nd-1)/2:(nd-1)/2))
     !$OMP PARALLEL
     do n=-(nd-1)/2,(nd-1)/2
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 P((j-jmin_)*nx_+(i-imin_)+1,k-kmin_+1,n) = Az(i,j,k,n)
              end do
           end do
        end do
        !$OMP END DO
     end do
     !$OMP END PARALLEL
     call polydiagonal((nd-1)/2,P,Rz,nz_,nx_*ny_,'z')
     deallocate(P)
  end select
  
  return
end subroutine implicit_solve_z
