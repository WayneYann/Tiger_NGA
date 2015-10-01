module sgs_lagrangian
  use sgsmodel
  use optdata
  implicit none
  
  ! Indices for OptData
  integer :: iLM_vel,iMM_vel
  integer :: iLM_var,iMM_var
  integer, dimension(:), pointer :: iLM_sc,iMM_sc
  logical :: sgs_vel_present,sgs_sc_present,sgs_var_present
  
  ! Start with Germano?
  logical :: precomp_vel,precomp_var
  logical, dimension(:), pointer :: precomp_sc
  
  ! Numerator and denominator for the Dynamic Smagorinsky
  real(WP), dimension(:,:,:),   pointer :: LM_vel,MM_vel
  real(WP), dimension(:,:,:,:), pointer :: LM_sc, MM_sc
  real(WP), dimension(:,:,:),   pointer :: LM_var,MM_var
  
  ! For "clipping" purposes
  real(WP), parameter :: eps_LM = 1.0E-34_WP
  real(WP), parameter :: Cs_typ = 0.16_WP
  
contains
  
  ! Get the indices of variables from their names
  ! ---------------------------------------------
  subroutine get_indices
    implicit none
    character(len=str_short) :: name,var
    integer :: isc,iod
    
    allocate(iLM_sc(nscalar))
    allocate(iMM_sc(nscalar))
    
    iLM_vel = -1
    iMM_vel = -1
    iLM_sc  = -1
    iMM_sc  = -1
    iLM_var = -1
    iMM_var = -1
    
    sgs_vel_present = .false.
    sgs_sc_present  = .false.
    sgs_var_present = .false.
    
    do iod=1,nod
       read(OD_name(iod),*) name
       
       ! Numerator  - LM
       if (name(1:3).eq.'LM_') then
          read(name(4:str_short),*) var
          select case(trim(var))
          case ('VEL')
             iLM_vel = iod
             sgs_vel_present = .true.
          case ('VAR')
             iLM_var = iod
             sgs_var_present = .true.
          case default
             loop1:do isc=1,nscalar
                if (trim(var).eq.trim(SC_name(isc))) exit loop1
             end do loop1
             if (isc.gt.nscalar) call die('sgs_lagrangian_get_indices: unknown scalar')
             iLM_sc(isc) = iod
             sgs_sc_present = .true.
          end select
       end if

       ! Denominator  - MM
       if (name(1:3).eq.'MM_') then
          read(name(4:str_short),*) var
          select case(trim(var))
          case ('VEL')
             iMM_vel = iod
             sgs_vel_present = .true.
          case ('VAR')
             iMM_var = iod
             sgs_var_present = .true.
          case default
             loop2:do isc=1,nscalar
                if (trim(var).eq.trim(SC_name(isc))) exit loop2
             end do loop2
             if (isc.gt.nscalar) call die('sgs_lagrangian_get_indices: unknown scalar')
             iMM_sc(isc) = iod
             sgs_sc_present = .true.
          end select
       end if
       
    end do
    
    ! Test if we have both num & den for velocity
    if (sgs_vel_present) then
       if (iLM_vel.eq.-1) call die('sgs_lagrangian_get_indices: needs both LM and MM for velocity')
       if (iMM_vel.eq.-1) call die('sgs_lagrangian_get_indices: needs both LM and MM for velocity')
    end if
    
    ! Test if we have all or none of the scalar num & den
    if (sgs_sc_present) then
       do isc=1,nscalar
          if (iLM_sc(isc).eq.-1) call die('sgs_lagrangian_get_indices: missing scalar quantity')
          if (iMM_sc(isc).eq.-1) call die('sgs_lagrangian_get_indices: missing scalar quantity')
       end do
    end if
    
    ! Test if we have both num & den for variance
    if (sgs_var_present) then
       if (iLM_var.eq.-1) call die('sgs_lagrangian_get_indices: needs both LM and MM for variance')
       if (iMM_var.eq.-1) call die('sgs_lagrangian_get_indices: needs both LM and MM for variance')
    end if

    return
  end subroutine get_indices

  ! Interpolate LM and MM at the given position
  ! -------------------------------------------
  subroutine sgs_lagrangian_interpolate(xp,yp,zp,A,Ap)
    implicit none
    
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
    real(WP), intent(inout) :: xp,yp,zp
    real(WP), intent(out)   :: Ap
    real(WP) :: wx1,wy1,wz1,wx2,wy2,wz2
    integer :: i,j,k
    
    ! Periodicity in z
    zp = mod(zp-z(kmin),zL) + z(kmin)
    if (zp.lt.z(kmin)) zp = zp + zL
    
    ! Find the index for the interpolation
    ! Compute the weights
    if (xp.gt.xm(imaxo_)) then
       i = imaxo_
    else
       i = imino_+1
       do while (xm(i).lt.xp)
          i = i+1
       end do
    end if
    wx1 = (xp-xm(i-1))/(xm(i)-xm(i-1))
    wx2 = 1.0_WP - wx1
    
    if (yp.gt.ym(jmaxo_)) then
       j = jmaxo_
    else
       j = jmino_+1
       do while (ym(j).lt.yp)
          j = j+1
       end do
    end if
    wy1 = (yp-ym(j-1))/(ym(j)-ym(j-1))
    wy2 = 1.0_WP - wy1
    
    if (zp.gt.zm(kmaxo_)) then
       k = kmaxo_
    else
       k = kmino_+1
       do while (zm(k).lt.zp)
          k = k+1
       end do
    end if
    wz1 = (zp-zm(k-1))/(zm(k)-zm(k-1))
    wz2 = 1.0_WP - wz1
    
    ! Interpolate
    Ap = wz1*(wy1*(wx1*A(i,j,k)     + wx2*A(i-1,j,k))    + &
              wy2*(wx1*A(i,j-1,k)   + wx2*A(i-1,j-1,k))) + &
         wz2*(wy1*(wx1*A(i,j,k-1)   + wx2*A(i-1,j,k-1))  + &
              wy2*(wx1*A(i,j-1,k-1) + wx2*A(i-1,j-1,k-1)))
    
    return
  end subroutine sgs_lagrangian_interpolate
  
  
  ! Perform the average in time
  ! ---------------------------
  subroutine sgs_lagrangian_average(LM_new,MM_new,LM_old,MM_old)
    use interpolate
    use filter
    use time_info
    use masks
    implicit none
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: LM_new,MM_new
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: LM_old,MM_old
    real(WP) :: xp,yp,zp,MM_p,LM_p,alpha,tau
    integer  :: i,j,k
    
    !$OMP PARALLEL PRIVATE(xp,yp,zp,LM_p,MM_p,tau,alpha)

    !$OMP DO
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             
             ! Lagrangian back tracking
             xp = xm(i) - Ui(i,j,k)*dt
             yp = ym(j) - Vi(i,j,k)*dt
             if (icyl.eq.0) then
                zp = zm(k) - Wi(i,j,k)*dt
             else
                zp = zm(k) - 2.0_WP*Wi(i,j,k)*dt/(yp+ym(j))
             end if
             
             ! Interpolation at old location
             call sgs_lagrangian_interpolate(xp,yp,zp,LM_old,LM_p)
             call sgs_lagrangian_interpolate(xp,yp,zp,MM_old,MM_p)
           
             ! Advance LM and MM
             tau = dt * (LM_vel(i,j,k)*MM_vel(i,j,k))**0.125_WP / (1.5_WP*delta_3D(i,j))
             alpha = tau / (1.0_WP + tau)
             LM_new(i,j,k) = alpha * LM_new(i,j,k) + (1.0_WP-alpha) * LM_p
             MM_new(i,j,k) = alpha * MM_new(i,j,k) + (1.0_WP-alpha) * MM_p
             if (LM_new(i,j,k).lt.eps_LM) LM_new(i,j,k) = eps_LM
             if (MM_new(i,j,k).lt.eps_LM/Cs_typ**2) MM_new(i,j,k) = eps_LM/Cs_typ**2
          end do
       end do
    end do
    !$OMP END DO
    
    ! Treat the inlet if necessary
    if ((iproc.eq.1) .and. (xper.eq.0)) then
       !$OMP DO
       do j=jmin_,jmax_
          LM_new(imin,j,:) = sum(LM_new(imin,j,kmin_:kmax_))/real(nz_,WP)
          MM_new(imin,j,:) = sum(MM_new(imin,j,kmin_:kmax_))/real(nz_,WP)
       end do
       !$OMP END DO
    end if

    !$OMP END PARALLEL
    
    ! Apply BCs
    call sgsmodel_apply_border(LM_new)
    call sgsmodel_apply_border(MM_new)
    
    ! Enforce non zero values
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             if (mask(i,j).ne.0 .and. mask(i,j).ne.2) then
                LM_new(i,j,k) = eps_LM
                MM_new(i,j,k) = eps_LM/Cs_typ**2
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine sgs_lagrangian_average
  
end module sgs_lagrangian


! =================================== !
! Initialize the Lagrangian SGS model !
! =================================== !
subroutine sgs_lagrangian_init
  use sgs_lagrangian
  use parser
  implicit none
  
  integer :: isc,i,j,k
  logical :: precomp
  
  ! Allocate the precomp arrays
  allocate(precomp_sc(nscalar))

  ! Get the indices from OptData
  call get_indices
  
  ! Viscosity
  allocate(LM_vel(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(MM_vel(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  if (.not.sgs_vel_present) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              LM_vel(i,j,k) = 0.0_WP
              MM_vel(i,j,k) = 0.0_WP
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     precomp_vel = .true.
  else
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              LM_vel(i,j,k) = OD(i,j,k,iLM_vel)
              MM_vel(i,j,k) = OD(i,j,k,iMM_vel)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call sgsmodel_apply_border(LM_vel)
     call sgsmodel_apply_border(MM_vel)
     call monitor_log('SGS_LAGRANGIAN: VELOCITY LM & MM DETECTED')
     precomp_vel = .false.
  end if
  
  if (nscalar.ge.1) then
     
     ! Diffusivity
     if (.not.sgs_sc_present) then
        allocate(LM_sc(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
        allocate(MM_sc(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))

        !$OMP PARALLEL

        do isc=1,nscalar
           !$OMP DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    LM_sc(i,j,k,isc) = 0.0_WP
                    MM_sc(i,j,k,isc) = 0.0_WP
                 end do
              end do
           end do
           !$OMP END DO
        end do

        !$OMP END PARALLEL

        precomp_sc = .true.
     else
        do isc=1,nscalar
           call sgsmodel_apply_border(OD(:,:,:,iLM_sc(isc)))
           call sgsmodel_apply_border(OD(:,:,:,iMM_sc(isc)))
        end do
        call monitor_log('SGS_LAGRANGIAN: SCALAR LM & MM DETECTED')
        precomp_sc = .false.
     end if
     
     ! For variance
     if (.not.sgs_var_present) then
        allocate(LM_var(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
        allocate(MM_var(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 LM_var(i,j,k) = 0.0_WP
                 MM_var(i,j,k) = 0.0_WP
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        precomp_var = .true.
     else
        call sgsmodel_apply_border(OD(:,:,:,iLM_var))
        call sgsmodel_apply_border(OD(:,:,:,iMM_var))
        call monitor_log('SGS_LAGRANGIAN: VARIANCE LM & MM DETECTED')
        precomp_var = .false.
     end if

  end if
  
  ! Override the initial field with Germano
  call parser_read('SGS Override Lag.',precomp,.true.)
  if (precomp) then
     precomp_vel = .true.
     precomp_sc  = .true.
     precomp_var = .true.
  end if

  return
end subroutine sgs_lagrangian_init


! ============================================== !
! Average LM and MM over time for Eddy Viscosity !
! ============================================== !
subroutine sgs_lagrangian_average_VISC
  use sgs_lagrangian
  implicit none

  integer :: i,j,k
  
  ! Compute the average with Germano if necessary
  if (precomp_vel) then
     call sgs_germano_average_VISC
     where (LM.lt.eps_LM) LM = eps_LM
     precomp_vel = .false.
  else
     call sgs_lagrangian_average(LM,MM,LM_vel,MM_vel)
  end if
  
  ! Save data for next time
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           LM_vel(i,j,k) = LM(i,j,k)
           MM_vel(i,j,k) = MM(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  if (sgs_vel_present) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              OD(i,j,k,iLM_vel) = LM(i,j,k)
              OD(i,j,k,iMM_vel) = MM(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine sgs_lagrangian_average_VISC


! ================================================ !
! Average LM and MM over time for Eddy Diffusivity !
! ================================================ !
subroutine sgs_lagrangian_average_DIFF(isc)
  use sgs_lagrangian
  implicit none
  integer, intent(in) :: isc
  integer :: i,j,k
  
  ! Compute the average with Germano if necessary
  if (precomp_sc(isc)) then
     call sgs_germano_average_DIFF
     precomp_sc(isc) = .false.
  else
     
     if (sgs_sc_present) then
        call sgs_lagrangian_average(LM,MM,OD(:,:,:,iLM_sc(isc)),OD(:,:,:,iMM_sc(isc)))
     else
        call sgs_lagrangian_average(LM,MM,LM_sc(:,:,:,isc),MM_sc(:,:,:,isc))
     end if
  end if
  
  ! Save data for next time
  if (sgs_sc_present) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              OD(i,j,k,iLM_sc(isc)) = LM(i,j,k)
              OD(i,j,k,iMM_sc(isc)) = MM(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              LM_sc(i,j,k,isc) = LM(i,j,k)
              MM_sc(i,j,k,isc) = MM(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine sgs_lagrangian_average_DIFF


! =============================================== !
! Average LM and MM over time for Scalar Variance !
! =============================================== !
subroutine sgs_lagrangian_average_ZVAR
  use sgs_lagrangian
  implicit none

  integer :: i,j,k
  
  ! Compute the average with Germano if necessary
  if (precomp_var) then
     call sgs_germano_average_ZVAR
     precomp_var = .false.
  else
     if (sgs_var_present) then
        call sgs_lagrangian_average(LM,MM,OD(:,:,:,iLM_var),OD(:,:,:,iMM_var))
     else
        call sgs_lagrangian_average(LM,MM,LM_var,MM_var)
     end if
  end if
  
  ! Save data for next time
  if (sgs_var_present) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              OD(i,j,k,iLM_var) = LM(i,j,k)
              OD(i,j,k,iMM_var) = MM(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              LM_var(i,j,k) = LM(i,j,k)
              MM_var(i,j,k) = MM(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine sgs_lagrangian_average_ZVAR
