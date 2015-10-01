module boundary
  use precision
  use geometry
  use partition
  implicit none

  ! Change of mass in the domain
  real(WP) :: masschange

  ! Momentum flux across boundaries
  real(WP) :: massflux,massflux_exit
  
  ! Added mass in the domain
  real(WP) :: massadded
  
  ! Mass correction
  real(WP) :: masscorrection
  
contains

  ! ========================================== !
  ! Compute change in mass flux across all BCs !
  ! ========================================== !
  subroutine boundary_massflux
    use data
    use parallel
    use masks
    implicit none

    real(WP) :: tmp
    integer :: i,j,k

    ! Return if periodic in x
    tmp = 0.0_WP
    massflux_exit = 0.0_WP
    if (xper.eq.1) return
    
    !$OMP PARALLEL

    ! Left BC
    if (iproc.eq.1) then
       !$OMP DO REDUCTION(+:tmp)
       do j=jmin_,jmax_
          do k=kmin_,kmax_
             tmp = tmp + rhoU(imin,j,k)*dA(j)
          end do
       end do
       !$OMP END DO
    end if
    ! Right BC
    if (iproc.eq.npx) then
       !$OMP DO REDUCTION(+:tmp,massflux_exit)
       do j=jmin_,jmax_
          do k=kmin_,kmax_
             tmp = tmp - rhoU(imax+1,j,k)*dA(j)
             massflux_exit = massflux_exit + rhoU(imax+1,j,k)*dA(j)
          end do
       end do
       !$OMP END DO
    end if
    ! Lower BC
    if (yper.ne.1 .and. jproc.eq.1) then
       !$OMP DO REDUCTION(+:tmp)
       do i=imin_,imax_
          do k=kmin_,kmax_
             tmp = tmp + rhoV(i,jmin,k)*dx(i)*dz_v(jmin)
          end do
       end do
       !$OMP END DO
    end if
    ! Upper BC
    if (yper.ne.1 .and. jproc.eq.npy) then
       !$OMP DO REDUCTION(+:tmp)
       do i=imin_,imax_
          do k=kmin_,kmax_
             tmp = tmp - rhoV(i,jmax+1,k)*dx(i)*dz_v(jmax+1)
          end do
       end do
       !$OMP END DO
    end if

    !$OMP END PARALLEL
    
    call parallel_sum(tmp,massflux)
    call parallel_sum(massflux_exit,tmp)
    massflux_exit=tmp
    
    return
  end subroutine boundary_massflux


  ! ======================================== !
  ! Compute mass of change inside the domain !
  ! ======================================== !
  subroutine boundary_masschange
    use combustion
    use time_info
    implicit none
    
    masschange = sum_dRHO/dt
    
    return
  end subroutine boundary_masschange
  
  
  ! ==================================== !
  ! Compute added mass inside the domain !
  ! ==================================== !
  subroutine boundary_massadded
    use velocity
    use time_info
    use parallel
    implicit none
    real(WP) :: tmp,sum_srcP
    integer :: i,j,k
    
    tmp = 0.0_WP
    !$OMP PARALLEL DO REDUCTION(+:tmp)
    do i=imin_,imax_
       do j=jmin_,jmax_
          do k=kmin_,kmax_
             tmp = tmp + srcP(i,j,k)*vol(i,j)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    call parallel_sum(tmp,sum_srcP)
    
    massadded = sum_srcP/dt
    
    return
  end subroutine boundary_massadded
  
  
end module boundary


! ============================== !
! Initialize the Boundary module !
! ============================== !
subroutine boundary_init
  use boundary
  use parallel
  use data
  use masks
  implicit none
  integer :: i,j
  
  ! Default treatment of the walls
  do j=jmino_,jmaxo_
     do i=imino_,imaxo_
        if (mask_u(i,j).eq.1) then
           U(i,j,:) = 0.0_WP
           rhoU(i,j,:) = 0.0_WP
        end if
        if (mask_v(i,j).eq.1) then
           V(i,j,:) = 0.0_WP
           rhoV(i,j,:) = 0.0_WP
        end if
        if (mask_w(i,j).eq.1) then
           W(i,j,:) = 0.0_WP
           rhoW(i,j,:) = 0.0_WP
        end if
        if (mask(i,j).eq.1) then
           if (nscalar.gt.0) SC(i,j,:,:) = 0.0_WP
           RHO(i,j,:) = 1.0_WP
           RHOold(i,j,:) = 1.0_WP
           dRHO(i,j,:) = 0.0_WP
        end if
     end do
  end do
  
  ! Inflow/Outflow if not periodic
  call inflow_init
  call outflow_init

  return
end subroutine boundary_init


! ========================== !
! Apply Dirichlet Conditions !
! ========================== !
subroutine boundary_dirichlet(vec,dir)
  use boundary
  use parallel
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec
  character(len=*), intent(in) :: dir
  integer :: i,j
  
  select case(trim(dir))
  case('-x')
     if (xper.eq.0 .and. iproc.eq.1) then
        do i=imino,imin
           vec(i,:,:) = 0.0_WP
        end do
     end if
  case('+x')
     if (xper.eq.0 .and. iproc.eq.npx) then
        do i=imax+1,imaxo
           vec(i,:,:) = 0.0_WP
        end do
     end if
  case('-y')
     if (yper.eq.0 .and. jproc.eq.1 .and. icyl.eq.0) then
        do j=jmino,jmin
           vec(:,j,:) = 0.0_WP
        end do
     end if
  case('+y')
     if (yper.eq.0 .and. jproc.eq.npy) then
        do j=jmax+1,jmaxo
           vec(:,j,:) = 0.0_WP
        end do
     end if
  case('-xm')
     if (xper.eq.0 .and. iproc.eq.1) then
        do i=imino,imin-1
           vec(i,:,:) = 0.0_WP
        end do
     end if
  case('+xm')
     if (xper.eq.0 .and. iproc.eq.npx) then
        do i=imax+1,imaxo
           vec(i,:,:) = 0.0_WP
        end do
     end if
  case('-ym')
     if (yper.eq.0 .and. jproc.eq.1 .and. icyl.eq.0) then
        do j=jmino,jmin-1
           vec(:,j,:) = 0.0_WP
        end do
     end if
  case('+ym')
     if (yper.eq.0 .and. jproc.eq.npy) then
        do j=jmax+1,jmaxo
           vec(:,j,:) = 0.0_WP
        end do
     end if
  end select
  
  return
end subroutine boundary_dirichlet


! ============================ !
! Apply Von Neumann Conditions !
! ============================ !
subroutine boundary_neumann(vec,dir)
  use boundary
  use parallel
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec
  character(len=*), intent(in) :: dir
  integer :: i,j,k
  
  select case(trim(dir))
  case('-x')
     if (xper.eq.0 .and. iproc.eq.1) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmax_
              do i=imino,imin
                 vec(i,j,k) = vec(imin+1,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('+x')
     if (xper.eq.0 .and. iproc.eq.npx) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imax+1,imaxo
                 vec(i,j,k) = vec(imax,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('-y')
     if (yper.eq.0 .and. jproc.eq.1 .and. icyl.eq.0) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do i=imino_,imaxo_
              do j=jmino,jmin
                 vec(i,j,k) = vec(i,jmin+1,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('+y')
     if (yper.eq.0 .and. jproc.eq.npy) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do i=imino_,imaxo_
              do j=jmax+1,jmaxo
                 vec(i,j,k) = vec(i,jmax,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('-xm')
     if (xper.eq.0 .and. iproc.eq.1) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino,imin-1
                 vec(i,j,k) = vec(imin,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('+xm')
     if (xper.eq.0 .and. iproc.eq.npx) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imax+1,imaxo
                 vec(i,j,k) = vec(imax,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('-ym')
     if (yper.eq.0 .and. jproc.eq.1 .and. icyl.eq.0) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do i=imino_,imaxo_
              do j=jmino,jmin-1
                 vec(i,j,k) = vec(i,jmin,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case('+ym')
     if (yper.eq.0 .and. jproc.eq.npy) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do i=imino_,imaxo_
              do j=jmax+1,jmaxo
                 vec(i,j,k) = vec(i,jmax,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  end select

  return
end subroutine boundary_neumann


! ======================================== !
! Update the ghost cells                   !
! Not corresponding to boundary conditions !
! -> Domain decomposition                  !
! -> Periodic directions                   !
! -> Centerline                            !
! ======================================== !
subroutine boundary_update_border(A,sym,axis)
  use boundary
  use parallel
  implicit none
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  character(len=*), intent(in) :: sym,axis
  real(WP) :: coeff

  ! Domain decomposition
  ! Periodic direction
  call communication_border(A)
  
  ! Centerline
  if (trim(sym).eq.'+') coeff = +1.0_WP
  if (trim(sym).eq.'-') coeff = -1.0_WP
  if (trim(axis).eq.'y')  call centerline_update_y (A,coeff)
  if (trim(axis).eq.'ym') call centerline_update_ym(A,coeff)

  return
end subroutine boundary_update_border

! ==================================== !
! Boundary update border for intergers !
! ==================================== !
subroutine boundary_update_border_int(A,sym,axis)
  use boundary
  use parallel
  implicit none
  integer, dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  character(len=*), intent(in) :: sym,axis
  integer :: coeff

  ! Domain decomposition
  ! Periodic direction
  call communication_border_int(A)
  
  ! Centerline
  if (trim(sym).eq.'+') coeff = +1
  if (trim(sym).eq.'-') coeff = -1
  if (trim(axis).eq.'y')  call centerline_update_y_int (A,coeff)
  if (trim(axis).eq.'ym') call centerline_update_ym_int(A,coeff)

  return
end subroutine boundary_update_border_int


! =========================================== !
! Compute Dirichlet BCs for the new time step !
! And apply them to the velocity components   !
! -> Inlet                                    !
! =========================================== !
subroutine boundary_velocity_dirichlet
  use boundary
  use parallel
  use data
  implicit none

  ! Inflow Condition (left)
  call inflow_velocity
  
  ! Open Boundary Conditions (OBC)
  ! Bottom
  call boundary_dirichlet(V,'-y')
  ! Top
  call boundary_dirichlet(V,'+y')
  
  return
end subroutine boundary_velocity_dirichlet


! =========================================== !
! Compute Dirichlet BCs for the new time step !
! And apply them to the scalars               !
! -> Inlet                                    !
! =========================================== !
subroutine boundary_scalar_dirichlet
  use boundary
  use parallel
  use data
  implicit none

  ! Inflow Condition (left)
  call inflow_scalar
  
  return
end subroutine boundary_scalar_dirichlet


! ======================================= !
! Compute Neumann BCs for each iteration  !
! And apply them to the momentum          !
! -> Open Boundaries                      !
! ======================================= !
subroutine boundary_momentum_neumann
  use boundary
  use parallel
  use data
  implicit none
  
  ! Open Boundary Conditions (OBC)
  ! Bottom
  call boundary_neumann(rhoU,'-ym')
  call boundary_neumann(rhoW,'-ym')
  ! Top
  call boundary_neumann(rhoU,'+ym')
  call boundary_neumann(rhoW,'+ym')

!!$  call boundary_neumann(rhoV,'-y')
!!$  call boundary_neumann(rhoV,'+y')
  
  return
end subroutine boundary_momentum_neumann


! ======================================= !
! Compute Neumann BCs for each iteration  !
! And apply them to the velocity          !
! -> Open Boundaries                      !
! ======================================= !
subroutine boundary_velocity_neumann
  use boundary
  use parallel
  use data
  implicit none
  
  ! Open Boundary Conditions (OBC)
  ! Bottom
  call boundary_neumann(U,'-ym')
  call boundary_neumann(W,'-ym')
  ! Top
  call boundary_neumann(U,'+ym')
  call boundary_neumann(W,'+ym')

!!$  call boundary_neumann(V,'-y')
!!$  call boundary_neumann(V,'+y')
  
  return
end subroutine boundary_velocity_neumann


! =========================================== !
! Compute Neumann BCs for each iteration      !
! And apply them to the scalar field          !
! -> Open Boundaries                          !
! =========================================== !
subroutine boundary_scalar_neumann
  use boundary
  use parallel
  use data
  implicit none
  integer :: isc
  
  ! Open Boundary Conditions (OBC)
  ! Bottom
  do isc=1,nscalar
     call boundary_neumann(SC(:,:,:,isc),'-ym')
  end do
  ! Top
  do isc=1,nscalar
     call boundary_neumann(SC(:,:,:,isc),'+ym')
  end do
  
  return
end subroutine boundary_scalar_neumann


! =========================================== !
! Compute Neumann BCs for each iteration      !
! And apply them to the density field         !
! -> Open Boundaries                          !
! =========================================== !
subroutine boundary_density_neumann
  use boundary
  use parallel
  use data
  implicit none
  
  ! Open Boundary Conditions (OBC)
  ! Bottom
  call boundary_neumann(RHO,'-ym')
  ! Top
  call boundary_neumann(RHO,'+ym')
  
  return
end subroutine boundary_density_neumann


! ========================================= !
! Compute Convective BCs for each iteration !
! And apply them to the velocity            !
! -> Outflow                                !
! ========================================= !
subroutine boundary_velocity_outflow
  use boundary
  use parallel
  implicit none
  
  ! Outflow Conditions (right)
  call outflow_velocity
  
  return
end subroutine boundary_velocity_outflow


! ========================================== !
! Compute Convective BCs for each iteration  !
! And apply them to the velocity             !
! -> Outflow                                 !
! ========================================== !
subroutine boundary_scalar_outflow
  use boundary
  use parallel
  implicit none
  
  ! Outflow Conditions (right)
  call outflow_scalar
  
  return
end subroutine boundary_scalar_outflow


! =========================================== !
! Compute and correct mass fluxes             !
! Ensures Poisson solver convergence          !
! =========================================== !
subroutine boundary_velocity_massflux
  use boundary
  implicit none

  ! Compute Total 
  !   - mass flux across boundaries
  !   - change of mass inside the domain
  ! Correct outflow for mass consistency
  !call outflow_clip_negative
  call boundary_massflux
  call boundary_masschange
  call boundary_massadded
  masscorrection = massflux-masschange+massadded
  call outflow_correction
  
  return
end subroutine boundary_velocity_massflux

