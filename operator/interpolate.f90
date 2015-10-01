! ====================== !
! Interpolation routines !
! ====================== !
module interpolate
  use geometry
  use partition
  implicit none
  
  ! Cell centered velocities
  real(WP), dimension(:,:,:), pointer :: Ui
  real(WP), dimension(:,:,:), pointer :: Vi
  real(WP), dimension(:,:,:), pointer :: Wi

  ! Cell centered velocities -- Scalar
  real(WP), dimension(:,:,:), pointer :: Uti
  real(WP), dimension(:,:,:), pointer :: Vti
  real(WP), dimension(:,:,:), pointer :: Wti
  
end module interpolate


subroutine interpolate_init
  use interpolate
  implicit none

  integer :: i,j,k
  
  allocate(Ui(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Vi(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Wi(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           Ui(i,j,k)=0.0_WP
           Vi(i,j,k)=0.0_WP
           Wi(i,j,k)=0.0_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  allocate(Uti(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Vti(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Wti(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           Uti(i,j,k)=0.0_WP
           Vti(i,j,k)=0.0_WP
           Wti(i,j,k)=0.0_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine interpolate_init


subroutine interpolate_velocities
  use data
  use interpolate
  use metric_generic
  implicit none

  integer :: i,j,k
  
  ! Interpolate to cell centers
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           Ui(i,j,k) = sum(interp_u_xm(i,j,:) * U(i-st1:i+st2,j,k))
           Vi(i,j,k) = sum(interp_v_ym(i,j,:) * V(i,j-st1:j+st2,k))
           Wi(i,j,k) = sum(interp_w_zm(i,j,:) * W(i,j,k-st1:k+st2))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Update ghost cells
  call boundary_update_border(Ui,'+','ym')
  call boundary_update_border(Vi,'-','ym')
  call boundary_update_border(Wi,'-','ym')
  
  ! Update geometric ghost cells
  call boundary_neumann(Ui,'+ym')
  call boundary_neumann(Ui,'-ym')
  call boundary_neumann(Ui,'+xm')
  call boundary_neumann(Ui,'-xm')
  call boundary_neumann(Vi,'+ym')
  call boundary_neumann(Vi,'-ym')
  call boundary_neumann(Vi,'+xm')
  call boundary_neumann(Vi,'-xm')
  call boundary_neumann(Wi,'+ym')
  call boundary_neumann(Wi,'-ym')
  call boundary_neumann(Wi,'+xm')
  call boundary_neumann(Wi,'-xm')
  
  return
end subroutine interpolate_velocities


subroutine interpolate_velocities_total(isc)
  use data
  use interpolate
  use metric_generic
  implicit none

  integer, intent(in) :: isc
  integer :: i,j,k
  
  ! Interpolate to cell centers
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           Uti(i,j,k) = sum(interp_u_xm(i,j,:) * Ut(i-st1:i+st2,j,k,isc))
           Vti(i,j,k) = sum(interp_v_ym(i,j,:) * Vt(i,j-st1:j+st2,k,isc))
           Wti(i,j,k) = sum(interp_w_zm(i,j,:) * Wt(i,j,k-st1:k+st2,isc))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Update ghost cells
  call boundary_update_border(Uti,'+','ym')
  call boundary_update_border(Vti,'-','ym')
  call boundary_update_border(Wti,'-','ym')
  
  ! Update geometric ghost cells
  call boundary_neumann(Uti,'+ym')
  call boundary_neumann(Uti,'-ym')
  call boundary_neumann(Uti,'+xm')
  call boundary_neumann(Uti,'-xm')
  call boundary_neumann(Vti,'+ym')
  call boundary_neumann(Vti,'-ym')
  call boundary_neumann(Vti,'+xm')
  call boundary_neumann(Vti,'-xm')
  call boundary_neumann(Wti,'+ym')
  call boundary_neumann(Wti,'-ym')
  call boundary_neumann(Wti,'+xm')
  call boundary_neumann(Wti,'-xm')
  
  return
end subroutine interpolate_velocities_total
