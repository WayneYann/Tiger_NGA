module gradient
  use config
  use partition
  use parallel
  use metric_generic
  implicit none
  ! Dummy module
end module gradient


! ================================== !
! Compute the cell centered gradient !
! Takes A as input                   !
! Return G=|grad(A)|^2               !
! ================================== !
subroutine gradient_squared(A,G)
  use gradient
  implicit none
  
  real(WP) :: G1,G2,G3
  integer :: i,j,k
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: G
  
  ! Compute the square of the gradient in the domain
  !$OMP PARALLEL DO PRIVATE(G1,G2,G3)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           G1 = sum(grad_xm(i,j,:)*A(i-stp:i+stp,j,k))
           G2 = sum(grad_ym(i,j,:)*A(i,j-stp:j+stp,k))
           G3 = sum(grad_zm(i,j,:)*A(i,j,k-stp:k+stp))
           
           G(i,j,k) = G1**2 + G2**2 + G3**2

        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Update ghost cells
  call boundary_update_border(G,'+','ym')
  
  ! Non periodic directions
  call boundary_neumann(G,'-xm')
  call boundary_neumann(G,'+xm')
  call boundary_neumann(G,'-ym')
  call boundary_neumann(G,'+ym')
  
  return
end subroutine gradient_squared


! ================================== !
! Compute the cell centered gradient !
! Takes A as input                   !
! Return (G1,G2,G3)=grad(A)          !
! ================================== !
subroutine gradient_vector(A,G1,G2,G3)
  use gradient
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: G1,G2,G3
  integer :: i,j,k
  
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           G1(i,j,k) = sum(grad_xm(i,j,:)*A(i-stp:i+stp,j,k))
           G2(i,j,k) = sum(grad_ym(i,j,:)*A(i,j-stp:j+stp,k))
           G3(i,j,k) = sum(grad_zm(i,j,:)*A(i,j,k-stp:k+stp))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Update ghost cells
  call boundary_update_border(G1,'+','ym')
  call boundary_update_border(G2,'-','ym')
  call boundary_update_border(G3,'-','ym')
  
  ! Non periodic directions
  call boundary_neumann(G1,'-xm')
  call boundary_neumann(G2,'-xm')
  call boundary_neumann(G3,'-xm')

  call boundary_neumann(G1,'+xm')
  call boundary_neumann(G2,'+xm')
  call boundary_neumann(G3,'+xm')

  call boundary_neumann(G1,'-ym')
  call boundary_neumann(G2,'-ym')
  call boundary_neumann(G3,'-ym')

  call boundary_neumann(G1,'+ym')
  call boundary_neumann(G2,'+ym')
  call boundary_neumann(G3,'+ym')
  
  return
end subroutine gradient_vector

