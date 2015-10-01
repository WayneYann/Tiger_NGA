module strainrate
  use config
  use partition
  use parallel
  use metric_generic
  implicit none
  
  ! Store the norm of the strain rate
  real(WP), dimension(:,:,:), pointer :: S

  ! Store the norm of the total strain rate
  real(WP), dimension(:,:,:), pointer :: St
  
  ! P,Q,R invariants
  real(WP), dimension(:,:,:), pointer :: Pcrit
  real(WP), dimension(:,:,:), pointer :: Qcrit
  real(WP), dimension(:,:,:), pointer :: Rcrit
  
contains
  
  ! Compute velocity gradient tensor
  subroutine vel_grad_local(i,j,k,dUdx)
    use interpolate
    use data
    implicit none
    
    real(WP), dimension(3,3), intent(out) :: dUdx
    integer, intent(in) :: i,j,k
    
    ! dU/dx
    dUdx(1,1)=sum(div_u(i,j,:)*U(i-st1:i+st2,j,k))
    ! dV/dx
    dUdx(2,1)=dxi(i)*( +sum(interp_uvw_x(i+1,j,:)*Vi(i-st2+1:i+st1+1,j,k)) &
                       -sum(interp_uvw_x(i,j,:)  *Vi(i-st2:i+st1,j,k))    )
    ! dW/dx
    dUdx(3,1)=dxi(i)*( +sum(interp_uvw_x(i+1,j,:)*Wi(i-st2+1:i+st1+1,j,k)) &
                       -sum(interp_uvw_x(i,j,:)  *Wi(i-st2:i+st1,j,k))    )
    ! dU/dy
    dUdx(1,2)=dyi(j)*( +sum(interp_uvw_y(i,j+1,:)*Ui(i,j-st2+1:j+st1+1,k)) &
                       -sum(interp_uvw_y(i,j,:)  *Ui(i,j-st2:j+st1,k))    )
    ! dV/dy
    dUdx(2,2)=sum(div_v(i,j,:)*V(i,j-st1:j+st2,k)) - icyl*ymi(j)*Vi(i,j,k)
    ! dW/dy
    dUdx(3,2)=dyi(j)*( +sum(interp_uvw_y(i,j+1,:)*Wi(i,j-st2+1:j+st1+1,k)) &
                       -sum(interp_uvw_y(i,j,:)  *Wi(i,j-st2:j+st1,k))    )
    ! dU/dz
    dUdx(1,3)=dzi_u(j)*( +sum(interp_uvw_z(i,j,:)*Ui(i,j,k-st2+1:k+st1+1)) &
                         -sum(interp_uvw_z(i,j,:)*Ui(i,j,k-st2:k+st1))    )
    ! dV/dz
    dUdx(2,3)=dzi_v(j)*( +sum(interp_uvw_z(i,j,:)*Vi(i,j,k-st2+1:k+st1+1)) &
                         -sum(interp_uvw_z(i,j,:)*Vi(i,j,k-st2:k+st1))    )&
                         -icyl * ymi(j) * Wi(i,j,k)
    ! dW/dz
    dUdx(3,3)=sum(div_w(i,j,:)*W(i,j,k-st1:k+st2)) + icyl*ymi(j)*Vi(i,j,k)
           
    return
  end subroutine vel_grad_local

  ! Compute total velocity gradient tensor
  subroutine vel_tot_grad_local(isc,i,j,k,dUtdx)
    use interpolate
    use data
    implicit none
    
    integer, intent(in) :: isc
    real(WP), dimension(3,3), intent(out) :: dUtdx
    integer, intent(in) :: i,j,k
    
    ! dUt/dx
    dUtdx(1,1)=sum(div_u(i,j,:)*Ut(i-st1:i+st2,j,k,isc))
    ! dVt/dx
    dUtdx(2,1)=dxi(i)*( +sum(interp_uvw_x(i+1,j,:)*Vti(i-st2+1:i+st1+1,j,k)) &
                       -sum(interp_uvw_x(i,j,:)  *Vti(i-st2:i+st1,j,k))    )
    ! dWt/dx
    dUtdx(3,1)=dxi(i)*( +sum(interp_uvw_x(i+1,j,:)*Wti(i-st2+1:i+st1+1,j,k)) &
                       -sum(interp_uvw_x(i,j,:)  *Wti(i-st2:i+st1,j,k))    )
    ! dUt/dy
    dUtdx(1,2)=dyi(j)*( +sum(interp_uvw_y(i,j+1,:)*Uti(i,j-st2+1:j+st1+1,k)) &
                       -sum(interp_uvw_y(i,j,:)  *Uti(i,j-st2:j+st1,k))    )
    ! dVt/dy
    dUtdx(2,2)=sum(div_v(i,j,:)*Vt(i,j-st1:j+st2,k,isc)) - icyl*ymi(j)*Vti(i,j,k)
    ! dWt/dy
    dUtdx(3,2)=dyi(j)*( +sum(interp_uvw_y(i,j+1,:)*Wti(i,j-st2+1:j+st1+1,k)) &
                       -sum(interp_uvw_y(i,j,:)  *Wti(i,j-st2:j+st1,k))    )
    ! dUt/dz
    dUtdx(1,3)=dzi_u(j)*( +sum(interp_uvw_z(i,j,:)*Uti(i,j,k-st2+1:k+st1+1)) &
                         -sum(interp_uvw_z(i,j,:)*Uti(i,j,k-st2:k+st1))    )
    ! dVt/dz
    dUtdx(2,3)=dzi_v(j)*( +sum(interp_uvw_z(i,j,:)*Vti(i,j,k-st2+1:k+st1+1)) &
                         -sum(interp_uvw_z(i,j,:)*Vti(i,j,k-st2:k+st1))    )&
                         -icyl * ymi(j) * Wti(i,j,k)
    ! dWt/dz
    dUtdx(3,3)=sum(div_w(i,j,:)*Wt(i,j,k-st1:k+st2,isc)) + icyl*ymi(j)*Vti(i,j,k)
           
    return
  end subroutine vel_tot_grad_local
  
end module strainrate


subroutine strainrate_init
  use strainrate
  implicit none
  
  allocate(S(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(St(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Pcrit(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Qcrit(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Rcrit(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  return
end subroutine strainrate_init


! Compute the strain rate
! Take U,V,W,Ui,Vi,Wi from ...
! Return the 6 components of the strain rate and its norm
subroutine strainrate_compute(Sij1,Sij2,Sij3,Sij4,Sij5,Sij6)
  use strainrate
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: Sij1,Sij2,Sij3,Sij4,Sij5,Sij6
  real(WP), dimension(3,3) :: dUdx
  integer :: i,j,k
  real(WP) :: div
  
  !     ( 1 4 6 )
  ! S = ( 4 2 5 )
  !     ( 6 5 3 )
  
  ! Enforce zero in the geometric ghost cells
  Sij1 = 0.0_WP
  Sij2 = 0.0_WP
  Sij3 = 0.0_WP
  Sij4 = 0.0_WP
  Sij5 = 0.0_WP
  Sij6 = 0.0_WP
  
  !$OMP PARALLEL DO PRIVATE(dUdx,div)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           call vel_grad_local(i,j,k,dUdx)
           div = dUdx(1,1) + dUdx(2,2) + dUdx(3,3)
           
           ! u-u: du/dx - 1/3*(du_k/dx_k)
           Sij1(i,j,k) = dUdx(1,1) - div/3.0_WP
           
           ! v-v: dv/dy - 1/3*(du_k/dx_k)
           Sij2(i,j,k) = dUdx(2,2) - div/3.0_WP
           
           ! w-w: dw/dz - 1/3*(du_k/dx_k)
           Sij3(i,j,k) = dUdx(3,3) - div/3.0_WP
           
           ! u-v: 1/2*(du/dy+dv/dx)
           Sij4(i,j,k) = 0.5_WP*(dUdx(1,2)+dUdx(2,1))
           
           ! v-w: 1/2*(dw/dy+dv/dz)
           Sij5(i,j,k) = 0.5_WP*(dUdx(2,3)+dUdx(3,2))
           
           ! u-w: 1/2*(du/dz+dw/dx)
           Sij6(i,j,k) = 0.5_WP*(dUdx(1,3)+dUdx(3,1))
           
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Non periodic directions
  call boundary_neumann(Sij1,'-xm')
  call boundary_neumann(Sij2,'-xm')
  call boundary_neumann(Sij3,'-xm')
  call boundary_neumann(Sij4,'-xm')
  call boundary_neumann(Sij5,'-xm')
  call boundary_neumann(Sij6,'-xm')

  call boundary_neumann(Sij1,'+xm')
  call boundary_neumann(Sij2,'+xm')
  call boundary_neumann(Sij3,'+xm')
  call boundary_neumann(Sij4,'+xm')
  call boundary_neumann(Sij5,'+xm')
  call boundary_neumann(Sij6,'+xm')

  call boundary_neumann(Sij1,'-ym')
  call boundary_neumann(Sij2,'-ym')
  call boundary_neumann(Sij3,'-ym')
  call boundary_neumann(Sij4,'-ym')
  call boundary_neumann(Sij5,'-ym')
  call boundary_neumann(Sij6,'-ym')

  call boundary_neumann(Sij1,'+ym')
  call boundary_neumann(Sij2,'+ym')
  call boundary_neumann(Sij3,'+ym')
  call boundary_neumann(Sij4,'+ym')
  call boundary_neumann(Sij5,'+ym')
  call boundary_neumann(Sij6,'+ym')

  ! Update ghost cells
  call boundary_update_border(Sij1,'+','ym')
  call boundary_update_border(Sij2,'+','ym')
  call boundary_update_border(Sij3,'+','ym')
  call boundary_update_border(Sij4,'-','ym')
  call boundary_update_border(Sij5,'-','ym')
  call boundary_update_border(Sij6,'+','ym')
  
  ! Compute the norm of the strain rate
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           S(i,j,k) = sqrt( Sij1(i,j,k)**2 + Sij2(i,j,k)**2 + Sij3(i,j,k)**2 + &
                2.0_WP*( Sij4(i,j,k)**2 + Sij5(i,j,k)**2 + Sij6(i,j,k)**2 ))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine strainrate_compute


! Compute the total strain rate for the scalar
! Take Ut,Vt,Wt,Uti,Vti,Wti from ...
! Return the 6 components of the strain rate and its norm
subroutine strainrate_total_compute(Stij1,Stij2,Stij3,Stij4,Stij5,Stij6,isc)
  use strainrate
  implicit none
  
  integer, intent(in) :: isc
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: Stij1,Stij2,Stij3,Stij4,Stij5,Stij6
  real(WP), dimension(3,3) :: dUtdx
  integer :: i,j,k
  real(WP) :: divt
  
  !      ( 1 4 6 )
  ! St = ( 4 2 5 )
  !      ( 6 5 3 )
  
  ! Enforce zero in the geometric ghost cells
  Stij1 = 0.0_WP
  Stij2 = 0.0_WP
  Stij3 = 0.0_WP
  Stij4 = 0.0_WP
  Stij5 = 0.0_WP
  Stij6 = 0.0_WP
  
  !$OMP PARALLEL DO PRIVATE(dUtdx,divt)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           call vel_tot_grad_local(isc,i,j,k,dUtdx)
           divt = dUtdx(1,1) + dUtdx(2,2) + dUtdx(3,3)
           
           ! u-u: du/dx - 1/3*(du_k/dx_k)
           Stij1(i,j,k) = dUtdx(1,1) - divt/3.0_WP
           
           ! v-v: dv/dy - 1/3*(du_k/dx_k)
           Stij2(i,j,k) = dUtdx(2,2) - divt/3.0_WP
           
           ! w-w: dw/dz - 1/3*(du_k/dx_k)
           Stij3(i,j,k) = dUtdx(3,3) - divt/3.0_WP
           
           ! u-v: 1/2*(du/dy+dv/dx)
           Stij4(i,j,k) = 0.5_WP*(dUtdx(1,2)+dUtdx(2,1))
           
           ! v-w: 1/2*(dw/dy+dv/dz)
           Stij5(i,j,k) = 0.5_WP*(dUtdx(2,3)+dUtdx(3,2))
           
           ! u-w: 1/2*(du/dz+dw/dx)
           Stij6(i,j,k) = 0.5_WP*(dUtdx(1,3)+dUtdx(3,1))
           
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Non periodic directions
  call boundary_neumann(Stij1,'-xm')
  call boundary_neumann(Stij2,'-xm')
  call boundary_neumann(Stij3,'-xm')
  call boundary_neumann(Stij4,'-xm')
  call boundary_neumann(Stij5,'-xm')
  call boundary_neumann(Stij6,'-xm')

  call boundary_neumann(Stij1,'+xm')
  call boundary_neumann(Stij2,'+xm')
  call boundary_neumann(Stij3,'+xm')
  call boundary_neumann(Stij4,'+xm')
  call boundary_neumann(Stij5,'+xm')
  call boundary_neumann(Stij6,'+xm')

  call boundary_neumann(Stij1,'-ym')
  call boundary_neumann(Stij2,'-ym')
  call boundary_neumann(Stij3,'-ym')
  call boundary_neumann(Stij4,'-ym')
  call boundary_neumann(Stij5,'-ym')
  call boundary_neumann(Stij6,'-ym')

  call boundary_neumann(Stij1,'+ym')
  call boundary_neumann(Stij2,'+ym')
  call boundary_neumann(Stij3,'+ym')
  call boundary_neumann(Stij4,'+ym')
  call boundary_neumann(Stij5,'+ym')
  call boundary_neumann(Stij6,'+ym')

  ! Update ghost cells
  call boundary_update_border(Stij1,'+','ym')
  call boundary_update_border(Stij2,'+','ym')
  call boundary_update_border(Stij3,'+','ym')
  call boundary_update_border(Stij4,'-','ym')
  call boundary_update_border(Stij5,'-','ym')
  call boundary_update_border(Stij6,'+','ym')
  
  ! Compute the norm of the strain rate
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           St(i,j,k) = sqrt( Stij1(i,j,k)**2 + Stij2(i,j,k)**2 + Stij3(i,j,k)**2 + &
                2.0_WP*( Stij4(i,j,k)**2 + Stij5(i,j,k)**2 + Stij6(i,j,k)**2 ))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine strainrate_total_compute


! ====================================================================== !
! Compute (P,Q,R) - velocity gradient tensor invariants                  !
! M.S. Chong, A.E. Perry and B.J. Cantwell, A general classification     !
! of three-dimensional flow fields, PoF A 2 (5), May 1990                !
! Aij=dui/dxj, Sij=0.5(Aij+Aji), Oij=0.5(Aij-Aji)                        !
! P=-Sii, Q=0.5(P^2-SijSji-OijOji), R=1/3(-P^3+3PQ-SijSjkOki-3OijOjkSki) !
! ====================================================================== !
subroutine strainrate_PQRcrit
  use strainrate
  use data
  use interpolate
  implicit none
  
  integer  :: i,j,k,n1,n2,n3
  real(WP) :: buf1,buf2
  real(WP), dimension(3,3) :: dUdx
  real(WP), dimension(3,3) :: Oij
  real(WP), dimension(3,3) :: Sij
  
  !        ( 1 2 3 )
  ! dUdx = ( 4 5 6 )
  !        ( 7 8 9 )
  
  ! Enforce zero in the geometric ghost cells
  Pcrit = 0.0_WP
  Qcrit = 0.0_WP
  Rcrit = 0.0_WP
  
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           call vel_grad_local(i,j,k,dUdx)

           ! Construct Oij and Sij
           do n1=1,3
              do n2=1,3
                 Oij(n1,n2)=0.5_WP*(dUdx(n1,n2)-dUdx(n2,n1))
                 Sij(n1,n2)=0.5_WP*(dUdx(n1,n2)+dUdx(n2,n1))
              end do
           end do
           
           ! Compute Pcrit
           buf1=0.0_WP;buf2=0.0_WP
           do n1=1,3
              buf1=buf1-Sij(n1,n1)
           end do
           Pcrit(i,j,k)=buf1
           
           ! Compute Qcrit
           buf1=0.0_WP;buf2=0.0_WP
           do n1=1,3
              do n2=1,3
                 buf1=buf1+Sij(n1,n2)*Sij(n2,n1)
                 buf2=buf2+Oij(n1,n2)*Oij(n2,n1)
              end do
           end do
           Qcrit(i,j,k)=1.0_WP/2.0_WP*(Pcrit(i,j,k)**2-buf1-buf2)
           
           ! Compute Rcrit
           buf1=0.0_WP;buf2=0.0_WP
           do n1=1,3
              do n2=1,3
                 do n3=1,3
                    buf1=buf1+Sij(n1,n2)*Sij(n2,n3)*Sij(n3,n1)
                    buf2=buf2+Oij(n1,n2)*Oij(n2,n3)*Sij(n3,n1)
                 end do
              end do
           end do
           Rcrit(i,j,k)=1.0_WP/3.0_WP*(-Pcrit(i,j,k)**3+3.0_WP*Pcrit(i,j,k)*Qcrit(i,j,k)-buf1-3.0_WP*buf2)
           
        end do
     end do
  end do
  
  ! Update BCs
  call boundary_dirichlet(Pcrit,'-xm')
  call boundary_neumann  (Pcrit,'+xm')
  call boundary_neumann  (Pcrit,'-ym')
  call boundary_neumann  (Pcrit,'+ym')
  call boundary_update_border(Pcrit,'+','ym')
  
  call boundary_dirichlet(Qcrit,'-xm')
  call boundary_neumann  (Qcrit,'+xm')
  call boundary_neumann  (Qcrit,'-ym')
  call boundary_neumann  (Qcrit,'+ym')
  call boundary_update_border(Qcrit,'+','ym')
  
  call boundary_dirichlet(Rcrit,'-xm')
  call boundary_neumann  (Rcrit,'+xm')
  call boundary_neumann  (Rcrit,'-ym')
  call boundary_neumann  (Rcrit,'+ym')
  call boundary_update_border(Rcrit,'+','ym')
  
  return
end subroutine strainrate_PQRcrit
