module sgs_gradient
  use sgsmodel
  implicit none
  
  ! Size of the mesh in the z directioin
  real(WP), dimension(:), pointer :: dz_c
  
  !                         ( 11 12 13 )
  ! Velocity Gradient : G = ( 21 22 23 )
  !                         ( 31 32 33 )
  real(WP), dimension(:,:,:), pointer :: Gij11,Gij12,Gij13
  real(WP), dimension(:,:,:), pointer :: Gij21,Gij22,Gij23
  real(WP), dimension(:,:,:), pointer :: Gij31,Gij32,Gij33

  ! Scalar Gradient : G = ( 1 2 3 )
  real(WP), dimension(:,:,:), pointer :: GSi1,GSi2,GSi3
  
contains
  
  ! Compute the gradient tensor of velocities
  ! -----------------------------------------
  subroutine compute_gradient_velocity
    use strainrate
    use filter
    implicit none
    
    real(WP), dimension(3,3) :: dUdx
    integer :: i,j,k
    
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             
             call vel_grad_local(i,j,k,dUdx)
             dUdx(:,1) = dUdx(:,1)*dx(i)
             dUdx(:,2) = dUdx(:,2)*dy(j)
             dUdx(:,3) = dUdx(:,3)*dz_c(j)
             
             Hij1(i,j,k) = -RHO(i,j,k)*sum(dUdx(1,:)*dUdx(1,:))/12.0_WP
             Hij2(i,j,k) = -RHO(i,j,k)*sum(dUdx(2,:)*dUdx(2,:))/12.0_WP
             Hij3(i,j,k) = -RHO(i,j,k)*sum(dUdx(3,:)*dUdx(3,:))/12.0_WP
             Hij4(i,j,k) = -RHO(i,j,k)*sum(dUdx(1,:)*dUdx(2,:))/12.0_WP
             Hij5(i,j,k) = -RHO(i,j,k)*sum(dUdx(2,:)*dUdx(3,:))/12.0_WP
             Hij6(i,j,k) = -RHO(i,j,k)*sum(dUdx(1,:)*dUdx(3,:))/12.0_WP
             
          end do
       end do
    end do
    
    ! Update Ghost Cells
    call boundary_update_border(Hij1,'+','n')
    call boundary_update_border(Hij2,'+','n')
    call boundary_update_border(Hij3,'+','n')
    call boundary_update_border(Hij4,'-','n')
    call boundary_update_border(Hij5,'+','n')
    call boundary_update_border(Hij6,'-','n')
    
    return
  end subroutine compute_gradient_velocity
  
  ! Compute the gradient vector of scalars
  ! -----------------------------------------
  subroutine compute_gradient_scalar(isc)
    use strainrate
    use filter
    implicit none
    
    integer, intent(in) :: isc
    real(WP), dimension(3,3) :: dUtdx
    real(WP), dimension(3) :: dZdx
    integer :: i,j,k
    
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             
             call vel_tot_grad_local(isc,i,j,k,dUtdx)
             dUtdx(:,1) = dUtdx(:,1)*dx(i)
             dUtdx(:,2) = dUtdx(:,2)*dy(j)
             dUtdx(:,3) = dUtdx(:,3)*dz_c(j)
             dZdx = sum(grad_xm(i,j,:)*SC(i-stp:i+stp,j,k,isc))*dx(i)
             dZdx = sum(grad_ym(i,j,:)*SC(i,j-stp:j+stp,k,isc))*dy(j)
             dZdx = sum(grad_zm(i,j,:)*SC(i,j,k-stp:k+stp,isc))*dz_c(j)

             Hi1(i,j,k) = -RHO(i,j,k)*sum(dZdx*dUtdx(1,:))/12.0_WP
             Hi2(i,j,k) = -RHO(i,j,k)*sum(dZdx*dUtdx(2,:))/12.0_WP
             Hi3(i,j,k) = -RHO(i,j,k)*sum(dZdx*dUtdx(3,:))/12.0_WP
             
          end do
       end do
    end do
    
    ! Update Ghost Cells
    call boundary_update_border(Hi1,'+','n')
    call boundary_update_border(Hi2,'-','n')
    call boundary_update_border(Hi3,'-','n')
    
    return
  end subroutine compute_gradient_scalar
  
end module sgs_gradient


! ======================== !
! Initialize the SGS model !
! ======================== !
subroutine sgs_gradient_init
  use sgs_gradient
  use memory
  implicit none
  
  ! Link temporary arrays
  Gij11=>tmp1;  Gij12=>tmp2;  Gij13=>tmp3
  Gij21=>tmp4;  Gij22=>tmp5;  Gij23=>tmp6
  Gij31=>tmp7;  Gij32=>tmp8;  Gij33=>tmp9
  GSi1 =>tmp13; GSi2 =>tmp14; GSi3 =>tmp15
  
  ! Additionnal notations
  allocate(dz_c(jmin:jmax)) 
  if (icyl .eq. 1) then
     dz_c = y(jmin:jmax)*dz
  else
     dz_c = dz
  end if
  
  return
end subroutine sgs_gradient_init


! ========================================= !
! Compute the gradient tensor of velocities !
! ========================================= !
subroutine sgs_gradient_compute_Hij
  use sgs_gradient
  use strainrate
  use filter
  implicit none
  
  real(WP), dimension(3,3) :: dUdx
  integer :: i,j,k
  
  ! Compute the Model Term
  call compute_gradient_velocity
  
  ! Filter the Model Term
  call filter_global_3D(Hij1,Hij1,'+','n')
  call filter_global_3D(Hij2,Hij2,'+','n')
  call filter_global_3D(Hij3,Hij3,'+','n')
  call filter_global_3D(Hij4,Hij4,'-','n')
  call filter_global_3D(Hij5,Hij5,'+','n')
  call filter_global_3D(Hij6,Hij6,'-','n')
  
  ! Compute the Model Term based on the Filtered Velocities
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_

           call vel_grad_local(i,j,k,dUdx)
           dUdx(:,1) = dUdx(:,1)*dx(i)
           dUdx(:,2) = dUdx(:,2)*dy(j)
           dUdx(:,3) = dUdx(:,3)*dz_c(j)
           
           Gij11(i,j,k) = RHO(i,j,k)*dUdx(1,1)
           Gij12(i,j,k) = RHO(i,j,k)*dUdx(1,2)
           Gij13(i,j,k) = RHO(i,j,k)*dUdx(1,3)
           Gij21(i,j,k) = RHO(i,j,k)*dUdx(2,1)
           Gij22(i,j,k) = RHO(i,j,k)*dUdx(2,2)
           Gij23(i,j,k) = RHO(i,j,k)*dUdx(2,3)
           Gij31(i,j,k) = RHO(i,j,k)*dUdx(3,1)
           Gij32(i,j,k) = RHO(i,j,k)*dUdx(3,2)
           Gij33(i,j,k) = RHO(i,j,k)*dUdx(3,3)
        end do
     end do
  end do
  
  ! Update Ghost Cells
  call boundary_update_border(Gij11,'+','n')
  call boundary_update_border(Gij12,'-','n')
  call boundary_update_border(Gij13,'-','n')
  call boundary_update_border(Gij21,'-','n')
  call boundary_update_border(Gij22,'+','n')
  call boundary_update_border(Gij23,'+','n')
  call boundary_update_border(Gij31,'-','n')
  call boundary_update_border(Gij32,'+','n')
  call boundary_update_border(Gij33,'+','n')
  
  ! Apply filter to the velocities and density
  call filter_global_3D(Gij11,Gij11,'+','n')
  call filter_global_3D(Gij12,Gij12,'-','n')
  call filter_global_3D(Gij13,Gij13,'-','n')
  call filter_global_3D(Gij21,Gij21,'-','n')
  call filter_global_3D(Gij22,Gij22,'+','n')
  call filter_global_3D(Gij23,Gij23,'+','n')
  call filter_global_3D(Gij31,Gij31,'-','n')
  call filter_global_3D(Gij32,Gij32,'+','n')
  call filter_global_3D(Gij33,Gij33,'+','n')
  
  ! Compute the full term
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           Hij1(i,j,k) = - Hij1(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (Gij11(i,j,k)*Gij11(i,j,k)+Gij12(i,j,k)*Gij12(i,j,k)+Gij13(i,j,k)*Gij13(i,j,k))/12.0_WP
           Hij2(i,j,k) = - Hij2(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (Gij21(i,j,k)*Gij21(i,j,k)+Gij22(i,j,k)*Gij22(i,j,k)+Gij23(i,j,k)*Gij23(i,j,k))/12.0_WP
           Hij3(i,j,k) = - Hij3(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (Gij31(i,j,k)*Gij31(i,j,k)+Gij32(i,j,k)*Gij32(i,j,k)+Gij33(i,j,k)*Gij33(i,j,k))/12.0_WP
           Hij4(i,j,k) = - Hij4(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (Gij11(i,j,k)*Gij21(i,j,k)+Gij12(i,j,k)*Gij22(i,j,k)+Gij13(i,j,k)*Gij23(i,j,k))/12.0_WP
           Hij5(i,j,k) = - Hij5(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (Gij21(i,j,k)*Gij31(i,j,k)+Gij22(i,j,k)*Gij32(i,j,k)+Gij23(i,j,k)*Gij33(i,j,k))/12.0_WP
           Hij6(i,j,k) = - Hij6(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (Gij11(i,j,k)*Gij31(i,j,k)+Gij12(i,j,k)*Gij32(i,j,k)+Gij13(i,j,k)*Gij33(i,j,k))/12.0_WP
        end do
     end do
  end do
  
  ! Update Ghost Cells
  call boundary_update_border(Hij1,'+','n')
  call boundary_update_border(Hij2,'+','n')
  call boundary_update_border(Hij3,'+','n')
  call boundary_update_border(Hij4,'-','n')
  call boundary_update_border(Hij5,'+','n')
  call boundary_update_border(Hij6,'-','n')
  
  return
end subroutine sgs_gradient_compute_Hij


! ========================================= !
! Compute the gradient tensor of velocities !
! ========================================= !
subroutine sgs_gradient_compute_Hi(isc)
  use sgs_gradient
  use strainrate
  use filter
  implicit none
  
  integer, intent(in) :: isc
  real(WP), dimension(3,3) :: dUtdx
  real(WP), dimension(3) :: dZdx
  integer :: i,j,k
  
  ! Compute the Model Term
  call compute_gradient_scalar(isc)
  
  ! Filter the Model Term
  call filter_global_3D(Hi1,Hi1,'+','n')
  call filter_global_3D(Hi2,Hi2,'-','n')
  call filter_global_3D(Hi3,Hi3,'-','n')
  
  ! Compute the Model Term based on the Filtered Velocities
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_

           call vel_tot_grad_local(isc,i,j,k,dUtdx)
           dUtdx(:,1) = dUtdx(:,1)*dx(i)
           dUtdx(:,2) = dUtdx(:,2)*dy(j)
           dUtdx(:,3) = dUtdx(:,3)*dz_c(j)
           
           Gij11(i,j,k) = RHO(i,j,k)*dUtdx(1,1)
           Gij12(i,j,k) = RHO(i,j,k)*dUtdx(1,2)
           Gij13(i,j,k) = RHO(i,j,k)*dUtdx(1,3)
           Gij21(i,j,k) = RHO(i,j,k)*dUtdx(2,1)
           Gij22(i,j,k) = RHO(i,j,k)*dUtdx(2,2)
           Gij23(i,j,k) = RHO(i,j,k)*dUtdx(2,3)
           Gij31(i,j,k) = RHO(i,j,k)*dUtdx(3,1)
           Gij32(i,j,k) = RHO(i,j,k)*dUtdx(3,2)
           Gij33(i,j,k) = RHO(i,j,k)*dUtdx(3,3)

           dZdx = sum(grad_xm(i,j,:)*SC(i-stp:i+stp,j,k,isc))*dx(i)
           dZdx = sum(grad_ym(i,j,:)*SC(i,j-stp:j+stp,k,isc))*dy(j)
           dZdx = sum(grad_zm(i,j,:)*SC(i,j,k-stp:k+stp,isc))*dz_c(j)

           GSi1(i,j,k) = RHO(i,j,k)*dZdx(1)
           GSi2(i,j,k) = RHO(i,j,k)*dZdx(2)
           GSi3(i,j,k) = RHO(i,j,k)*dZdx(3)
        end do
     end do
  end do
  
  ! Update Ghost Cells
  call boundary_update_border(Gij11,'+','n')
  call boundary_update_border(Gij12,'-','n')
  call boundary_update_border(Gij13,'-','n')
  call boundary_update_border(Gij21,'-','n')
  call boundary_update_border(Gij22,'+','n')
  call boundary_update_border(Gij23,'+','n')
  call boundary_update_border(Gij31,'-','n')
  call boundary_update_border(Gij32,'+','n')
  call boundary_update_border(Gij33,'+','n')

  call boundary_update_border(GSi1,'+','n')
  call boundary_update_border(GSi2,'-','n')
  call boundary_update_border(GSi3,'-','n')
  
  ! Apply filter to the velocities and density
  call filter_global_3D(Gij11,Gij11,'+','n')
  call filter_global_3D(Gij12,Gij12,'-','n')
  call filter_global_3D(Gij13,Gij13,'-','n')
  call filter_global_3D(Gij21,Gij21,'-','n')
  call filter_global_3D(Gij22,Gij22,'+','n')
  call filter_global_3D(Gij23,Gij23,'+','n')
  call filter_global_3D(Gij31,Gij31,'-','n')
  call filter_global_3D(Gij32,Gij32,'+','n')
  call filter_global_3D(Gij33,Gij33,'+','n')
  
  call filter_global_3D(GSi1,GSi1,'+','n')
  call filter_global_3D(GSi2,GSi2,'-','n')
  call filter_global_3D(GSi3,GSi3,'-','n')

  ! Compute the full term
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           Hi1(i,j,k) = - Hi1(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (GSi1(i,j,k)*Gij11(i,j,k)+GSi2(i,j,k)*Gij12(i,j,k)+GSi3(i,j,k)*Gij13(i,j,k))/12.0_WP
           Hi2(i,j,k) = - Hi2(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (GSi1(i,j,k)*Gij21(i,j,k)+GSi2(i,j,k)*Gij22(i,j,k)+GSi3(i,j,k)*Gij23(i,j,k))/12.0_WP
           Hi3(i,j,k) = - Hi3(i,j,k) + ratio_3D(i,j)**2/FRHO(i,j,k)* &
                (GSi1(i,j,k)*Gij31(i,j,k)+GSi2(i,j,k)*Gij32(i,j,k)+GSi3(i,j,k)*Gij33(i,j,k))/12.0_WP
        end do
     end do
  end do
  
  ! Update Ghost Cells
  call boundary_update_border(Hi1,'+','n')
  call boundary_update_border(Hi2,'-','n')
  call boundary_update_border(Hi3,'-','n')
  
  return
end subroutine sgs_gradient_compute_Hi


! ======================= !
! Create SGS source terms !
! ======================= !
subroutine sgs_gradient_src_vel(srcU,srcV,srcW)
  use sgs_gradient
  use memory
  use filter
  use metric_velocity_visc
  use metric_generic
  use time_info
  use masks
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: srcU,srcV,srcW
  integer  :: i,j,k,ii,jj,kk

  ! Compute the tensor for the gradient model
  call compute_gradient_velocity
  
  ! ==== U ====
  ! Compute the fluxes
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           FX(i,j,k) = Hij1(i,j,k)

           i = ii; j = jj; k = kk;
           FY(i,j,k) = sum(interp_sc_xy(i,j,:,:)*Hij4(i-st2:i+st1,j-st2:j+st1,k))
           FZ(i,j,k) = sum(interp_sc_xz(i,j,:,:)*Hij6(i-st2:i+st1,j,k-st2:k+st1))
        end do
     end do
  end do
  
  ! Compute the source terms
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_           
           srcU(i,j,k) = srcU(i,j,k) + dt_uvw*( &
                +sum(divv_xx(i,j,:)*FX(i-stv2:i+stv1,j,k)) &
                +sum(divv_xy(i,j,:)*FY(i,j-stv1:j+stv2,k)) &
                +sum(divv_xz(i,j,:)*FZ(i,j,k-stv1:k+stv2)) )
        end do
     end do
  end do
  
  
  ! ==== V ====
  ! Compute the fluxes
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           FY(i,j,k)    = Hij2(i,j,k)
           Fcylv(i,j,k) = Hij3(i,j,k)
           
           i = ii; j = jj; k = kk;
           FX(i,j,k) = sum(interp_sc_xy(i,j,:,:)*Hij4(i-st2:i+st1,j-st2:j+st1,k))
           FZ(i,j,k) = sum(interp_sc_yz(i,j,:,:)*Hij5(i,j-st2:j+st1,k-st2:k+st1))
        end do
     end do
  end do
  
  ! Compute the source terms
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           srcV(i,j,k) = srcV(i,j,k) + dt_uvw*( &
                +sum(divv_yx(i,j,:)*FX(i-stv1:i+stv2,j,k)) &
                +sum(divv_yy(i,j,:)*FY(i,j-stv2:j+stv1,k)) &
                +sum(divv_yz(i,j,:)*FZ(i,j,k-stv1:k+stv2)) &
                -yi(j)*sum(interpv_cyl_F_y(i,j,:)*Fcylv(i,j-stv2:j+stv1,k)))
        end do
     end do
  end do
  

  ! ==== W ====
  ! Compute the fluxes
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           FZ(i,j,k) = Hij3(i,j,k)
           
           i = ii; j = jj; k = kk;
           FX(i,j,k) = sum(interp_sc_xz(i,j,:,:)*Hij6(i-st2:i+st1,j,k-st2:k+st1))
           FY(i,j,k) = sum(interp_sc_yz(i,j,:,:)*Hij5(i,j-st2:j+st1,k-st2:k+st1))
        end do
     end do
  end do
  
  ! Compute the source terms
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           srcW(i,j,k) = srcW(i,j,k) + dt_uvw*( &
                +sum(divv_zx(i,j,:)*FX(i-stv1:i+stv2,j,k)) &
                +sum(divv_zy(i,j,:)*FY(i,j-stv1:j+stv2,k)) &
                +sum(divv_zz(i,j,:)*FZ(i,j,k-stv2:k+stv1)) &
                +ymi(j)*sum(interpv_cyl_F_ym(i,j,:)*FY(i,j-stv1:j+stv2,k)))
        end do
     end do
  end do

  ! Set the source terms to zero in the walls
  do j=jmin_,jmax_
     do i=imin_,imax_
        if (mask_u(i,j).ne.0) srcU(i,j,:) = 0.0_WP
        if (mask_v(i,j).ne.0) srcV(i,j,:) = 0.0_WP
        if (mask_w(i,j).ne.0) srcW(i,j,:) = 0.0_WP
     end do
  end do
  
  return
end subroutine sgs_gradient_src_vel


! ======================= !
! Create SGS source terms !
! ======================= !
subroutine sgs_gradient_src_sc(srcSC)
  use sgs_gradient
  use memory
  use filter
  use metric_velocity_visc
  use metric_generic
  use time_info
  use masks
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: srcSC
  integer  :: i,j,k,isc
  
  do isc=1,nscalar
     
     ! Compute the tensor for the gradient model
     call compute_gradient_scalar(isc)
     
     ! Compute the fluxes
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              FX(i,j,k) = sum(interp_sc_x(i,j,:)*Hi1(i-st2:i+st1,j,k))
              FY(i,j,k) = sum(interp_sc_y(i,j,:)*Hi2(i,j-st2:j+st1,k))
              FZ(i,j,k) = sum(interp_sc_z(i,j,:)*Hi3(i,j,k-st2:k+st1))
           end do
        end do
     end do
     
     ! Compute the source terms
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_           
              srcSC(i,j,k,isc) = srcSC(i,j,k,isc) + dt_uvw*( &
                   +sum(div_u(i,j,:)*FX(i-st1:i+st2,j,k)) &
                   +sum(div_v(i,j,:)*FY(i,j-st1:j+st2,k)) &
                   +sum(div_w(i,j,:)*FZ(i,j,k-st1:k+st2)) )
           end do
        end do
     end do

     ! Set the source terms to zero in the walls
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.0) srcSC(i,j,:,isc) = 0.0_WP
        end do
     end do
     
  end do
  
  return
end subroutine sgs_gradient_src_sc


