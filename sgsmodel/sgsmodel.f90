module sgsmodel
  use precision
  use geometry
  use partition
  use data
  use string
  implicit none
  
  ! SGS modeling
  character(len=str_medium) :: avg_model
  logical :: grad_model
  logical :: pope_crit
  
  !                   ( 1 4 6 )
  ! Strain rate : S = ( 4 2 5 )
  !                   ( 6 5 3 )
  real(WP), dimension(:,:,:), pointer :: Sij1
  real(WP), dimension(:,:,:), pointer :: Sij2
  real(WP), dimension(:,:,:), pointer :: Sij3
  real(WP), dimension(:,:,:), pointer :: Sij4
  real(WP), dimension(:,:,:), pointer :: Sij5
  real(WP), dimension(:,:,:), pointer :: Sij6
  real(WP), dimension(:,:,:), pointer :: FS

  !                         ( 1 4 6 )
  ! Total strain rate : S = ( 4 2 5 )
  !                         ( 6 5 3 )
  real(WP), dimension(:,:,:), pointer :: Stij1
  real(WP), dimension(:,:,:), pointer :: Stij2
  real(WP), dimension(:,:,:), pointer :: Stij3
  real(WP), dimension(:,:,:), pointer :: Stij4
  real(WP), dimension(:,:,:), pointer :: Stij5
  real(WP), dimension(:,:,:), pointer :: Stij6
  real(WP), dimension(:,:,:), pointer :: FSt
  
  !                       ( 1 )
  ! Scalar gradient : G = ( 2 )
  !                       ( 3 )
  real(WP), dimension(:,:,:), pointer :: Gi1
  real(WP), dimension(:,:,:), pointer :: Gi2
  real(WP), dimension(:,:,:), pointer :: Gi3
  real(WP), dimension(:,:,:), pointer :: Gnorm
  
  ! Filtered quantities
  real(WP), dimension(:,:,:), pointer :: FUi
  real(WP), dimension(:,:,:), pointer :: FVi
  real(WP), dimension(:,:,:), pointer :: FWi
  real(WP), dimension(:,:,:), pointer :: FUti
  real(WP), dimension(:,:,:), pointer :: FVti
  real(WP), dimension(:,:,:), pointer :: FWti
  real(WP), dimension(:,:,:), pointer :: FRHO
  real(WP), dimension(:,:,:), pointer :: FrhoSC
  
  ! Additional quantities for Mixed Model
  real(WP), dimension(:,:,:), pointer :: Hij1
  real(WP), dimension(:,:,:), pointer :: Hij2
  real(WP), dimension(:,:,:), pointer :: Hij3
  real(WP), dimension(:,:,:), pointer :: Hij4
  real(WP), dimension(:,:,:), pointer :: Hij5
  real(WP), dimension(:,:,:), pointer :: Hij6
  real(WP), dimension(:,:,:), pointer :: Hi1
  real(WP), dimension(:,:,:), pointer :: Hi2
  real(WP), dimension(:,:,:), pointer :: Hi3
  
  ! Numerator (LM) and denominator (MM)
  real(WP), dimension(:,:,:), pointer :: LM
  real(WP), dimension(:,:,:), pointer :: MM
  
  ! Smagorinsky constant for plotting
  real(WP), dimension(:,:,:), pointer :: Cs_visc
  real(WP), dimension(:,:,:), pointer :: Cs_diff
  
  ! Pope's criterion for resolution
  real(WP), dimension(:,:,:), pointer :: k_sgs
  
contains  
  
  ! Compute the Favre averaged velocities at the cell centers (velocity mid-point)
  ! ------------------------------------------------------------------------------
  subroutine compute_filtered_velocity
    use interpolate
    implicit none

    integer :: i,j,k
    
    ! Obtain cell centered U/V/W => Ui/Vi/Wi
    call interpolate_velocities
    
    ! Compute rhoU - cell centered
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FUi(i,j,k) = Ui(i,j,k)*RHOold(i,j,k)
             FVi(i,j,k) = Vi(i,j,k)*RHOold(i,j,k)
             FWi(i,j,k) = Wi(i,j,k)*RHOold(i,j,k)
             FRHO(i,j,k) = RHOold(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Apply filter to the velocities and density
    call filter_global_3D(FUi,FUi,  '+','d')
    call filter_global_3D(FVi,FVi,  '-','d')
    call filter_global_3D(FWi,FWi,  '-','d')
    call filter_global_3D(FRHO,FRHO,'+','n')
    
    ! Divide by filtered density to get Favre average
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FUi(i,j,k) = FUi(i,j,k)/FRHO(i,j,k)
             FVi(i,j,k) = FVi(i,j,k)/FRHO(i,j,k)
             FWi(i,j,k) = FWi(i,j,k)/FRHO(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine compute_filtered_velocity


  ! Compute the Favre averaged velocities at the cell centers (scalar mid-point)
  ! ----------------------------------------------------------------------------
  subroutine compute_filtered_velocity_mid
    use interpolate
    use velocity
    implicit none

    integer :: i,j,k
    
    ! Obtain cell centered U/V/W => Ui/Vi/Wi
    call interpolate_velocities
    
    ! Compute rhoU - cell centered
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FUi(i,j,k) = Ui(i,j,k)*RHOmid(i,j,k)
             FVi(i,j,k) = Vi(i,j,k)*RHOmid(i,j,k)
             FWi(i,j,k) = Wi(i,j,k)*RHOmid(i,j,k)
             FRHO(i,j,k) = RHOmid(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Apply filter to the velocities and density
    call filter_global_3D(FUi,FUi,  '+','d')
    call filter_global_3D(FVi,FVi,  '-','d')
    call filter_global_3D(FWi,FWi,  '-','d')
    call filter_global_3D(FRHO,FRHO,'+','n')
    
    ! Divide by filtered density to get Favre average
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FUi(i,j,k) = FUi(i,j,k)/FRHO(i,j,k)
             FVi(i,j,k) = FVi(i,j,k)/FRHO(i,j,k)
             FWi(i,j,k) = FWi(i,j,k)/FRHO(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine compute_filtered_velocity_mid


  ! Compute the Favre averaged total velocities at the cell centers
  ! ---------------------------------------------------------------
  subroutine compute_filtered_velocity_total(isc)
    use interpolate
    implicit none

    integer, intent(in) :: isc
    integer :: i,j,k
    
    ! Obtain cell centered Ut/Vt/Wt => Uti/Vti/Wti
    call interpolate_velocities_total(isc)
    
    ! Compute rhoUt - cell centered
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FUti(i,j,k) = Uti(i,j,k)*RHO(i,j,k)
             FVti(i,j,k) = Vti(i,j,k)*RHO(i,j,k)
             FWti(i,j,k) = Wti(i,j,k)*RHO(i,j,k)
             FRHO(i,j,k) = RHO(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Apply filter to the velocities and density
    call filter_global_3D(FUti,FUti,  '+','d')
    call filter_global_3D(FVti,FVti,  '-','d')
    call filter_global_3D(FWti,FWti,  '-','d')
    call filter_global_3D(FRHO,FRHO,'+','n')
    
    ! Divide by filtered density to get Favre average
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FUti(i,j,k) = FUti(i,j,k)/FRHO(i,j,k)
             FVti(i,j,k) = FVti(i,j,k)/FRHO(i,j,k)
             FWti(i,j,k) = FWti(i,j,k)/FRHO(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine compute_filtered_velocity_total
  
  
  ! Compute the Favre averaged scalar at the cell centers
  ! -----------------------------------------------------
  subroutine compute_filtered_scalar(isc)
    implicit none
    
    integer, intent(in) :: isc
    integer :: i,j,k
    
    ! Compute rhoSC - cell centered
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FrhoSC(i,j,k) = RHO(i,j,k)*SC(i,j,k,isc)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Apply filter to the scalar
    call filter_global_3D(FrhoSC,FrhoSC,'+','n')
    
    return
  end subroutine compute_filtered_scalar
  
  
  ! Compute the gradient tensor of velocities
  ! -----------------------------------------
  subroutine compute_Hij
    implicit none
    
    Hij1 = 0.0_WP
    Hij2 = 0.0_WP
    Hij3 = 0.0_WP
    Hij4 = 0.0_WP
    Hij5 = 0.0_WP
    Hij6 = 0.0_WP
    
    ! If gradient model do something
    if (grad_model) call sgs_gradient_compute_Hij
    
    return
  end subroutine compute_Hij
  
  
  ! Compute the gradient vector of scalar
  ! -----------------------------------------
  subroutine compute_Hi(isc)
    implicit none
    integer, intent(in) :: isc
    
    Hi1 = 0.0_WP
    Hi2 = 0.0_WP
    Hi3 = 0.0_WP
    
    ! If gradient model do something
    if (grad_model) call sgs_gradient_compute_Hi(isc)
    
    return
  end subroutine compute_Hi
  
  
  ! Compute local LM and MM for Eddy Viscosity
  ! ------------------------------------------
  subroutine sgsmodel_LM_MM_VISC
    use interpolate
    use strainrate
    use filter
    use masks
    implicit none
    
    real(WP), dimension(-1:+1,-1:+1,-1:+1) :: buffer1,buffer2,buffer3,buffer4
    real(WP), dimension(6) :: Mij,Lij,FrhoUiUj,FrhoSSij,FSij
    integer :: i,j,k
    
    ! Compute the Favre averaged velocities at the cell centers
    call compute_filtered_velocity
    
    ! Compute Hij for Mixed Gradient Model
    call compute_Hij
    
    ! Compute the strain rate
    call strainrate_compute(Sij1,Sij2,Sij3,Sij4,Sij5,Sij6)
    
    !$OMP PARALLEL DO PRIVATE(buffer1,buffer2,buffer3,buffer4,FSij,FrhoSSij,Mij,FrhoUiUj,Lij)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             
             ! Filter Sij
             buffer1 = Sij1(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FSij(1),i,j,'n')
             buffer1 = Sij2(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FSij(2),i,j,'n')
             buffer1 = Sij3(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FSij(3),i,j,'n')
             buffer1 = Sij4(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FSij(4),i,j,'n')
             buffer1 = Sij5(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FSij(5),i,j,'n')
             buffer1 = Sij6(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FSij(6),i,j,'n')
             
             ! Norm of FSij
             FS(i,j,k) = sqrt(sum(FSij(1:3)**2 + 2.0_WP*FSij(4:6)**2))
             
             ! Filter the product rho.S.Sij
             buffer1 = S(i-1:i+1,j-1:j+1,k-1:k+1) * RHOold(i-1:i+1,j-1:j+1,k-1:k+1)
             buffer2 = buffer1*Sij1(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2,FrhoSSij(1),i,j,'n')
             buffer2 = buffer1*Sij2(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2,FrhoSSij(2),i,j,'n')
             buffer2 = buffer1*Sij3(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2,FrhoSSij(3),i,j,'n')
             buffer2 = buffer1*Sij4(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2,FrhoSSij(4),i,j,'n')
             buffer2 = buffer1*Sij5(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2,FrhoSSij(5),i,j,'n')
             buffer2 = buffer1*Sij6(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2,FrhoSSij(6),i,j,'n')
             
             ! Compute Mij = <rho.S.Sij> - ratio^2.<rho>.<S>.<Sij>
             Mij = 2.0_WP*delta_3D(i,j)**2 * (FrhoSSij - ratio_3D(i,j)**2 *FRHO(i,j,k)*FS(i,j,k)*FSij)
             
             ! Filter the product rho.Ui.Uj
             buffer1 = RHOold(i-1:i+1,j-1:j+1,k-1:k+1) * Ui(i-1:i+1,j-1:j+1,k-1:k+1)
             buffer2 = RHOold(i-1:i+1,j-1:j+1,k-1:k+1) * Vi(i-1:i+1,j-1:j+1,k-1:k+1)
             buffer3 = RHOold(i-1:i+1,j-1:j+1,k-1:k+1) * Wi(i-1:i+1,j-1:j+1,k-1:k+1)
             buffer4 = buffer1 * Ui(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer4, FrhoUiUj(1),i,j,'n')
             buffer4 = buffer2 * Vi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer4, FrhoUiUj(2),i,j,'n')
             buffer4 = buffer3 * Wi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer4, FrhoUiUj(3),i,j,'n')
             buffer4 = buffer1 * Vi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer4, FrhoUiUj(4),i,j,'n')
             buffer4 = buffer2 * Wi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer4, FrhoUiUj(5),i,j,'n')
             buffer4 = buffer1 * Wi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer4, FrhoUiUj(6),i,j,'n')
             
             ! Compute Lij = <rho.Ui.Uj> - <rho>.<Ui>.<Uj>
             Lij(1) = FrhoUiUj(1) - FRHO(i,j,k) * FUi(i,j,k) * FUi(i,j,k)
             Lij(2) = FrhoUiUj(2) - FRHO(i,j,k) * FVi(i,j,k) * FVi(i,j,k)
             Lij(3) = FrhoUiUj(3) - FRHO(i,j,k) * FWi(i,j,k) * FWi(i,j,k)
             Lij(4) = FrhoUiUj(4) - FRHO(i,j,k) * FUi(i,j,k) * FVi(i,j,k)
             Lij(5) = FrhoUiUj(5) - FRHO(i,j,k) * FVi(i,j,k) * FWi(i,j,k)
             Lij(6) = FrhoUiUj(6) - FRHO(i,j,k) * FUi(i,j,k) * FWi(i,j,k)
             
             ! Remove Hij from Lij
             Lij(1) = Lij(1) - Hij1(i,j,k)
             Lij(2) = Lij(2) - Hij2(i,j,k)
             Lij(3) = Lij(3) - Hij3(i,j,k)
             Lij(4) = Lij(4) - Hij4(i,j,k)
             Lij(5) = Lij(5) - Hij5(i,j,k)
             Lij(6) = Lij(6) - Hij6(i,j,k)
             
             ! Compute LijMij = sum Lij.Mij
             LM(i,j,k) = sum (Lij(1:3)*Mij(1:3) + 2.0_WP*Lij(4:6)*Mij(4:6) )
             
             ! Compute MijMij = sum Mij.Mij
             MM(i,j,k) = sum (Mij(1:3)*Mij(1:3) + 2.0_WP*Mij(4:6)*Mij(4:6) )
             
             ! Force 0 in walls
             if (mask(i,j).ne.0) then
                LM(i,j,k) = 0.0_WP
                MM(i,j,k) = 0.0_WP
             end if
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine sgsmodel_LM_MM_VISC
  
  
  ! Compute local LM and MM for Eddy Diffusivity
  ! --------------------------------------------
  subroutine sgsmodel_LM_MM_DIFF(isc)
    use interpolate
    use strainrate
    use filter
    use masks
    use velocity
    implicit none
    
    integer, intent(in) :: isc
    real(WP), dimension(-1:+1,-1:+1,-1:+1) :: buffer1,buffer2
    real(WP), dimension(3) :: Li,Mi,FGi,FrhoStGi,FrhoUtiSC
    real(WP), dimension(6) :: FStij
    integer :: i,j,k
    
    ! Compute the Favre averaged scalar at the cell centers
    call compute_filtered_scalar(isc)

    ! Compute the Favre averaged total velocity at the cell centers
    if (isc.eq.1) call compute_filtered_velocity_mid

    ! Compute Hi for Mixed Gradient Model
    call compute_Hi(isc)

    ! Compute the total strain rate
    if (isc.eq.1) call strainrate_compute(Stij1,Stij2,Stij3,Stij4,Stij5,Stij6)
!!$    call strainrate_compute(Stij1,Stij2,Stij3,Stij4,Stij5,Stij6,isc)
    
    ! Compute the scalar gradient
    call gradient_vector(SC(:,:,:,isc),Gi1,Gi2,Gi3)
    
    !$OMP PARALLEL DO PRIVATE(buffer1,buffer2,FStij,FGi,FrhoStGi,Mi,FrhoUtiSC,Li)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             
             ! Filter Stij
             buffer1 = Stij1(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FStij(1),i,j,'n')
             buffer1 = Stij2(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FStij(2),i,j,'n')
             buffer1 = Stij3(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FStij(3),i,j,'n')
             buffer1 = Stij4(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FStij(4),i,j,'n')
             buffer1 = Stij5(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FStij(5),i,j,'n')
             buffer1 = Stij6(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FStij(6),i,j,'n')

             ! Norm of FStij
             FSt(i,j,k) = sqrt(sum(FStij(1:3)**2 + 2.0_WP*FStij(4:6)**2))

             ! Filter Gi
             buffer1 = Gi1(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FGi(1),i,j,'d')
             buffer1 = Gi2(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FGi(2),i,j,'d')
             buffer1 = Gi3(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FGi(3),i,j,'d')
             
             ! Filter the product rho.St.Gi
!!$             buffer1 = St(i-1:i+1,j-1:j+1,k-1:k+1) * RHOmid(i-1:i+1,j-1:j+1,k-1:k+1)
             buffer1 = S(i-1:i+1,j-1:j+1,k-1:k+1) * RHOmid(i-1:i+1,j-1:j+1,k-1:k+1)
             buffer2 = buffer1*Gi1(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoStGi(1),i,j,'d')
             buffer2 = buffer1*Gi2(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoStGi(2),i,j,'d')
             buffer2 = buffer1*Gi3(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoStGi(3),i,j,'d')
             
             ! Compute Mi = <rho.St.Gi> - ratio^2.<rho>.<St>.<Gi>
             Mi = 2.0_WP * delta_3D(i,j)**2 * (FrhoStGi - ratio_3D(i,j)**2  * FRHO(i,j,k) * FSt(i,j,k) * FGi)
!!$             Mi = 2.0_WP * delta_3D(i,j)**2 * (FrhoStGi - ratio_3D(i,j)**2  * FRHO(i,j,k) * FS(i,j,k) * FGi)
             
             ! Filter the product rho.Uti.SC
             ! -> Has to be Dirichlet here for consistency with FUti,FVti,FWti
             buffer1 = RHOmid(i-1:i+1,j-1:j+1,k-1:k+1) * SC(i-1:i+1,j-1:j+1,k-1:k+1,isc)
!!$             buffer2 = buffer1*Uti(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoUtiSC(1),i,j,'d')
!!$             buffer2 = buffer1*Vti(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoUtiSC(2),i,j,'d')
!!$             buffer2 = buffer1*Wti(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoUtiSC(3),i,j,'d')
             buffer2 = buffer1*Ui(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoUtiSC(1),i,j,'d')
             buffer2 = buffer1*Vi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoUtiSC(2),i,j,'d')
             buffer2 = buffer1*Wi(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer2, FrhoUtiSC(3),i,j,'d')
             
             ! Compute Li = <rho.Uti.SC> - <rho>.<Uti>.<SC>
!!$             Li(1) = FrhoUtiSC(1) - FrhoSC(i,j,k) * FUti(i,j,k)
!!$             Li(2) = FrhoUtiSC(2) - FrhoSC(i,j,k) * FVti(i,j,k)
!!$             Li(3) = FrhoUtiSC(3) - FrhoSC(i,j,k) * FWti(i,j,k)
             Li(1) = FrhoUtiSC(1) - FrhoSC(i,j,k) * FUi(i,j,k)
             Li(2) = FrhoUtiSC(2) - FrhoSC(i,j,k) * FVi(i,j,k)
             Li(3) = FrhoUtiSC(3) - FrhoSC(i,j,k) * FWi(i,j,k)
             
             ! Remove Hi from Li
             Li(1) = Li(1) - Hi1(i,j,k)
             Li(2) = Li(2) - Hi2(i,j,k)
             Li(3) = Li(3) - Hi3(i,j,k)
             
             ! Compute LiMi = sum Li.Mi
             LM(i,j,k) = sum(Li*Mi)
             
             ! Compute MiMi = sum Mi.Mi
             MM(i,j,k) = sum(Mi*Mi)

             ! Force 0 in walls
             if (mask(i,j).ne.0) then
                LM(i,j,k) = 0.0_WP
                MM(i,j,k) = 0.0_WP
             end if
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine sgsmodel_LM_MM_DIFF
  
  
  ! Compute local LM and MM for Scalar Variance
  ! --------------------------------------------
  subroutine sgsmodel_LM_MM_ZVAR(isc)
    use interpolate
    use strainrate
    use filter
    use masks
    implicit none
    
    integer, intent(in) :: isc
    real(WP), dimension(-1:+1,-1:+1,-1:+1) :: buffer1
    real(WP), dimension(3) :: FGi
    real(WP) :: Li,Mi,FrhoSC2,FrhoFSC2,FSC
    real(WP) :: FGnorm,FrhoGnorm
    integer  :: i,j,k
    
    ! Compute the Favre averaged scalar at the cell centers
    call compute_filtered_scalar(isc)
    call filter_global_3D(RHO,FRHO,'+','n')
    
    ! Compute the gradient
    call gradient_vector(SC(:,:,:,isc),Gi1,Gi2,Gi3)

    !$OMP PARALLEL PRIVATE(buffer1,FGi,FGnorm,FrhoGnorm,FSC,FrhoFSC2,FrhoSC2,Li,Mi)

    !$OMP DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             Gnorm(i,j,k) = Gi1(i,j,k)**2 + Gi2(i,j,k)**2 + Gi3(i,j,k)**2
          end do
       end do
    end do
    !$OMP END DO
    
    !$OMP DO
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             
             ! Filter Gi
             buffer1 = Gi1(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FGi(1),i,j,'n')
             buffer1 = Gi2(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FGi(2),i,j,'n')
             buffer1 = Gi3(i-1:i+1,j-1:j+1,k-1:k+1); call filter_local_3D(buffer1,FGi(3),i,j,'n')
             
             ! Compute norm^2
             FGnorm = sum(FGi**2)
             
             ! Filter the product rho.Gnorm
             buffer1 = Gnorm(i-1:i+1,j-1:j+1,k-1:k+1) * RHO(i-1:i+1,j-1:j+1,k-1:k+1)
             call filter_local_3D(buffer1,FrhoGnorm,i,j,'n')
             
             ! Calculate the model term
             !Mi = delta_3D(i,j)**2*(ratio_3D(i,j)**2*FRHO(i,j,k)*FGnorm-FrhoGnorm)
             Mi = delta_3D(i,j)**2*ratio_3D(i,j)**2*FRHO(i,j,k)*FGnorm!-FrhoGnorm)
             
             ! Build Frho*FSC**2
             FSC = FrhoSC(i,j,k)/FRHO(i,j,k)
             FrhoFSC2 = FRHO(i,j,k)*FSC**2
             
             ! Filter rho*SC**2
             buffer1 = RHO(i-1:i+1,j-1:j+1,k-1:k+1) * SC(i-1:i+1,j-1:j+1,k-1:k+1,isc)**2
             call filter_local_3D(buffer1,FrhoSC2,i,j,'n')
             
             ! Calculate the Leonard term
             Li = FrhoSC2-FrhoFSC2
             
             ! Compute LiMi = Li.Mi
             LM(i,j,k) = Li*Mi
             
             ! Compute MiMi = Mi.Mi
             MM(i,j,k) = Mi*Mi
             
             ! Force 0 in walls
             if (mask(i,j).ne.0) then
                LM(i,j,k) = 0.0_WP
                MM(i,j,k) = 0.0_WP
             end if
             
          end do
       end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL
    
    return
  end subroutine sgsmodel_LM_MM_ZVAR
  
  
  ! Apply the boundary conditions
  ! -----------------------------
  subroutine sgsmodel_apply_border(A)
    implicit none
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
    
    ! Update Ghost Cells
    call boundary_update_border(A,'+','ym')
    
    ! Take care of non periodic BCs
    call boundary_neumann(A,'-xm')
    call boundary_neumann(A,'+xm')
    call boundary_neumann(A,'-ym')
    call boundary_neumann(A,'+ym')
    
    return
  end subroutine sgsmodel_apply_border
  
end module sgsmodel


! ======================== !
! Initialize the SGS model !
! ======================== !
subroutine sgsmodel_init
  use sgsmodel
  use parser
  use memory
  implicit none
  
  logical :: default_pope
  
  ! Determine what to do
  call parser_read('Use SGS model',use_sgs)
  if (.not.use_sgs) then
     default_pope=.false.
  else
     default_pope=.true.
  end if
  call parser_read('Pope criterion',pope_crit,default_pope)

  ! Nothing to do return
  if (.not.use_sgs .and. .not.pope_crit) return

  ! Do we use a term for scale similarity?
  call parser_read('Use gradient model',grad_model,.false.)
  
  ! Read in the averaging method
  call parser_read('SGS averaging',avg_model)
  
  ! Create & Start the timer
  call timing_create('sgsmodel')
  call timing_start ('sgsmodel')
  
  ! Allocate work arrays
  allocate(FUi   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FVi   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FWi   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FUti  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FVti  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FWti  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FRHO  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FS    (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FSt   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FRHOSC(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(LM    (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(MM    (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Cs_visc(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Cs_diff(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  if (nscalar.ge.1) allocate(Gnorm(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  if (pope_crit) allocate(k_sgs(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Link temporary arrays
  Sij1 =>tmp1;  Sij2 =>tmp2;  Sij3 =>tmp3
  Sij4 =>tmp4;  Sij5 =>tmp5;  Sij6 =>tmp6
  Gi1  =>tmp1;  Gi2  =>tmp2;  Gi3  =>tmp3
  Stij1=>tmp4;  Stij2=>tmp5;  Stij3=>tmp6
  Stij4=>tmp7;  Stij5=>tmp8;  Stij6=>tmp9
  Hij1 =>tmp10; Hij2 =>tmp11; Hij3 =>tmp12
  Hij4 =>tmp13; Hij5 =>tmp14; Hij6 =>tmp15
  Hi1  =>tmp10; Hi2  =>tmp11; Hi3  =>tmp12

  ! Specific initialization for each model
  select case(trim(avg_model))
  case ('Germano')
     call sgs_germano_init
  case ('Ad-Hoc')
     call sgs_adhoc_init
  case ('Lagrangian')
     call sgs_lagrangian_init
  case default
     call die ('sgsmodel_init : unknown SGS model')
  end select
  
  ! Initialize the mixed models
  if (grad_model) call sgs_gradient_init
  
  ! Create a new file to monitor at each timestep
  call monitor_create_file_step('sgsmodel',1)
  call monitor_set_header (1,'K_SGS','r')
  
  ! Stop the timer
  call timing_stop('sgsmodel')
  
  return
end subroutine sgsmodel_init


! ========================== !
! Compute the Eddy Viscosity !
! ========================== !
subroutine sgsmodel_eddyVISC
  use sgsmodel
  use strainrate
  use filter
  use velocity
  implicit none
  integer  :: i,j,k
  
  ! Nothing to be done if no SGS
  if (.not.use_sgs .and. .not.pope_crit) return
  
  ! Start the timer
  call timing_start ('sgsmodel')

  ! Compute velocity at mid-point
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           U(i,j,k) = 0.5_WP * (U(i,j,k) + Uold(i,j,k))
           V(i,j,k) = 0.5_WP * (V(i,j,k) + Vold(i,j,k))
           W(i,j,k) = 0.5_WP * (W(i,j,k) + Wold(i,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Compute the local LM and MM
  call sgsmodel_LM_MM_VISC
  
  ! Perform the average
  select case(trim(avg_model))
  case ('Germano')
     call sgs_germano_average_VISC
  case ('Ad-Hoc')
     call sgs_adhoc_average_VISC
  case ('Lagrangian')
     call sgs_lagrangian_average_VISC
  end select
  
  ! Compute the Eddy Viscosity
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           if (MM(i,j,k).ne.0.0_WP) then
              Cs_visc(i,j,k) = LM(i,j,k)/MM(i,j,k)
           else
              Cs_visc(i,j,k) = 0.0_WP
           end if
           !VISC(i,j,k) = RHO(i,j,k) * S(i,j,k) * Cs_visc(i,j,k) * delta_3D(i,j)**2
           VISC(i,j,k) = RHOold(i,j,k) * S(i,j,k) * Cs_visc(i,j,k) * delta_3D(i,j)**2
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Compute Pope's criterion for LES
  if (pope_crit) then
     call sgsmodel_pope_visc
  end if
  
  ! Clip and apply
  if (use_sgs) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              if (VISC(i,j,k).lt.-VISCmol(i,j,k)) VISC(i,j,k) = -VISCmol(i,j,k)
              VISC(i,j,k) = VISCmol(i,j,k) + VISC(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              VISC(i,j,k) = VISCmol(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Second clipping
  !call sgsmodel_pdf_clip(VISC)

  ! Return the correct velocity
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           U(i,j,k) = 2.0_WP * U(i,j,k) - Uold(i,j,k)
           V(i,j,k) = 2.0_WP * V(i,j,k) - Vold(i,j,k)
           W(i,j,k) = 2.0_WP * W(i,j,k) - Wold(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Stop the timer
  call timing_stop('sgsmodel')
  
  return
end subroutine sgsmodel_eddyVISC


! ============================ !
! Compute the Eddy Diffusivity !
! ============================ !
subroutine sgsmodel_eddyDIFF
  use sgsmodel
  use strainrate
  use filter
  use scalar
  implicit none
  integer  :: i,j,k,isc
  
  ! Nothing to be done if no SGS
  if (.not.use_sgs .or. nscalar.eq.0) return
  
  ! Start the timer
  call timing_start ('sgsmodel')
  
  ! Development Note:
  ! -- For the diffusivity, we need to use the total velocities, which requires an update prior to this
  ! -- computation.  I need to think about the structure of this more...

  ! Compute scalars at mid-point
  !$OMP PARALLEL
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = 0.5_WP * (SC(i,j,k,isc) + SCold(i,j,k,isc))
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL

  do isc=1,nscalar
     
     ! Compute the local LM and MM
     call sgsmodel_LM_MM_DIFF(isc)
     
     ! Perform the average
     select case(trim(avg_model))
     case ('Germano')
        call sgs_germano_average_DIFF
     case ('Ad-Hoc')
        call sgs_adhoc_average_DIFF
     case ('Lagrangian')
        call sgs_lagrangian_average_DIFF(isc)
     end select
     
     ! Compute the Eddy Diffusivity
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              if (MM(i,j,k).ne.0.0_WP) then
                 Cs_diff(i,j,k) = LM(i,j,k)/MM(i,j,k)
              else
                 Cs_diff(i,j,k) = 0.0_WP
              end if
!!$              DIFF(i,j,k,isc) = RHO(i,j,k) * St(i,j,k) * Cs_diff(i,j,k) * delta_3D(i,j)**2
              DIFF(i,j,k,isc) = RHO(i,j,k) * S(i,j,k) * Cs_diff(i,j,k) * delta_3D(i,j)**2
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     
  end do
  
  !$OMP PARALLEL

  ! Clip and apply
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              if (DIFF(i,j,k,isc).lt.-DIFFmol(i,j,k,isc)) DIFF(i,j,k,isc) = -DIFFmol(i,j,k,isc)
              DIFF(i,j,k,isc) = DIFFmol(i,j,k,isc) + DIFF(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL

  ! Second clipping
  !do isc=1,nscalar
  !   call sgsmodel_pdf_clip(DIFF(:,:,:,isc))
  !end do

  ! Return the correct scalars
  !$OMP PARALLEL
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = 2.0_WP * SC(i,j,k,isc) - SCold(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL
  
  ! Stop the timer
  call timing_stop('sgsmodel')
  
  return
end subroutine sgsmodel_eddyDIFF


! ======================== !
! Compute the SGS Variance !
! ======================== !
subroutine sgsmodel_ZVAR(ZVAR,isc)
  use sgsmodel
  use filter
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: ZVAR
  integer, intent(in) :: isc
  integer  :: i,j,k
  real(WP) :: Cs
  
  ! Nothing to be done if no SGS
  if (.not.use_sgs) return
  
  ! Start the timer
  call timing_start('sgsmodel')
  
  ! Compute the local LM and MM
  call sgsmodel_LM_MM_ZVAR(isc)
  
  ! Perform the average
  select case(trim(avg_model))
  case ('Germano')
     call sgs_germano_average_ZVAR
  case ('Ad-Hoc')
     call sgs_adhoc_average_ZVAR
  case ('Lagrangian')
     call sgs_lagrangian_average_ZVAR
  end select
  
  ! Compute the Variance
  !$OMP PARALLEL DO PRIVATE(Cs)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           if (MM(i,j,k).ne.0.0_WP) then
              Cs = LM(i,j,k)/MM(i,j,k)
           else
              Cs = 0.0_WP
           end if
           Cs = max(Cs,0.0_WP)
           ZVAR(i,j,k) = Gnorm(i,j,k) * Cs * delta_3D(i,j)**2
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Stop the timer
  call timing_stop('sgsmodel')
  
  return
end subroutine sgsmodel_ZVAR


! ======================= !
! Create SGS source terms !
! ======================= !
subroutine sgsmodel_src_vel(srcU,srcV,srcW)
  use sgsmodel
  implicit none
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: srcU,srcV,srcW

  ! Nothing to be done if no SGS
  if (.not.use_sgs) return
  
  ! Start the timer
  call timing_start ('sgsmodel')
  
  ! Compute the source terms for the mixed model
  if (grad_model) call sgs_gradient_src_vel(srcU,srcV,srcW)
  
  ! Stop the timer
  call timing_stop('sgsmodel')
  
  return
end subroutine sgsmodel_src_vel


! ======================= !
! Create SGS source terms !
! ======================= !
subroutine sgsmodel_src_sc(srcSC)
  use sgsmodel
  implicit none
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: srcSC

  ! Nothing to be done if no SGS
  if (.not.use_sgs) return
  
  ! Start the timer
  call timing_start ('sgsmodel')
  
  ! Compute the source terms for the mixed model
  if (grad_model) call sgs_gradient_src_sc(srcSC)
  
  ! Stop the timer
  call timing_stop('sgsmodel')
  
  return
end subroutine sgsmodel_src_sc


! ================================== !
! Clip the quantity based on its PDF !
! ================================== !
subroutine sgsmodel_pdf_clip(A)
  use sgsmodel
  use masks
  implicit none
  
  integer, parameter :: nbins=100
  real(WP), parameter :: clip=0.01_WP
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), dimension(nbins) :: bins,my_pdf,pdf
  real(WP) :: width,mini,maxi,tot,limit,clipped
  integer :: i,j,k,ibin,npoints,my_npoints,nclipped,my_nclipped
  
  ! Prepare pdf of A
  call parallel_min(minval(A),mini)
  call parallel_max(maxval(A),maxi)
  width = (maxi-mini)/real(nbins,WP)
  if (width.eq.0.0_WP) return
  do i=1,nbins
     bins(i) = mini+real(i-1,WP)*width
  end do
  
  ! Compute pdf of A
  my_pdf = 0.0_WP
  my_npoints = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              ! Count the points
              my_npoints = my_npoints+1
              ! Find bin
              ibin = min(floor((A(i,j,k)-mini)/width)+1,nbins)
              ! Increment it
              my_pdf(ibin) = my_pdf(ibin) + 1.0_WP
           end if
        end do
     end do
  end do
  call parallel_sum(my_pdf,pdf)
  call parallel_sum(my_npoints,npoints)
  pdf = pdf/real(npoints,WP)
  
  ! Find the value of A
  tot = 0.0_WP
  i = nbins
  do while (tot<clip)
     tot = tot+pdf(i)
     i = i-1
  end do
  limit = bins(i+1)
  
  ! Clip above limit
  my_nclipped = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (A(i,j,k).gt.limit) then
              A(i,j,k) = limit
              my_nclipped = my_nclipped + 1
           end if
        end do
     end do
  end do
  call parallel_sum(my_nclipped,nclipped)
  clipped = real(nclipped,WP)/real(npoints,WP)
  
  return
end subroutine sgsmodel_pdf_clip


! =============================== !
! Pope's criterion for resolution !
! =============================== !
subroutine sgsmodel_pope_visc
  use sgsmodel
  use math
  use filter
  implicit none
  
  integer  :: i,j,k
  real(WP) :: Ck,Cs,ksgs,buf1,buf2
  real(WP), parameter :: CKol=1.5_WP
  
  ! Compute Ksgs
  !$OMP PARALLEL DO PRIVATE(Cs,Ck)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           ! Compute coefficient
           if (MM(i,j,k).ne.0.0_WP) then
              Cs = max(LM(i,j,k)/MM(i,j,k),0.0_WP)
              Ck = pi**(1.0_WP/3.0_WP) * (2.0_WP/(3.0_WP*CKol))**(1.0_WP/2.0_WP) * Cs**(4.0_WP/3.0_WP)
              k_sgs(i,j,k) = VISC(i,j,k)**2 / ( RHO(i,j,k)**2 * Ck**2 * delta_3D(i,j)**2 + epsilon(1.0_WP) )
           else
              Cs = 0.0_WP
              k_sgs(i,j,k) = -1.0_WP
           end if
           ! Account for negative VISC
           if (VISC(i,j,k).le.0.0_WP) k_sgs(i,j,k) = -1.0_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Monitoring of ksgs
  buf1 = 0.0_WP
  buf2 = 0.0_WP
  !$OMP PARALLEL DO REDUCTION(+:buf1,buf2)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + vol(i,j) * k_sgs(i,j,k) * RHO(i,j,k)
           buf2 = buf2 + RHO(i,j,k) * vol(i,j)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call parallel_sum(buf1,ksgs)
  call parallel_sum(buf2,buf1)
  ksgs = ksgs/buf1
  call monitor_select_file('sgsmodel')
  call monitor_set_single_value(1,ksgs)
  
  return
end subroutine sgsmodel_pope_visc
