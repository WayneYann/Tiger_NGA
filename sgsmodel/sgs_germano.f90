module sgs_germano
  use sgsmodel
  implicit none
  
  ! Dummy module
  
contains
  
  ! Average both arrays over the periodic directions
  ! ------------------------------------------------
  subroutine sgs_germano_average(A)
    implicit none
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
    integer  :: i,j
    real(WP) :: A_sum
    
    if (xper.eq.1) then
       if (yper.eq.1) then
          ! XYZ averaging
          call parallel_sum(sum(A(imin_:imax_,jmin_:jmax_,kmin_:kmax_)),A_sum)
          A = A_sum/real(nx*ny*nz,WP)
       else
          ! XZ averaging
          do j=jmin_,jmax_
             call parallel_sum_dir(sum(A(imin_:imax_,j,kmin_:kmax_)),A_sum,"xz")
             A(:,j,:) = A_sum/real(nx*nz,WP)
          end do
       end if
    else
       if (yper.eq.1) then
          ! YZ averaging
          do i=imin_,imax_
             call parallel_sum_dir(sum(A(i,jmin_:jmax_,kmin_:kmax_)),A_sum,"yz")
             A(i,:,:) = A_sum/real(ny*nz,WP)
          end do
       else
          ! Z averaging
          do j=jmin_,jmax_
             do i=imin_,imax_
                call parallel_sum_dir(sum(A(i,j,kmin_:kmax_)),A_sum,"z")
                A(i,j,:) = A_sum/real(nz,WP)
             end do
          end do
       end if
    end if
    
    ! Apply BCs
    call sgsmodel_apply_border(A)
    
    return
  end subroutine sgs_germano_average
  
end module sgs_germano


! ======================== !
! Initialize the SGS model !
! ======================== !
subroutine sgs_germano_init
  use sgs_germano
  implicit none
  
  ! Dummy init
  
  return
end subroutine sgs_germano_init


! ============================================== !
! Average LM and MM over time for Eddy Viscosity !
! ============================================== !
subroutine sgs_germano_average_VISC
  use sgs_germano
  implicit none
  
  call sgs_germano_average(LM)
  call sgs_germano_average(MM)
  
  return
end subroutine sgs_germano_average_VISC


! ================================================ !
! Average LM and MM over time for Eddy Diffusivity !
! ================================================ !
subroutine sgs_germano_average_DIFF
  use sgs_germano
  implicit none
  
  call sgs_germano_average(LM)
  call sgs_germano_average(MM)
  
  return
end subroutine sgs_germano_average_DIFF


! =============================================== !
! Average LM and MM over time for Scalar Variance !
! =============================================== !
subroutine sgs_germano_average_ZVAR
  use sgs_germano
  implicit none
  
  call sgs_germano_average(LM)
  call sgs_germano_average(MM)
  
  return
end subroutine sgs_germano_average_ZVAR
