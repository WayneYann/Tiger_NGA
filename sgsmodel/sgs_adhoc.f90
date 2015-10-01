module sgs_adhoc
  use sgsmodel
  implicit none
  
  integer :: nfilter
  
end module sgs_adhoc


! ============================== !
! Initialize the Adhoc SGS model !
! ============================== !
subroutine sgs_adhoc_init
  use sgs_adhoc
  use parser
  implicit none
  
  ! Read the number of filtering operations
  call parser_read('SGS Ad-Hoc filtering',nfilter,3)  
  
  return
end subroutine sgs_adhoc_init


! ============================================== !
! Average LM and MM over time for Eddy Viscosity !
! ============================================== !
subroutine sgs_adhoc_average_VISC
  use sgs_adhoc
  implicit none
  integer :: n

  ! The filter operation requires values in the ghost cells
  call boundary_update_border(LM,'ym','+')
  call boundary_update_border(MM,'ym','+')
  
  do n=1,nfilter
     call filter_global_3D(LM,LM,'+','n')
     call filter_global_3D(MM,MM,'+','n')
  end do
  
  return
end subroutine sgs_adhoc_average_VISC


! ================================================ !
! Average LM and MM over time for Eddy Diffusivity !
! ================================================ !
subroutine sgs_adhoc_average_DIFF
  use sgs_adhoc
  implicit none
  integer :: n
  
  ! The filter operation requires values in the ghost cells
  call boundary_update_border(LM,'ym','+')
  call boundary_update_border(MM,'ym','+')
  
  do n=1,nfilter
     call filter_global_3D(LM,LM,'+','d')
     call filter_global_3D(MM,MM,'+','d')
  end do
  
  return
end subroutine sgs_adhoc_average_DIFF


! =============================================== !
! Average LM and MM over time for Scalar Variance !
! =============================================== !
subroutine sgs_adhoc_average_ZVAR
  use sgs_adhoc
  implicit none
  integer :: n
  
  ! The filter operation requires values in the ghost cells
  call boundary_update_border(LM,'ym','+')
  call boundary_update_border(MM,'ym','+')
  
  do n=1,nfilter
     call filter_global_3D(LM,LM,'+','n')
     call filter_global_3D(MM,MM,'+','n')
  end do
  
  return
end subroutine sgs_adhoc_average_ZVAR
