module chemtable
  use combustion
  use string
  implicit none

  ! Chemtable dimensionality
  integer :: combDim

  ! Combustion model
  character(len=str_medium) :: combModel

  ! Density method
  character(len=str_medium) :: density_model

end module chemtable


! ================================================= !
! Initialize tabulated chemistry                    !
! -- Read table with specified number of dimensions !
! -- Make sure the requisite scalars are present    !
! ================================================= !
subroutine chemtable_init(filename)
  use chemtable
  use parser
  use parallel
  implicit none

  character(len=str_medium), intent(inout) :: filename

  ! Get the dimensionality of the chemtable
  call parser_read('Chemtable dimensions',combDim)

  ! Initialize the appropriate chemtable type
  select case(combDim)
  case (1)
     call chemtable_1D_init(filename,combModel)
  case (2)
     call chemtable_2D_init(filename,combModel)
  case (3)
     call chemtable_3D_init(filename,combModel)
  case (4)
     call chemtable_4D_init(filename,combModel)
  case default
     call die('chemtable_init: Wrong number of chemtable dimensions')
  end select

  ! Check that the appropriate scalars are present
  select case(trim(combModel))
  case ('Single Flamelet')
     if (isc_ZMIX.eq.0) call die('chemtable_init: Single flamelet model requires scalar ZMIX')
  case ('Steady Flamelet')
     if (isc_ZMIX.eq.0) call die('chemtable_init: Steady flamelet model requires scalar ZMIX')
  case ('FPVA')
     if ((isc_ZMIX.eq.0) .or. (isc_PROG.eq.0)) call die('chemtable_init: FPV model requires scalars ZMIX and PROG')
  case ('PFPVA')
     if ((isc_PROG.eq.0) .or. (isc_ZMIX.eq.0)) call die('chemtable_init: PFPV model requires scalars PROG and ZMIX')
  case ('UFPVA')
     if ((isc_ZMIX.eq.0) .or. (isc_PROG.eq.0)) call die('chemtable_init: Unsteady FPV model requires scalars ZMIX and PROG')
  case ('RFPVA')
     if ((isc_ZMIX.eq.0) .or. (isc_PROG.eq.0) .or. (isc_ENTH.eq.0)) call die('chemtable_init: Radiation FPV model requires scalars ZMIX, PROG, and H')
  case ('RPFPVA')
     if ((isc_PROG.eq.0) .or. (isc_ZMIX.eq.0) .or. (isc_ENTH.eq.0)) call die('chemtable_init: Radiation FPV model requires scalars ZMIX, PROG, and H')
  end select

  ! Get the density model
  call parser_read('Density model',density_model,'ARTS')

  return
end subroutine chemtable_init


! ========================= !
! Chemtable lookup routines !
! ========================= !
subroutine chemtable_lookup(tag,Ain)
  use chemtable
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: Ain

  select case(combDim)
  case (1)
     call chemtable_1D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),nxo_*nyo_*nzo_)
  case (2)
     select case(trim(combModel))
     case ('Single Flamelet')
        call chemtable_2D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),ZVAR              ,nxo_*nyo_*nzo_)
     case ('Steady Flamelet')
        call chemtable_2D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),CHI               ,nxo_*nyo_*nzo_)
     case ('FPVA')
        call chemtable_2D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),SC(:,:,:,isc_PROG),nxo_*nyo_*nzo_)
     case ('PFPVA')
        call chemtable_2D_lookup(tag,Ain,SC(:,:,:,isc_PROG),SC(:,:,:,isc_ZMIX),nxo_*nyo_*nzo_)
     end select
  case (3)
     select case(trim(combModel))
     case ('Steady Flamelet')
        call chemtable_3D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),ZVAR              ,CHI               ,nxo_*nyo_*nzo_)
     case ('FPVA')
        call chemtable_3D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),ZVAR              ,SC(:,:,:,isc_PROG),nxo_*nyo_*nzo_)
     case ('UFPVA')
        call chemtable_3D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),SC(:,:,:,isc_PROG),CHI               ,nxo_*nyo_*nzo_)
     case ('RFPVA')
        call chemtable_3D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),SC(:,:,:,isc_PROG),SC(:,:,:,isc_ENTH)   ,nxo_*nyo_*nzo_)
     case ('RPFPVA')
        call chemtable_3D_lookup(tag,Ain,SC(:,:,:,isc_PROG),SC(:,:,:,isc_ZMIX),SC(:,:,:,isc_ENTH)   ,nxo_*nyo_*nzo_)
     end select
  case (4)
     select case(trim(combModel))
     case ('UFPVA')
        call chemtable_4D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),ZVAR,SC(:,:,:,isc_PROG),CHI            ,nxo_*nyo_*nzo_)
     case ('RFPVA')
        call chemtable_4D_lookup(tag,Ain,SC(:,:,:,isc_ZMIX),ZVAR,SC(:,:,:,isc_PROG),SC(:,:,:,isc_ENTH),nxo_*nyo_*nzo_)
     end select
  end select

  return
end subroutine chemtable_lookup


! ================================================== !
! Chemtable lookup routines for single density value !
! -> For use with Newton inversion                   !
! -> Currently, only invert with mixture fraction    !
! -> Chucks on other scalars                         !
! ================================================== !
function chemtable_lookup_rho_local(i,j,k,Z1)
  use chemtable
  implicit none

  real(WP) :: chemtable_lookup_rho_local
  integer, intent(in) :: i,j,k
  real(WP), intent(in) :: Z1

  real(WP) :: chemtable_1D_lookup_rho_local, chemtable_2D_lookup_rho_local
  real(WP) :: chemtable_3D_lookup_rho_local, chemtable_4D_lookup_rho_local

  select case(combDim)
  case (1)
     chemtable_lookup_rho_local = chemtable_1D_lookup_rho_local(Z1)
  case (2)
     select case(trim(combModel))
     case ('Single Flamelet')
        chemtable_lookup_rho_local = chemtable_2D_lookup_rho_local(Z1,ZVAR(i,j,k))
     case ('Steady Flamelet')
        chemtable_lookup_rho_local = chemtable_2D_lookup_rho_local(Z1,CHI(i,j,k))
     case ('FPVA')
        chemtable_lookup_rho_local = chemtable_2D_lookup_rho_local(Z1,SC(i,j,k,isc_PROG))
     case ('PFPVA')
        chemtable_lookup_rho_local = chemtable_2D_lookup_rho_local(Z1,SC(i,j,k,isc_ZMIX))
     end select
  case (3)
     select case(trim(combModel))
     case ('Steady Flamelet')
        chemtable_lookup_rho_local = chemtable_3D_lookup_rho_local(Z1,ZVAR(i,j,k),CHI(i,j,k))
     case ('FPVA')
        chemtable_lookup_rho_local = chemtable_3D_lookup_rho_local(Z1,ZVAR(i,j,k),SC(i,j,k,isc_PROG))
     case ('UFPVA')
        chemtable_lookup_rho_local = chemtable_3D_lookup_rho_local(Z1,SC(i,j,k,isc_PROG),CHI(i,j,k))
     case ('RFPVA')
        chemtable_lookup_rho_local = chemtable_3D_lookup_rho_local(Z1,SC(i,j,k,isc_PROG),SC(i,j,k,isc_ENTH))
     case ('RPFPVA')
        chemtable_lookup_rho_local = chemtable_3D_lookup_rho_local(Z1,SC(i,j,k,isc_ZMIX),SC(i,j,k,isc_ENTH))
     end select
  case (4)
     select case(trim(combModel))
     case ('UFPVA')
        chemtable_lookup_rho_local = chemtable_4D_lookup_rho_local(Z1,ZVAR(i,j,k),SC(i,j,k,isc_PROG),CHI(i,j,k))
     case ('RFPVA')
        chemtable_lookup_rho_local = chemtable_4D_lookup_rho_local(Z1,ZVAR(i,j,k),SC(i,j,k,isc_PROG),SC(i,j,k,isc_ENTH))
     end select
  end select
     
  return
end function chemtable_lookup_rho_local


function chemtable_lookup_rho_der(dir)
  use chemtable
  implicit none

  real(WP) :: chemtable_lookup_rho_der
  integer, intent(in) :: dir

  real(WP) :: chemtable_1D_lookup_rho_der, chemtable_2D_lookup_rho_der
  real(WP) :: chemtable_3D_lookup_rho_der, chemtable_4D_lookup_rho_der

  select case(combDim)
  case (1)
     chemtable_lookup_rho_der = chemtable_1D_lookup_rho_der(dir)
  case (2)
     chemtable_lookup_rho_der = chemtable_2D_lookup_rho_der(dir)
  case (3)
     chemtable_lookup_rho_der = chemtable_3D_lookup_rho_der(dir)
  case (4)
     chemtable_lookup_rho_der = chemtable_4D_lookup_rho_der(dir)
  end select

  return
end function chemtable_lookup_rho_der


! ============================================ !
! Source terms for chemtable mapping variables !
! ============================================ !
subroutine chemtable_source(srcSC)
  use chemtable
  use memory
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: srcSC
  integer :: i,j,k

  ! Progress variable: chemistry
  if (isc_PROG.ne.0) then
     call chemtable_lookup('SRC_PROG',tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              srcSC(i,j,k,isc_PROG) = srcSC(i,j,k,isc_PROG) + dt_*RHO(i,j,k)*tmp1(i,j,k)
              ! SD force ignition
!!$              if (i.eq.101) then
!!$                 srcSC(i,j,k,isc_PROG) = srcSC(i,j,k,isc_PROG) + dt_*RHO(i,j,k)*max(0.3_WP - tmp1(i,j,k), 0.0_WP)/2e-7_WP
!!$              end if
              
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Enthalpy: radiation
  if (isc_ENTH.ne.0) then
     call chemtable_lookup('SRC_RAD',tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              srcSC(i,j,k,isc_ENTH) = srcSC(i,j,k,isc_ENTH) + dt_*RHO(i,j,k)*tmp1(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine chemtable_source


! ============================= !
! Determine density given rhoSC !
! ============================= !
subroutine chemtable_invert_density
  use chemtable
  implicit none

  integer :: i,j,k,n
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: RHOtmp

  select case (trim(density_model))
  case ('chucks non-rescaled')
     ! Get new density
     call chemtable_lookup('RHO',RHO)
  case ('chucks')
     ! Get new density
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              RHOtmp(i,j,k) = RHO(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call chemtable_lookup('RHO',RHO)
     ! Rescale
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,:) = RHOtmp(i,j,k) * SC(i,j,k,:) / RHO(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  case ('ARTS')
     ! Invert with only mixture fraction (chucks on others)
     ! Need to look at full matrix inversion again (especially for DNS)
     call newton_invert
  end select
     
  return
end subroutine chemtable_invert_density
        
 
