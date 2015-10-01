module centerline
  use precision
  use geometry
  use partition
  use parallel
  implicit none

end module centerline

! ========================================= !
! Update any quantity across the centerline !
! Assumes that the values are cell centered !
! Assume there is no decomposition in z     !
! ========================================= !
subroutine centerline_update_ym(A,coeff)
  use centerline
  use data
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), intent(in) :: coeff
  integer :: i,j,k,k0,kPi
  
  ! Nothing to do
  if (icyl.ne.1 .or. jproc.ne.1) return

  ! Update the opposite points
  !$OMP PARALLEL DO PRIVATE(k0,kPi)
  do k=kmino_,kmaxo_
     k0 = kmin+modulo(k-kmin,nz)
     kPi = k0+nz/2
     if (kPi.gt.kmax) kPi = kPi-nz
     
     do i=imino_,imaxo_
        do j=1,nover
           A(i,jmin-j,k) = coeff*A(i,jmin+j-1,kPi)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine centerline_update_ym


! ============ !
! FOR INTEGERS !
! ============ !
subroutine centerline_update_ym_int(A,coeff)
  use centerline
  use data
  implicit none
  
  integer, dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  integer, intent(in) :: coeff
  integer :: i,j,k,k0,kPi
  
  ! Nothing to do
  if (icyl.ne.1 .or. jproc.ne.1) return

  ! Update the opposite points
  !$OMP PARALLEL DO PRIVATE(k0,kPi)
  do k=kmino_,kmaxo_
     k0 = kmin+modulo(k-kmin,nz)
     kPi = k0+nz/2
     if (kPi.gt.kmax) kPi = kPi-nz
        
     do i=imino_,imaxo_
        do j=1,nover
           A(i,jmin-j,k) = coeff*A(i,jmin+j-1,kPi)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine centerline_update_ym_int


! ========================================= !
! FOR V FOR V FOR V FOR V FOR V FOR V FOR V !
! Update only V/rhoV across the centerline  !
! Assume there is no decomposition in z     !
! ========================================= !
subroutine centerline_update_y(A,coeff)
  use centerline
  use data
  use masks
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), intent(in) :: coeff
  integer  :: i,j,k,k0,kPi
  real(WP) :: Acos,Asin
  
  ! Nothing to do
  if (icyl.ne.1 .or. jproc.ne.1) return
  
  ! Update the opposite points
  !$OMP PARALLEL DO PRIVATE(k0,kPi)
  do k=kmino_,kmaxo_
     k0 = kmin+modulo(k-kmin,nz)
     kPi = k0+nz/2
     if (kPi.gt.kmax) kPi = kPi-nz
     
     do i=imino_,imaxo_
        do j=1,nover
           A(i,jmin-j,k) = coeff*A(i,jmin+j,kPi)       
        end do
     end do
     
  end do
  !$OMP END PARALLEL DO
  
  ! Points on the axis
  if (coeff.lt.0.0_WP) then
     !$OMP PARALLEL DO PRIVATE(Acos,Asin)
     do i=imino_,imaxo_
        if (mask(i,jmin).ne.1) then
           Acos = sum((A(i,jmin-1,kmin:kmax)+A(i,jmin+1,kmin:kmax))*cos(zm(kmin:kmax)))/real(nz,WP)
           Asin = sum((A(i,jmin-1,kmin:kmax)+A(i,jmin+1,kmin:kmax))*sin(zm(kmin:kmax)))/real(nz,WP)
           A(i,jmin,kmino:kmaxo) = Acos*cos(zm(kmino:kmaxo)) + Asin*sin(zm(kmino:kmaxo))
        end if
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO
     do i=imino_,imaxo_
        if (mask(i,jmin).ne.1) then
           A(i,jmin,:) = sum(A(i,jmin,kmin:kmax))/real(nz,WP)
        end if
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine centerline_update_y


! ============ !
! FOR INTEGERS !
! ============ !
subroutine centerline_update_y_int(A,coeff)
  use centerline
  use data
  implicit none
  
  integer, dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  integer, intent(in) :: coeff
  integer :: i,j,k,k0,kPi
  
  ! Nothing to do
  if (icyl.ne.1 .or. jproc.ne.1) return
  
  ! Update the opposite points
  !$OMP PARALLEL DO PRIVATE(k0,kPi)
  do k=kmino_,kmaxo_
     k0 = kmin+modulo(k-kmin,nz)
     kPi = k0+nz/2
     if (kPi.gt.kmax) kPi = kPi-nz
     
     do i=imino_,imaxo_
        do j=1,nover
           A(i,jmin-j,k) = coeff*A(i,jmin+j,kPi)       
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine centerline_update_y_int

