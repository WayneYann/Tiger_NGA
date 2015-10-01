module newton
  use precision
  use combustion
  implicit none

  integer, parameter :: nniter = 10
  real(WP), external :: chemtable_lookup_rho_local,chemtable_lookup_rho_der

contains
  
  ! Invert a linear system of 'size' eqs
  subroutine invert(rhs,jac,size)
    implicit none
    
    integer, intent(in) :: size
    real(WP), dimension(size) :: rhs
    real(WP), dimension(size,size) :: jac
    real(WP), dimension(size) :: sol

    select case(size)
    case (1)
       rhs = rhs(1) / jac(1,1)
    case (2)
       rhs = rhs / (jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)) 
       sol(1) = + jac(2,2)*rhs(1) - jac(1,2)*rhs(2)
       sol(2) = - jac(2,1)*rhs(1) + jac(1,1)*rhs(2)
       rhs = sol
    case (3)
       rhs(1) = rhs(1) / jac(1,1)
       jac(1,:) = jac(1,:) / jac(1,1)
       rhs(2) = rhs(2) - rhs(1)*jac(2,1)
       jac(2,:) = jac(2,:) - jac(1,:)*jac(2,1)
       rhs(3) = rhs(3) - rhs(1)*jac(3,1)
       jac(3,:) = jac(3,:) - jac(1,:)*jac(3,1)

       rhs(2) = rhs(2) / jac(2,2)
       jac(2,:) = jac(2,:) / jac(2,2)

       rhs(3) = rhs(3) - rhs(2)*jac(3,2)
       jac(3,:) = jac(3,:) - jac(2,:)*jac(3,2)

       sol(3) = rhs(3)/jac(3,3)
       sol(2) = rhs(2) - jac(2,3)*sol(3)
       sol(1) = rhs(1) - jac(1,2)*sol(2) - jac(1,3)*sol(3)

       rhs = sol
    case default
       call die('newton/invert: not yet implemented')
    end select

    return
  end subroutine invert

  ! ---------------------------------------------------------------------

  ! Invert the 1D nonlinear system
  ! -> rho increases with ZMIX
  subroutine invert_ZMIX_inc(i,j,k,rho1,der)
    implicit none

    real(WP), intent(out) :: rho1,der
    integer, intent(in) :: i,j,k
    integer :: n
    real(WP) :: rhs(1),jac(1,1)
    real(WP) :: rhoZ,Z1

    ! Value we want to invert
    rhoZ = SC(i,j,k,isc_ZMIX)*RHO(i,j,k)
    
    ! Initial guess = old value
    Z1 = SC(i,j,k,isc_ZMIX)
    
    ! Clip the values
    Z1 = min(1.0_WP,max(0.0_WP,Z1))
    
    ! Get the guess density and the derivative
    rho1 = chemtable_lookup_rho_local(i,j,k,Z1)
    
    ! Newton method
    do n=1,nniter
       rhs(1) = rho1*Z1 - rhoZ
       jac(1,1) = rho1 + Z1*chemtable_lookup_rho_der(1)
       
       call invert(rhs,jac,1)
       Z1 = Z1 - rhs(1)
       
       Z1 = min(1.0_WP,max(0.0_WP,Z1))
       rho1 = chemtable_lookup_rho_local(i,j,k,Z1)
    end do

    ! Return the derivative
    der = chemtable_lookup_rho_der(1)
    
    return
  end subroutine invert_ZMIX_inc

  ! Invert the 1D nonlinear system
  ! -> rho decreases with ZMIX
  subroutine invert_ZMIX_dec(i,j,k,rho1,der)
    implicit none

    real(WP), intent(out) :: rho1,der
    integer, intent(in) :: i,j,k
    integer :: n
    real(WP) :: rhs(1),jac(1,1)
    real(WP) :: rhoY,Y1

    ! Value we want to invert
    rhoY = (1.0_WP-SC(i,j,k,isc_ZMIX))*RHO(i,j,k)
    
    ! Initial guess = old value
    Y1 = 1.0_WP - SC(i,j,k,isc_ZMIX)
    
    ! Clip the values
    Y1 = min(1.0_WP,max(0.0_WP,Y1))
    
    ! Get the guess density and the derivative
    rho1 = chemtable_lookup_rho_local(i,j,k,1.0_WP-Y1)
    
    ! Newton method
    do n=1,nniter
       rhs(1) = rho1*Y1 - rhoY
       jac(1,1) = rho1 - Y1*chemtable_lookup_rho_der(1)
       
       call invert(rhs,jac,1)
       Y1 = Y1 - rhs(1)
       
       Y1 = min(1.0_WP,max(0.0_WP,Y1))
       rho1 = chemtable_lookup_rho_local(i,j,k,1.0_WP-Y1)
    end do
    
    ! Return the derivative
    der = chemtable_lookup_rho_der(1)

    return
  end subroutine invert_ZMIX_dec
  
  
  ! Invert the MD nonlinear system

   
end module newton


! =============================================== !
! Newton inversion                                !
! -> Currently, only invert with mixture fraction !
! -> Chucks on other scalars                      !
! =============================================== !
subroutine newton_invert
  use newton
  implicit none

  integer :: i,j,k
  real(WP) :: rho_old,rho_new,der
  real(WP) :: Z1

  !$OMP PARALLEL DO PRIVATE(Z1,rho_old,rho_new,der)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           ! Initial guess = old value
           Z1 = SC(i,j,k,isc_ZMIX)
           
           ! Clip the values
           Z1 = min(1.0_WP,max(0.0_WP,Z1))
           
           ! Get the guess density and save old density
           rho_old = RHO(i,j,k)
           rho_new = chemtable_lookup_rho_local(i,j,k,Z1)
           der     = chemtable_lookup_rho_der(1)
           
           ! Apply newton on the right quantity
           if (der.gt.0.0_WP) then
              call invert_ZMIX_inc(i,j,k,rho_new,der)
              if (der.lt.0.0_WP) call invert_ZMIX_dec(i,j,k,rho_new,der)
           else
              call invert_ZMIX_dec(i,j,k,rho_new,der)
              if (der.gt.0.0_WP) call invert_ZMIX_inc(i,j,k,rho_new,der)
           end if
           
           ! Now get all the scalars from conservation of rho*sc
           RHO(i,j,k) = rho_new
           SC(i,j,k,:) = rho_old * SC(i,j,k,:) / rho_new
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine newton_invert
