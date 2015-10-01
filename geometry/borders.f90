module borders
  use precision
  use geometry
  use masks
  implicit none


  ! Definition of type bc
  ! -> location (of walls)
  ! -> area of the cross section
  type bc
     integer :: jmin, jmax
     integer :: jmino, jmaxo
     real(WP) :: A
  end type bc

  ! Number of inlets and array of inlets
  integer :: ninlet
  type(bc), dimension(:), pointer :: inlet

  ! Number of outlets and array of outlets
  integer :: noutlet
  type(bc), dimension(:), pointer :: outlet

  ! Totral outlet area
  real(WP) :: outlet_area

contains

  ! ====================================== !
  ! Detect the inlets on the left boundary !
  ! ====================================== !
  subroutine borders_inlet
    implicit none
    integer :: nflow,j

    ! Allocate
    allocate(inlet(ny))

    ! Locate the inlet boundary points
    ninlet = 0
    if (mask(imino,jmin-1).eq.2) then
       ninlet = 1
       inlet(ninlet)%jmin = jmin
       inlet(ninlet)%jmax = 0
    end if

    do j=jmin,jmax     
       if ((mask(imino,j).eq.2) .and. (mask(imino,j+1).eq.1)) then
          inlet(ninlet)%jmax = j
       else if ((mask(imino,j).eq.1) .and. (mask(imino,j+1).eq.2)) then
          ninlet = ninlet + 1
          inlet(ninlet)%jmin = j+1
          inlet(ninlet)%jmax = 0
       end if
    end do

    ! If no inlet return 
    if (ninlet.eq.0) return

    if (inlet(ninlet)%jmax.eq.0) then
       inlet(ninlet)%jmax = jmax
    end if

    ! Compute cross section area of the inlets
    if (icyl.eq.0) then
       do nflow=1,ninlet
          inlet(nflow)%A = (y(inlet(nflow)%jmax+1)-y(inlet(nflow)%jmin)) * (z(kmax+1)-z(kmin))
       end do
    else
       do nflow=1,ninlet
          inlet(nflow)%A = 0.5_WP* &
               (y(inlet(nflow)%jmax+1)**2-y(inlet(nflow)%jmin)**2 ) * (z(kmax+1)-z(kmin))
       end do
    end if

    ! Get the min and max index with ghost cells
    do nflow=1,ninlet
       ! Min
       inlet(nflow)%jmino = inlet(nflow)%jmin
       loop1:do j=1,nover
          if (mask(imino,inlet(nflow)%jmin-j).eq.2) then
             inlet(nflow)%jmino = inlet(nflow)%jmin-j
          else
             exit loop1
          end if
       end do loop1
       ! Max
       inlet(nflow)%jmaxo = inlet(nflow)%jmax
       loop2:do j=1,nover
          if (mask(imino,inlet(nflow)%jmax+j).eq.2) then
             inlet(nflow)%jmaxo = inlet(nflow)%jmax+j
          else
             exit loop2
          end if
       end do loop2
    end do

    return
  end subroutine borders_inlet


  ! ======================================== !
  ! Detect the outlets on the right boundary !
  ! ======================================== !
  subroutine borders_outlet
    implicit none
    integer :: nflow,j

    ! Allocate
    allocate(outlet(ny))

    ! Locate the outlet boundary points
    noutlet = 0
    if (mask(imaxo,jmin-1).EQ.2) then
       noutlet = 1
       outlet(noutlet)%jmin = jmin
       outlet(noutlet)%jmax = 0
    end if

    do j=jmin,jmax     
       if ((mask(imaxo,j).eq.2) .and. (mask(imaxo,j+1).eq.1)) then
          outlet(noutlet)%jmax = j
       else if ((mask(imaxo,j).eq.1) .and. (mask(imaxo,j+1).eq.2)) then
          noutlet = noutlet + 1
          outlet(noutlet)%jmin = j+1
          outlet(noutlet)%jmax = 0
       end if
    end do
    
    ! If no outlet return 
    if (noutlet.eq.0) return

    if (outlet(noutlet)%jmax.eq.0) then
       outlet(noutlet)%jmax = jmax
    end if
    
    ! Compute cross section area of the outlets
    outlet_area = 0.0_WP
    if (icyl.eq.0) then
       do nflow=1,noutlet
          outlet(nflow)%A = (y(outlet(nflow)%jmax+1)-y(outlet(nflow)%jmin)) * (z(kmax+1)-z(kmin))
          outlet_area = outlet_area + outlet(nflow)%A
       end do
    else
       do nflow=1,noutlet
          outlet(nflow)%A = 0.5_WP* (y(outlet(nflow)%jmax+1)**2-y(outlet(nflow)%jmin)**2 ) * (z(kmax+1)-z(kmin))
          outlet_area = outlet_area + outlet(nflow)%A
       end do
    end if

    ! Get the min and max index with ghost cells
    do nflow=1,noutlet
       ! Min
       outlet(nflow)%jmino = outlet(nflow)%jmin
       loop1:do j=1,nover
          if (mask(imaxo,outlet(nflow)%jmin-j).eq.2) then
             outlet(nflow)%jmino = outlet(nflow)%jmin-j
          else
             exit loop1
          end if
       end do loop1
       ! Max
       outlet(nflow)%jmaxo = outlet(nflow)%jmax
       loop2:do j=1,nover
          if (mask(imaxo,outlet(nflow)%jmax+j).eq.2) then
             outlet(nflow)%jmaxo = outlet(nflow)%jmax+j
          else
             exit loop2
          end if
       end do loop2
    end do

    return
  end subroutine borders_outlet

end module borders



! ====================== !
! Initialize the borders !
! -> Detect inlets       !
! -> Detect outlets      !
! ====================== !
subroutine borders_init
  use borders
  implicit none
  
  call borders_inlet
  call borders_outlet
  
  return
end subroutine borders_init
