module chemtable_2D
  use string
  use precision
  implicit none

  ! Mapping variables of the chemtable
  integer :: n1,n2
  real(WP) :: x1min, x1max, x2min, x2max
  real(WP), dimension(:), pointer :: x1,x2

  ! Number of variables tabulated
  integer :: nvar
  character(len=str_medium), dimension(:), pointer :: chem_name

  ! Arrays of mapped variables
  real(WP), dimension(:,:,:), pointer :: table
  
  ! Store the values for interpolation for speedup in newton
  integer :: index_rho
  integer :: i1, i2
  real(WP) :: w11, w12, w21, w22

  !$OMP THREADPRIVATE(i1,i2,w11,w12,w21,w22)

end module chemtable_2D


! ===================== !
! Read in the chemtable !
! ===================== !
subroutine chemtable_2D_init(filename,combModel)
  use chemtable_2D
  use parser
  use parallel
  implicit none

  character(len=str_medium), intent(inout) :: filename, combModel
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer  :: ierr, var, iunit

  ! Open the chemtable file
  filename = trim(mpiiofs) // ":" // trim(filename)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,mpi_info,iunit,ierr)
  if (ierr .ne. 0) call die("chemtable_2D_init: error opening the chemtable")

  ! Read the headers
  call MPI_FILE_READ_ALL(iunit,n1,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,n2,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,nvar,1,MPI_INTEGER,status,ierr)

  ! Allocate the corresponding arrays
  allocate (x1(n1), x2(n2))
  allocate (chem_name(nvar))
  allocate (table(n1,n2,nvar))
  
  ! Read the mapping variables
  call MPI_FILE_READ_ALL(iunit,x1,n1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(iunit,x2,n2,MPI_REAL_WP,status,ierr)

  ! Read the combustion model used
  call MPI_FILE_READ_ALL(iunit,combModel,str_medium,MPI_CHARACTER,status,ierr)

  ! Read the names of the mapped variables
  do var=1,nvar
     call MPI_FILE_READ_ALL(iunit,chem_name(var),str_medium,MPI_CHARACTER,status,ierr)
  end do

  ! Read the mapped variables
  do var=1,nvar
     call MPI_FILE_READ_ALL(iunit,table(:,:,var),n1*n2,MPI_REAL_WP,status,ierr)
     ! Store index for density
     if (trim(chem_name(var)).eq.'RHO') index_rho = var
  end do
  
  ! Get some properties of the mapping
  x1min = minval(x1)
  x1max = maxval(x1)
  x2min = minval(x2)
  x2max = maxval(x2)

  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)

  return
end subroutine chemtable_2D_init


! ====================================================== !
! Look in the chemtable for the variable named 'tag'     !
! with the value A1 and A2 for the two mapping variables !
! ====================================================== !
subroutine chemtable_2D_lookup(tag, R, A1, A2, n)
  use chemtable_2D
  implicit none

  integer, intent(in) :: n
  character(len=*), intent(in) :: tag
  real(WP), dimension(n), intent(in)  :: A1, A2
  real(WP), dimension(n), intent(out) :: R

  integer :: var, i, j
  
  ! If density call another routine
  if (trim(tag).eq.'RHO') then
     call chemtable_2D_lookup_rho(R,A1,A2,n)
     return
  end if
  
  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_2D_lookup: unknown variable : '//trim(tag))
  end if

  ! Bilinear interpolation
  !$OMP PARALLEL DO
  do i=1,n

     ! First direction
     if (A1(i).lt.x1min) then
        i1 = 1
        w11 = 1.0_WP
     else if (A1(i).ge.x1max) then
        i1 = n1-1
        w11 = 0.0_WP
     else
        loop1:do j=1,n1-1
           if (A1(i).lt.x1(j+1)) then
              i1 = j
              exit loop1
           end if
        end do loop1
        w11 = (x1(i1+1)-A1(i))/(x1(i1+1)-x1(i1))
     end if
     w12 = 1.0_WP - w11

     ! Second direction
     if (A2(i).lt.x2min) then
        i2 = 1
        w21 = 1.0_WP
     else if (A2(i).ge.x2max) then
        i2 = n2-1
        w21 = 0.0_WP
     else
        loop2:do j=1,n2-1
           if (A2(i).lt.x2(j+1)) then
              i2 = j
              exit loop2
           end if
        end do loop2
        w21 = (x2(i2+1)-A2(i))/(x2(i2+1)-x2(i2))
     end if
     w22 = 1.0_WP - w21

     ! Interpolation
     R(i) =  w21*( w11*table(i1  ,i2  ,var)     &
                  +w12*table(i1+1,i2  ,var) )   &
            +w22*( w11*table(i1  ,i2+1,var)     &
                  +w12*table(i1+1,i2+1,var) )
  end do
  !$OMP END PARALLEL DO

  return
end subroutine chemtable_2D_lookup


! ====================================================== !
! Look in the chemtable for the density                  !
! with the value A1 and A2 for the two mapping variables !
! Different interpolation than for other variables       !
! ====================================================== !
subroutine chemtable_2D_lookup_rho(R, A1, A2, n)
  use chemtable_2D
  implicit none
  
  integer, intent(in) :: n
  real(WP), dimension(n), intent(in)  :: A1, A2
  real(WP), dimension(n), intent(out) :: R

  integer :: i,j

  ! Bilinear interpolation
  !$OMP PARALLEL DO
  do i=1,n

     ! First direction
     if (A1(i).lt.x1min) then
        i1 = 1
        w11 = 1.0_WP
     else if (A1(i).ge.x1max) then
        i1 = n1-1
        w11 = 0.0_WP
     else
        loop1:do j=1,n1-1
           if (A1(i).lt.x1(j+1)) then
              i1 = j
              exit loop1
           end if
        end do loop1
        w11 = (x1(i1+1)-A1(i))/(x1(i1+1)-x1(i1))
     end if
     w12 = 1.0_WP - w11

     ! Second direction
     if (A2(i).lt.x2min) then
        i2 = 1
        w21 = 1.0_WP
     else if (A2(i).ge.x2max) then
        i2 = n2-1
        w21 = 0.0_WP
     else
        loop2:do j=1,n2-1
           if (A2(i).lt.x2(j+1)) then
              i2 = j
              exit loop2
           end if
        end do loop2
        w21 = (x2(i2+1)-A2(i))/(x2(i2+1)-x2(i2))
     end if
     w22 = 1.0_WP - w21

     ! Interpolation of 1/rho
     R(i) =  w21*( w11/table(i1  ,i2  ,index_rho)   &
                  +w12/table(i1+1,i2  ,index_rho) ) &
            +w22*( w11/table(i1  ,i2+1,index_rho)   &
                  +w12/table(i1+1,i2+1,index_rho) )
     R(i) = 1.0_WP/R(i)
  end do
  !$OMP END PARALLEL DO

  return
end subroutine chemtable_2D_lookup_rho


! ====================================================== !
! Look in the chemtable for the density                  !
! with the value A1 and A2 for the two mapping variables !
! ====================================================== !
subroutine chemtable_2D_lookup_local(tag, R, A1, A2)
  use chemtable_2D
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), intent(out) :: R
  real(WP), intent(in)  :: A1, A2

  integer :: var, j

  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_2D_lookup: unknown variable : '//trim(tag))
  end if

  ! First direction
  if (A1.lt.x1min) then
     i1 = 1
     w11 = 1.0_WP
  else if (A1.ge.x1max) then
     i1 = n1-1
     w11 = 0.0_WP
  else
     loop1:do j=1,n1-1
        if (A1.lt.x1(j+1)) then
           i1 = j
           exit loop1
        end if
     end do loop1
     w11 = (x1(i1+1)-A1)/(x1(i1+1)-x1(i1))
  end if
  w12 = 1.0_WP - w11
  
  ! Second direction
  if (A2.lt.x2min) then
     i2 = 1
     w21 = 1.0_WP
  else if (A2.ge.x2max) then
     i2 = n2-1
     w21 = 0.0_WP
  else
     loop2:do j=1,n2-1
        if (A2.lt.x2(j+1)) then
           i2 = j
           exit loop2
        end if
     end do loop2
     w21 = (x2(i2+1)-A2)/(x2(i2+1)-x2(i2))
  end if
  w22 = 1.0_WP - w21
  
  ! Interpolation
  R =  w21*( w11*table(i1  ,i2  ,var)   &
            +w12*table(i1+1,i2  ,var) ) &
      +w22*( w11*table(i1  ,i2+1,var)   &
            +w12*table(i1+1,i2+1,var) )

  return
end subroutine chemtable_2D_lookup_local


! ====================================================== !
! Look in the chemtable for the density                  !
! with the value A1 and A2 for the two mapping variables !
! ====================================================== !
function chemtable_2D_lookup_rho_local(A1, A2)
  use chemtable_2D
  implicit none

  real(WP) :: chemtable_2D_lookup_rho_local
  real(WP), intent(in)  :: A1, A2

  integer :: var, j

  ! The variable is rho
  var = index_rho

  ! First direction
  if (A1.lt.x1min) then
     i1 = 1
     w11 = 1.0_WP
  else if (A1.ge.x1max) then
     i1 = n1-1
     w11 = 0.0_WP
  else
     loop1:do j=1,n1-1
        if (A1.lt.x1(j+1)) then
           i1 = j
           exit loop1
        end if
     end do loop1
     w11 = (x1(i1+1)-A1)/(x1(i1+1)-x1(i1))
  end if
  w12 = 1.0_WP - w11
  
  ! Second direction
  if (A2.lt.x2min) then
     i2 = 1
     w21 = 1.0_WP
  else if (A2.ge.x2max) then
     i2 = n2-1
     w21 = 0.0_WP
  else
     loop2:do j=1,n2-1
        if (A2.lt.x2(j+1)) then
           i2 = j
           exit loop2
        end if
     end do loop2
     w21 = (x2(i2+1)-A2)/(x2(i2+1)-x2(i2))
  end if
  w22 = 1.0_WP - w21
  
  ! Interpolation
  chemtable_2D_lookup_rho_local =  1.0_WP / ( &
        w21*( w11/table(i1  ,i2  ,var)   &
             +w12/table(i1+1,i2  ,var) ) &
       +w22*( w11/table(i1  ,i2+1,var)   &
             +w12/table(i1+1,i2+1,var) ) )

  return
end function chemtable_2D_lookup_rho_local


! ======================================================= !
! Look in the chemtable for the derivative of the density !
! with the value A1 and A2 for the two mapping variables  !
! ======================================================= !
function chemtable_2D_lookup_rho_der(dir)
  use chemtable_2D
  implicit none

  real(WP) :: chemtable_2D_lookup_rho_der
  integer, intent(in) :: dir

  real(WP) :: c11, c12, c21, c22
  integer :: var

  ! Get the index of the variable
  var = index_rho

  ! Compute the coefficients
  c11 = w11
  c12 = w12
  c21 = w21
  c22 = w22

  ! Update the coefficients to account for derivatives
  select case(dir)
  case (1)
     c12 = 1.0_WP/(x1(i1+1)-x1(i1))
     c11 = -c12
  case (2)
     c22 = 1.0_WP/(x2(i2+1)-x2(i2))
     c21 = -c22
  end select

  ! Compute derivative in the direction 'dir' 
  ! and interpolation in the two others
  chemtable_2D_lookup_rho_der =  &
        c21*( c11*table(i1  ,i2  ,var)   &
             +c12*table(i1+1,i2  ,var) ) &
       +c22*( c11*table(i1  ,i2+1,var)   &
             +c12*table(i1+1,i2+1,var) )
  return
end function chemtable_2D_lookup_rho_der


! =============================================== !
! Find the maximum of a variable in the chemtable !
! =============================================== !
subroutine chemtable_2D_lookup_max(tag, R)
  use chemtable_2D
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), intent(out) :: R

  integer :: var

  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_2D_lookup: unknown variable : '//trim(tag))
  end if

  ! Return the max
  R = maxval(table(:,:,var))

  return
end subroutine chemtable_2D_lookup_max

! =============================================== !
! Find the minimum of a variable in the chemtable !
! =============================================== !
subroutine chemtable_2D_lookup_min(tag, R)
  use chemtable_2D
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), intent(out) :: R

  integer :: var

  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_2D_lookup: unknown variable : '//trim(tag))
  end if

  ! Return the max
  R = minval(table(:,:,var))

  return
end subroutine chemtable_2D_lookup_min

