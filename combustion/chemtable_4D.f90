module chemtable_4D
  use string
  use precision
  implicit none

  ! Mapping variables of the chemtable
  integer :: n1,n2,n3,n4
  real(WP) :: x1min, x1max, x2min, x2max, x3min, x3max, x4min, x4max
  real(WP), dimension(:), pointer :: x1,x2,x3,x4

  ! Number of variables tabulated
  integer :: nvar
  character(len=str_medium), dimension(:), pointer :: chem_name

  ! Arrays of mapped variables
  real(WP), dimension(:,:,:,:,:), pointer :: table
  
  ! Store the values for interpolation for speedup in newton
  integer :: index_rho
  integer :: i1, i2, i3, i4
  real(WP) :: w11, w12, w21, w22, w31, w32, w41, w42

  !$OMP THREADPRIVATE(i1,i2,i3,i4,w11,w12,w21,w22,w31,w32,w41,w42)

end module chemtable_4D


! ===================== !
! Read in the chemtable !
! ===================== !
subroutine chemtable_4D_init(filename,combModel)
  use chemtable_4D
  use parser
  use parallel
  implicit none

  character(len=str_medium), intent(inout) :: filename, combModel
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer  :: ierr, var, iunit

  ! Open the chemtable file
  filename = trim(filename)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,mpi_info,iunit,ierr)
  if (ierr .ne. 0) call die("chemtable_4D_init: error opening the chemtable")

  ! Read the headers
  call MPI_FILE_READ_ALL(iunit,n1,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,n2,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,n3,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,n4,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,nvar,1,MPI_INTEGER,status,ierr)

  ! Allocate the corresponding arrays
  allocate (x1(n1),x2(n2),x3(n3),x4(n4))
  allocate (chem_name(nvar))
  allocate (table(n1,n2,n3,n4,nvar))
  
  ! Read the mapping variables
  call MPI_FILE_READ_ALL(iunit,x1,n1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(iunit,x2,n2,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(iunit,x3,n3,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(iunit,x4,n4,MPI_REAL_WP,status,ierr)

  ! Read the combustion model used
  call MPI_FILE_READ_ALL(iunit,combModel,str_medium,MPI_CHARACTER,status,ierr)

  ! Read the names of the mapped variables
  do var=1,nvar
     call MPI_FILE_READ_ALL(iunit,chem_name(var),str_medium,MPI_CHARACTER,status,ierr)
  end do

  ! Read the mapped variables
  do var=1,nvar
     call MPI_FILE_READ_ALL(iunit,table(:,:,:,:,var),n1*n2*n3*n4,MPI_REAL_WP,status,ierr)
     ! Store index for density
     if (trim(chem_name(var)).eq.'RHO') index_rho = var
  end do
  
  ! Get some properties of the mapping
  x1min = minval(x1)
  x1max = maxval(x1)
  x2min = minval(x2)
  x2max = maxval(x2)
  x3min = minval(x3)
  x3max = maxval(x3)
  x4min = minval(x4)
  x4max = maxval(x4)

  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)

  return
end subroutine chemtable_4D_init


! ================================================================ !
! Look in the chemtable for the variable named 'tag'               !
! with the value A1, A2, A3, and A4 for the four mapping variables !
! ================================================================ !
subroutine chemtable_4D_lookup(tag, R, A1, A2, A3, A4, n)
  use chemtable_4D
  implicit none

  integer, intent(in) :: n
  character(len=*), intent(in) :: tag
  real(WP), dimension(n), intent(in)  :: A1, A2, A3, A4
  real(WP), dimension(n), intent(out) :: R

  integer :: var, i, j
  
  ! If density call another routine
  if (trim(tag).eq.'RHO') then
     call chemtable_4D_lookup_rho(R,A1,A2,A3,A4,n)
     return
  end if
  
  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_4D_lookup: unknown variable : '//trim(tag))
  end if

  ! Quadrilinear interpolation
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

     ! Third direction
     if (A3(i).lt.x3min) then
        i3 = 1
        w31 = 1.0_WP
     else if (A3(i).ge.x3max) then
        i3 = n3-1
        w31 = 0.0_WP
     else
        loop3:do j=1,n3-1
           if (A3(i).lt.x3(j+1)) then
              i3 = j
              exit loop3
           end if
        end do loop3
        w31 = (x3(i3+1)-A3(i))/(x3(i3+1)-x3(i3))
     end if
     w32 = 1.0_WP - w31

     ! Fourth direction
     if (A4(i).lt.x4min) then
        i4 = 1
        w41 = 1.0_WP
     else if (A4(i).ge.x4max) then
        i4 = n4-1
        w41 = 0.0_WP
     else
        loop4:do j=1,n4-1
           if (A4(i).lt.x4(j+1)) then
              i4 = j
              exit loop4
           end if
        end do loop4
        w41 = (x4(i4+1)-A4(i))/(x4(i4+1)-x4(i4))
     end if
     w42 = 1.0_WP - w41

     ! Interpolation
     R(i) =  w41* (w31*( w21*(  w11*table(i1  ,i2  ,i3  ,i4  ,var)       &
                              + w12*table(i1+1,i2  ,i3  ,i4  ,var) )     &
                        +w22*(  w11*table(i1  ,i2+1,i3  ,i4  ,var)       &
                              + w12*table(i1+1,i2+1,i3  ,i4  ,var) ) )   &
                  +w32*( w21*(  w11*table(i1  ,i2  ,i3+1,i4  ,var)       &
                              + w12*table(i1+1,i2  ,i3+1,i4  ,var) )     &
                        +w22*(  w11*table(i1  ,i2+1,i3+1,i4  ,var)       &
                              + w12*table(i1+1,i2+1,i3+1,i4  ,var) ) ) ) &
           + w42* (w31*( w21*(  w11*table(i1  ,i2  ,i3  ,i4+1,var)       &
                              + w12*table(i1+1,i2  ,i3  ,i4+1,var) )     &
                        +w22*(  w11*table(i1  ,i2+1,i3  ,i4+1,var)       &
                              + w12*table(i1+1,i2+1,i3  ,i4+1,var) ) )   &
                 +w32*( w21*(   w11*table(i1  ,i2  ,i3+1,i4+1,var)       &
                              + w12*table(i1+1,i2  ,i3+1,i4+1,var) )     &
                        +w22*(  w11*table(i1  ,i2+1,i3+1,i4+1,var)       &
                              + w12*table(i1+1,i2+1,i3+1,i4+1,var) ) ) )
  end do
  !$OMP END PARALLEL DO

  return
end subroutine chemtable_4D_lookup


! ================================================================ !
! Look in the chemtable for the density                            !
! with the value A1, A2, A3, and A4 for the four mapping variables !
! Different interpolation than for other variables                 !
! ================================================================ !
subroutine chemtable_4D_lookup_rho(R, A1, A2, A3, A4, n)
  use chemtable_4D
  implicit none
  
  integer, intent(in) :: n
  real(WP), dimension(n), intent(in)  :: A1, A2, A3, A4
  real(WP), dimension(n), intent(out) :: R

  integer :: i,j

  ! Quadrilinear interpolation
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

     ! Third direction
     if (A3(i).lt.x3min) then
        i3 = 1
        w31 = 1.0_WP
     else if (A3(i).ge.x3max) then
        i3 = n3-1
        w31 = 0.0_WP
     else
        loop3:do j=1,n3-1
           if (A3(i).lt.x3(j+1)) then
              i3 = j
              exit loop3
           end if
        end do loop3
        w31 = (x3(i3+1)-A3(i))/(x3(i3+1)-x3(i3))
     end if
     w32 = 1.0_WP - w31

     ! Fourth direction
     if (A4(i).lt.x4min) then
        i4 = 1
        w41 = 1.0_WP
     else if (A4(i).ge.x4max) then
        i4 = n4-1
        w41 = 0.0_WP
     else
        loop4:do j=1,n4-1
           if (A4(i).lt.x4(j+1)) then
              i4 = j
              exit loop4
           end if
        end do loop4
        w41 = (x4(i4+1)-A4(i))/(x4(i4+1)-x4(i4))
     end if
     w42 = 1.0_WP - w41

     ! Interpolation of 1/rho
     R(i) =  w41*(w31*( w21*( w11/table(i1  ,i2  ,i3  ,i4  ,index_rho)       &
                             +w12/table(i1+1,i2  ,i3  ,i4  ,index_rho) )     &
                       +w22*( w11/table(i1  ,i2+1,i3  ,i4  ,index_rho)       &
                             +w12/table(i1+1,i2+1,i3  ,i4  ,index_rho) ) )   &
                 +w32*( w21*( w11/table(i1  ,i2  ,i3+1,i4  ,index_rho)       &
                             +w12/table(i1+1,i2  ,i3+1,i4  ,index_rho) )     &
                       +w22*( w11/table(i1  ,i2+1,i3+1,i4  ,index_rho)       &
                             +w12/table(i1+1,i2+1,i3+1,i4  ,index_rho) ) ) ) &
            +w42*(w31*( w21*( w11/table(i1  ,i2  ,i3  ,i4+1,index_rho)       &
                             +w12/table(i1+1,i2  ,i3  ,i4+1,index_rho) )     &
                       +w22*( w11/table(i1  ,i2+1,i3  ,i4+1,index_rho)       &
                             +w12/table(i1+1,i2+1,i3  ,i4+1,index_rho) ) )   &
                 +w32*( w21*( w11/table(i1  ,i2  ,i3+1,i4+1,index_rho)       &
                             +w12/table(i1+1,i2  ,i3+1,i4+1,index_rho) )     &
                       +w22*( w11/table(i1  ,i2+1,i3+1,i4+1,index_rho)       &
                             +w12/table(i1+1,i2+1,i3+1,i4+1,index_rho) ) ) )    
     R(i) = 1.0_WP/R(i)
  end do
  !$OMP END PARALLEL DO

  return
end subroutine chemtable_4D_lookup_rho


! ================================================================ !
! Look in the chemtable for the density                            !
! with the value A1, A2, A3, and A4 for the four mapping variables !
! ================================================================ !
subroutine chemtable_4D_lookup_local(tag, R, A1, A2, A3, A4)
  use chemtable_4D
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), intent(out) :: R
  real(WP), intent(in)  :: A1, A2, A3, A4

  integer :: var, j

  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_4D_lookup: unknown variable : '//trim(tag))
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
  
  ! Third direction
  if (A3.lt.x3min) then
     i3 = 1
     w31 = 1.0_WP
  else if (A3.ge.x3max) then
     i3 = n3-1
     w31 = 0.0_WP
  else
     loop3:do j=1,n3-1
        if (A3.lt.x3(j+1)) then
           i3 = j
           exit loop3
        end if
     end do loop3
     w31 = (x3(i3+1)-A3)/(x3(i3+1)-x3(i3))
  end if
  w32 = 1.0_WP - w31


  ! Fourth direction
  if (A4.lt.x4min) then
     i4 = 1
     w41 = 1.0_WP
  else if (A4.ge.x4max) then
     i4 = n4-1
     w41 = 0.0_WP
  else
     loop4:do j=1,n4-1
        if (A4.lt.x4(j+1)) then
           i4 = j
           exit loop4
        end if
     end do loop4
     w41 = (x4(i4+1)-A4)/(x4(i4+1)-x4(i4))
  end if
  w42 = 1.0_WP - w41

  
  ! Interpolation
  R =  w41*(w31*( w21*( w11*table(i1  ,i2  ,i3  ,i4  ,var)       &
                       +w12*table(i1+1,i2  ,i3  ,i4  ,var) )     &
                 +w22*( w11*table(i1  ,i2+1,i3  ,i4  ,var)       &
                       +w12*table(i1+1,i2+1,i3  ,i4  ,var) ) )   &
           +w32*( w21*( w11*table(i1  ,i2  ,i3+1,i4  ,var)       &
                       +w12*table(i1+1,i2  ,i3+1,i4  ,var) )     &
                 +w22*( w11*table(i1  ,i2+1,i3+1,i4  ,var)       &
                       +w12*table(i1+1,i2+1,i3+1,i4  ,var) ) ) ) &
      +w42*(w31*( w21*( w11*table(i1  ,i2  ,i3  ,i4+1,var)       &
                       +w12*table(i1+1,i2  ,i3  ,i4+1,var) )     &
                 +w22*( w11*table(i1  ,i2+1,i3  ,i4+1,var)       &
                       +w12*table(i1+1,i2+1,i3  ,i4+1,var) ) )   &
           +w32*( w21*( w11*table(i1  ,i2  ,i3+1,i4+1,var)       &
                       +w12*table(i1+1,i2  ,i3+1,i4+1,var) )     &
                 +w22*( w11*table(i1  ,i2+1,i3+1,i4+1,var)       &
                       +w12*table(i1+1,i2+1,i3+1,i4+1,var) ) ) )

  return
end subroutine chemtable_4D_lookup_local


! ================================================================ !
! Look in the chemtable for the density                            !
! with the value A1, A2, A3, and A4 for the four mapping variables !
! ================================================================ !
function chemtable_4D_lookup_rho_local(A1, A2, A3, A4)
  use chemtable_4D
  implicit none

  real(WP) :: chemtable_4D_lookup_rho_local
  real(WP), intent(in)  :: A1, A2, A3, A4

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
  
  ! Third direction
  if (A3.lt.x3min) then
     i3 = 1
     w31 = 1.0_WP
  else if (A3.ge.x3max) then
     i3 = n3-1
     w31 = 0.0_WP
  else
     loop3:do j=1,n3-1
        if (A3.lt.x3(j+1)) then
           i3 = j
           exit loop3
        end if
     end do loop3
     w31 = (x3(i3+1)-A3)/(x3(i3+1)-x3(i3))
  end if
  w32 = 1.0_WP - w31


  ! Fourth direction
  if (A4.lt.x4min) then
     i4 = 1
     w41 = 1.0_WP
  else if (A4.ge.x4max) then
     i4 = n4-1
     w41 = 0.0_WP
  else
     loop4:do j=1,n4-1
        if (A4.lt.x4(j+1)) then
           i4 = j
           exit loop4
        end if
     end do loop4
     w41 = (x4(i4+1)-A4)/(x4(i4+1)-x4(i4))
  end if
  w42 = 1.0_WP - w41

  
  ! Interpolation
  chemtable_4D_lookup_rho_local =  1.0_WP / ( &
       w41*(w31*( w21*( w11/table(i1  ,i2  ,i3  ,i4  ,var)       &
                       +w12/table(i1+1,i2  ,i3  ,i4  ,var) )     &
                 +w22*( w11/table(i1  ,i2+1,i3  ,i4  ,var)       &
                       +w12/table(i1+1,i2+1,i3  ,i4  ,var) ) )   &
           +w32*( w21*( w11/table(i1  ,i2  ,i3+1,i4  ,var)       &
                       +w12/table(i1+1,i2  ,i3+1,i4  ,var) )     &
                 +w22*( w11/table(i1  ,i2+1,i3+1,i4  ,var)       &
                       +w12/table(i1+1,i2+1,i3+1,i4  ,var) ) ) ) &
      +w42*(w31*( w21*( w11/table(i1  ,i2  ,i3  ,i4+1,var)       &
                       +w12/table(i1+1,i2  ,i3  ,i4+1,var) )     &
                 +w22*( w11/table(i1  ,i2+1,i3  ,i4+1,var)       &
                       +w12/table(i1+1,i2+1,i3  ,i4+1,var) ) )   &
           +w32*( w21*( w11/table(i1  ,i2  ,i3+1,i4+1,var)       &
                       +w12/table(i1+1,i2  ,i3+1,i4+1,var) )     &
                 +w22*( w11/table(i1  ,i2+1,i3+1,i4+1,var)       &
                       +w12/table(i1+1,i2+1,i3+1,i4+1,var) ) ) ) )

  return
end function chemtable_4D_lookup_rho_local


! ================================================================ !
! Look in the chemtable for the derivative of the density          !
! with the value A1, A2, A3, and A4 for the four mapping variables !
! ================================================================ !
function chemtable_4D_lookup_rho_der(dir)
  use chemtable_4D
  implicit none

  real(WP) :: chemtable_4D_lookup_rho_der
  integer, intent(in) :: dir

  real(WP) :: c11, c12, c21, c22, c31, c32, c41, c42
  integer :: var

  ! Get the index of the variable
  var = index_rho

  ! Compute the coefficients
  c11 = w11
  c12 = w12
  c21 = w21
  c22 = w22
  c31 = w31
  c32 = w32
  c41 = w41
  c42 = w42

  ! Update the coefficients to account for derivatives
  select case(dir)
  case (1)
     c12 = 1.0_WP/(x1(i1+1)-x1(i1))
     c11 = -c12
  case (2)
     c22 = 1.0_WP/(x2(i2+1)-x2(i2))
     c21 = -c22
  case (3)
     c32 = 1.0_WP/(x3(i3+1)-x3(i3))
     c31 = -c32
  case (4)
     c42 = 1.0_WP/(x4(i4+1)-x4(i4))
     c41 = -c42
  end select

  ! Compute derivative in the direction 'dir' 
  ! and interpolation in the two others
  chemtable_4D_lookup_rho_der =  &
       c41*( c31*( c21*( c11*table(i1  ,i2  ,i3  ,i4  ,var)       &
                        +c12*table(i1+1,i2  ,i3  ,i4  ,var) )     &
                  +c22*( c11*table(i1  ,i2+1,i3  ,i4  ,var)       &
                        +c12*table(i1+1,i2+1,i3  ,i4  ,var) ) )   &
            +c32*( c21*( c11*table(i1  ,i2  ,i3+1,i4  ,var)       &
                        +c12*table(i1+1,i2  ,i3+1,i4  ,var) )     &
                  +c22*( c11*table(i1  ,i2+1,i3+1,i4  ,var)       &
                        +c12*table(i1+1,i2+1,i3+1,i4  ,var) ) ) ) &
      +c42*( c31*( c21*( c11*table(i1  ,i2  ,i3  ,i4+1,var)       &
                        +c12*table(i1+1,i2  ,i3  ,i4+1,var) )     &
                  +c22*( c11*table(i1  ,i2+1,i3  ,i4+1,var)       &
                        +c12*table(i1+1,i2+1,i3  ,i4+1,var) ) )   &
            +c32*( c21*( c11*table(i1  ,i2  ,i3+1,i4+1,var)       &
                        +c12*table(i1+1,i2  ,i3+1,i4+1,var) )     &
                  +c22*( c11*table(i1  ,i2+1,i3+1,i4+1,var)       &
                        +c12*table(i1+1,i2+1,i3+1,i4+1,var) ) ) )
  return
end function chemtable_4D_lookup_rho_der


! =============================================== !
! Find the maximum of a variable in the chemtable !
! =============================================== !
subroutine chemtable_4D_lookup_max(tag, R)
  use chemtable_4D
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), intent(out) :: R

  integer :: var

  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_lookup: unknown variable : '//trim(tag))
  end if

  ! Return the max
  R = maxval(table(:,:,:,:,var))

  return
end subroutine chemtable_4D_lookup_max


! =============================================== !
! Find the minimum of a variable in the chemtable !
! =============================================== !
subroutine chemtable_4D_lookup_min(tag, R)
  use chemtable_4D
  implicit none

  character(len=*), intent(in) :: tag
  real(WP), intent(out) :: R

  integer :: var

  ! Get the index of the variable
  do var=1, nvar
     if (trim(chem_name(var)).eq.trim(tag)) exit
  end do
  if (var.gt.nvar) then
     call die('chemtable_lookup: unknown variable : '//trim(tag))
  end if

  ! Return the max
  R = minval(table(:,:,:,:,var))

  return
end subroutine chemtable_4D_lookup_min
