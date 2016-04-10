module bluffbodyrlb
  use precision
  use param
  implicit none
  
  ! Length and diameter of the domain
  real(WP) :: L,diam
  ! Length and diameter of the pipe
  real(WP) :: pipe_L,pipe_diam
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P

end module bluffbodyrlb

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine bluffbodyrlb_grid
  use bluffbodyrlb
  use parser
  use math
  implicit none

  integer :: i,j,k
  real(WP) :: Pi2,Ax,Bx,Cx,Ay,By,Cy,xref,ytilde,s,init_r,height

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

!!$  call parser_read('Full length',L)
!!$  call parser_read('Full diameter',diam)
!!$
!!$  call parser_read('Pipe length',pipe_L)
!!$  call parser_read('Pipe diameter',pipe_diam)

  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 1

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

!Choosing grid sizing based on N points in x
IF (nx.eq.384) THEN
Ax=.003199
Bx=1.014234
Cx=.0079
ELSE IF (nx.eq.512)THEN
Ax=.014667
Bx=1.00767
ELSE IF (nx.eq.768) THEN
Ax=.0255681
Bx=1.00440
END IF

xref=0.0001125

  ! Create the grid
!create x spacing
x(1)=-.0072_WP
  do i=2,nx+1
IF (i.le.65) THEN
x(i) = x(i-1) + xref
ELSE     
x(i) = Ax*Bx**(i-1)-Cx
END IF
end do
do i=1,nx+1
print*, i,x(i)
end do

!create y spacing
  do j=1,ny+1
 !    y(j) = (j-1)*diam/(2.0_WP*ny)
IF (j.le.17) THEN
 Ay=-.0065169_WP
 By=.98_WP
 y(j) =  Ay*By**(j-1)-Ay
ELSE IF (j.le.122) THEN
 ytilde = 2.0_WP*real(j-17,WP)/real(105,WP)-1.0_WP
 s=1.21_WP
 init_r= .0018_WP
 height = .0232_WP
 y(j) = .5_WP*height*tanh(s*ytilde)/tanh(s)+ init_r + 0.5_WP*height
ELSE IF (j.ge.123) THEN
 Ay=1.523153766e-6_WP !1.4065E-7_WP
 By=1.059673_WP !1.0777483_WP
 Cy=-2.3307444e-2_WP !-2.3696E-2_WP
 y(j) =  Ay*By**(j-1)-Cy
END IF
  end do
do j=1,ny+1
print*, j,y(j)
end do

!create z spacing
  Pi2 = 2.0_WP*acos(-1.0_WP)
  do k=1,nz+1
     z(k) = real(k-1,WP)*Pi2/real(nz,WP)
  end do
!!$  do k=1,nz
!!$     zm(k)= 0.5_WP*(z(k)+z(k+1))
!!$  end do
!!$
!!$  ! Create the mid points
!!$  do i=1,nx
!!$     xm(i)= 0.5_WP*(x(i)+x(i+1))
!!$  end do
!!$  do j=1,ny
!!$     ym(j)= 0.5_WP*(y(j)+y(j+1))
!!$  end do
!!$  do k=1,nz
!!$     zm(k)= 0.5_WP*(z(k)+z(k+1))
!!$  end do

  ! Create the masks
  mask = 0
  do i=1,nx-1
     do j=1,ny
        if ((y(j).ge.(0.0018_WP-1.0e-6_WP)).and.(x(i+1).le.0.0_WP).and.(y(j+1).le.0.025_WP) ) then
           mask(i,j) = 1
        end if
     end do
  end do
  mask(nx,:) = mask(nx-1,:)
  
  return
end subroutine bluffbodyrlb_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine bluffbodyrlb_data
  use bluffbodyrlb
  implicit none

  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Create them
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  
  return
end subroutine bluffbodyrlb_data

