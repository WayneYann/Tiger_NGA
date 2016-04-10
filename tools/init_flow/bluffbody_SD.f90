module bluffbody_SD
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

end module bluffbody_SD

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine bluffbody_SD_grid
  use bluffbody_SD
  use parser
  use math
  implicit none

  integer :: i,j,k
  real(WP) :: Pi2,Ax,Bx,Cx,Ay,By,Cy,ytilde,s
  integer :: N_in,N_J,N_B
  real(WP) :: L_x,L_in,r_J,r_B,r_C 

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

  !---------------------------!
  ! Set up geometry parameters!
  ! --------------------------!
  
  ! Combustion length
  L_x = 0.72_WP ! 20X jet diameter
  ! Inlet length
  L_in = 7.2e-3_WP ! 2X jet diameter
  ! Jet radius
  r_J = 1.8e-3_WP
  ! Bluff-body radius
  r_B = 2.5e-2_WP
  ! Coflow radius
  r_C = 1.25e-1_WP ! 5X r_D

  !-------------------!
  ! Set up fixed grids!
  !-------------------!

  ! Number of points in inlet
  N_in = 64 !96
  ! Number of points in jet
  N_J = 16 !22
  ! Number of points in bluff-body
  N_B = 105   !160  
  
  print*,'nx,ny,nz =',nx,ny,nz
  print*,'N_in,N_J,N_B = ',N_in,N_J,N_B

!!$  !Choosing grid sizing based on N points in x
!!$  IF (nx.eq.384) THEN
!!$     Ax=.003199
!!$     Bx=1.014234
!!$     Cx=.0079
!!$  ELSE IF (nx.eq.512)THEN
!!$     Ax=.014667
!!$     Bx=1.00767
!!$  ELSE IF (nx.eq.768) THEN
!!$     Ax=.0255681
!!$     Bx=1.00440
!!$  END IF


  ! Create the grid
  ! Create x spacing
  ! In inlet
  do i=1,N_in
     x(i) = -real(N_in-i+1,WP) / real(N_in,WP) * L_in
     print*,'Inlet:',i,x(i)
  end do
  ! Stretch out
  Bx = 1.0142_WP
  Ax = L_x / (Bx**nx - BX**N_in)
  Cx = -Ax * Bx**N_in
  do i=N_in+1,nx+1
     x(i) = Ax * Bx**(i-1) + Cx
     x(N_in+1) = 0.0_WP
     print*,'x direction:',i,x(i)
  end do
  ! Check stretch rates
  print*,'Bx = ',Bx
  print*,'CHECK: X stretch rate Ok? (0.97-1.03)', (x(N_in+2)-x(N_in+1)) / (x(N_in+1)-x(N_in)) 


  !create y spacing
  ! In jet
  By = 0.98_WP
  Ay = r_J / (By**N_J - 1)
  Cy = -Ay
  do j=1,N_J
     y(j) = Ay * By**(j-1) + Cy
     print*,'In jet:',j,y(j)
  end do
  ! In bluff-body
  s = 1.21_WP
  do j=N_J+1,N_J+N_B
     ytilde = 2.0_WP * real(j-N_J-1,WP) / real(N_B,WP) - 1.0_WP
     y(j) = 0.5_WP * (r_B - r_J) * tanh(s * ytilde) / tanh(s) + r_J + 0.5_WP * (r_B - r_J)
     y(N_J+1) = r_J  
     print*,'In bluff-body:',j,y(j)
  end do
  ! Check stretch rates
  print*,'By = ',By
  print*,'s = ',s
  print*,'CHECK: Y Stretch rate in jet Ok? (0.97-1.03)', (y(N_J+2)-y(N_J+1)) / (y(N_J+1)-y(N_J)) 
  ! In coflow
  By = 1.05967_WP
  Ay = (r_C - r_B) / (By**ny - By**(N_J + N_B))
  Cy = r_C - Ay * By**ny
  do j=N_J+N_B+1,ny+1
     y(j) = Ay * By**(j-1) + Cy
     y(N_J+N_B+1) = r_B
     print*,'In coflow:',j,y(j)
  end do
  ! Check stretch rates
  print*,'By = ',By
  print*,'CHECK: Y stretch rate in coflow Ok? (0.97-1.03)', (y(N_J+N_B+2)-y(N_J+N_B+1)) / (y(N_J+N_B+1)-y(N_J+N_B)) 

!!$     if (j.le.N_J) then
!!$        Ay=-.0065169_WP
!!$        By=.98_WP
!!$        Cy = 
!!$        y(j) =  Ay*By**(j-1)+Cy
!!$     else if (j.le.N_) then
!!$        ytilde = 2.0_WP*real(j-17,WP)/real(105,WP)-1.0_WP
!!$        s=1.21_WP
!!$        init_r= .0018_WP
!!$        height = .0232_WP
!!$        y(j) = .5_WP*height*tanh(s*ytilde)/tanh(s)+ init_r + 0.5_WP*height
!!$     else if (j.ge.123) then
!!$        Ay=1.523153766e-6_WP !1.4065E-7_WP
!!$        By=1.059673_WP !1.0777483_WP
!!$        Cy=-2.3307444e-2_WP !-2.3696E-2_WP
!!$        y(j) =  Ay*By**(j-1)-Cy
!!$     end if

  !create z spacing
  Pi2 = 2.0_WP * acos(-1.0_WP)
  do k=1,nz+1
     z(k) = real(k-1,WP) * Pi2 / real(nz,WP)
  end do

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
        if ((j.gt.N_J).and.(i.le.N_in).and.(j.le.N_J+N_B) ) then
           mask(i,j) = 1
        end if
     end do
  end do
  mask(nx,:) = mask(nx-1,:)
 
  return
end subroutine bluffbody_SD_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine bluffbody_SD_data
  use bluffbody_SD
  implicit none

  ! Allocate the array data
  nvar = 18
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  data(:,:,:, 5) = 0.0_WP;      names( 5) = 'ZMIX'
  data(:,:,:, 6) = 1.179640_WP; names( 6) = 'RHO'
  data(:,:,:, 7) = 0.0_WP;      names( 7) = 'dRHO'
  data(:,:,:, 8) = 0.0_WP;      names( 8) = 'PROG'
  data(:,:,:, 9) = 0.0_WP;      names( 9) = 'ENTH'
  data(:,:,:,10) = 0.0_WP;      names(10) = 'ZMIX2'
  data(:,:,:,11) = 1.0e-20_WP;      names(11) = 'S_M00'
  data(:,:,:,12) = 4.0e-19_WP;      names(12) = 'S_M10'
  data(:,:,:,13) = 1.1696e-19_WP;      names(13) = 'S_M01'
  data(:,:,:,14) = 1.0e-20_WP;      names(14) = 'S_N00'
  data(:,:,:,15) = 0.0_WP;      names(15) = 'PAH'
  data(:,:,:,16) = 1.0e-40_WP;      names(16) = 'S_MSQ'
  data(:,:,:,17) = 0.0_WP;      names(17) = 'VISC'
  data(:,:,:,18) = 0.0_WP;      names(18) = 'DIFF'  


  
  ! Create them
  U = 23.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  
  return
end subroutine bluffbody_SD_data

