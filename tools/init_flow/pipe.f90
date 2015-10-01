module pipe
  use precision
  use param
  implicit none
  
  ! Number of pipes
  integer :: npipes
  
  ! Number of points on y
  integer, dimension(:), pointer :: ny_pipe
  integer, dimension(:), pointer :: j_start
  integer, dimension(:), pointer :: j_end
  
  ! Length of the domain
  real(WP) :: L,amp
  real(WP), dimension(:,:), pointer :: radius
  real(WP), dimension(:),   pointer :: height
  
  ! Stretching
  real(WP), dimension(:), pointer :: r
  
  ! Mean U/W velocity
  real(WP), dimension(:), pointer :: Umean,Wmean
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  
end module pipe

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine pipe_grid
  use pipe
  use parser
  use math
  implicit none
  
  integer :: i,j,k,pp
  real(WP) :: ytilde
  
  ! Mesh points
  call parser_read('nx',nx)
  call parser_read('nz',nz)
  
  ! Length
  call parser_read('Length',L)
  
  ! Number of pipes
  call parser_getsize('ny',npipes)
  allocate(ny_pipe(npipes))
  allocate(j_start(npipes))
  allocate(j_end(npipes))
  allocate(radius(npipes,2))
  allocate(height(npipes))
  allocate(r(npipes))
  allocate(Umean(npipes))
  allocate(Wmean(npipes))
  
  ! Read each pipe
  call parser_read('ny',ny_pipe)
  call parser_read('Radius',radius)
  height(:) = radius(:,2)-radius(:,1)
  call parser_read('Stretching',r)
  
  do pp=1,npipes-1
     if (radius(pp,2).ge.radius(pp+1,1)) print*,"Error in pipe definition"
  end do
  
  ! Set the periodicity
  xper = 1
  yper = 0
  zper = 1
  
  ! Cylindrical
  icyl = 1
  
  ! Create the indices in y
  ! First pipe
  if (radius(1,1).eq.0.0_WP) then
     ! First pipe is a real pipe (not annular)
     j_start(1) = 1
  else
     ! First pipe is annular, centerline is a wall
     j_start(1) = 2
  end if
  do pp=1,npipes-1
     j_end(pp) = j_start(pp)+ny_pipe(pp)
     j_start(pp+1) = j_end(pp)+1
  end do
  ! Last pipe
  j_end(npipes) = j_start(npipes)+ny_pipe(npipes)
  ! Total number of points
  ny = j_end(npipes)
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid in x and z
  do i=1,nx+1
     x(i) = (i-1.0_WP)*L/real(nx)
  end do
  if (nz.eq.1) then
     do k=1,nz+1
        z(k) = (k-1)*twoPi/128.0_WP
     end do
  else  
     do k=1,nz+1
        z(k) = (k-1)*twoPi/real(nz)
     end do
  end if

  ! Create the grid in y
  do pp=1,npipes
     if (radius(pp,1).eq.0.0_WP) then
        do j=j_start(pp),j_end(pp)
           ytilde = (j-j_start(pp))/real(ny_pipe(pp),WP)
           y(j) = height(pp)*tanh(r(pp)*ytilde)/tanh(r(pp))+radius(pp,1)
        end do
     else
        do j=j_start(pp),j_end(pp)
           ytilde = 2.0_WP*(j-j_start(pp))/real(ny_pipe(pp),WP)-1.0_WP
           y(j) = 0.5_WP*height(pp)*tanh(r(pp)*ytilde)/tanh(r(pp))+radius(pp,1)+0.5_WP*height(pp)
        end do
     end if
  end do
  ! Update points if necessary
  if (j_start(1).ne.1) y(1)=2.0_WP*y(2)-y(3)
  y(ny+1) = 2.0_WP*y(ny)-y(ny-1)
  
  ! Create the mid points
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do
  
  ! Create the masks
  mask = 0
  if (radius(1,1).ne.0.0_WP) mask(:,1)=1
  do pp=1,npipes
     mask(:,j_end(pp)) = 1
  end do
  
  return
end subroutine pipe_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine pipe_data
  use pipe
  use parser
  use math
  implicit none
  
  integer  :: i,j,k,pp
  real(WP) :: rnd,y1,y2,coeff,Unorm
  integer  :: laminar
  
  ! Initialize the random number generator
  call random_init
  
  ! Read the means
  call parser_read('Mean U Velocity',Umean)
  call parser_read('Mean W Velocity',Wmean)
  
  ! Read in the amplitude of random numbers
  call parser_read('Laminar initial flow',laminar)
  call parser_read('Fluctuation rel. amp.',amp)
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U    => data(:,:,:,1); names(1) = 'U'
  V    => data(:,:,:,2); names(2) = 'V'
  W    => data(:,:,:,3); names(3) = 'W'
  P    => data(:,:,:,4); names(4) = 'P'
    
  ! Create them
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP

  do pp=1,npipes
     y2 = y(j_end(pp))
     
     ! Laminar profile
     if (radius(pp,1).eq.0.0_WP) then
        y1 = -y(j_end(pp))
        coeff = 8.0_WP
     else
        y1 = y(j_start(pp))
        coeff = 6.0_WP
     end if
     do j=j_start(pp),j_end(pp)
        U(:,j,:) = Umean(pp) * coeff * (ym(j)-y1)*(y2-ym(j))/(y2-y1)**2
        W(:,j,:) = Wmean(pp) * coeff * (ym(j)-y1)*(y2-ym(j))/(y2-y1)**2
     end do

     ! For faster transition
     if (laminar.ne.1) then
        do j=j_start(pp),j_end(pp)
           Unorm = sqrt(Umean(pp)**2+Wmean(pp)**2)
           ! Fluctuations in X for W
           do i=1,nx
              W(i,j,:) = W(i,j,:) + amp * Unorm * cos(8.0_WP*twoPi*xm(i)/L)
           end do
           ! Fluctuations in Z for U
           do k=1,nz
              U(:,j,k) = U(:,j,k) + amp * Unorm * cos(8.0_WP*zm(k))
           end do
        end do
     end if
     
     ! Random values
     do i=1,nx
        do j=j_start(pp),j_end(pp)
           do k=1,nz
              call random_number(rnd)
              U(i,j,k) = U(i,j,k) + (rnd-0.5_WP)*amp*U(i,j,k)
              call random_number(rnd)
              V(i,j,k) = V(i,j,k) + (rnd-0.5_WP)*amp*V(i,j,k)
              call random_number(rnd)
              W(i,j,k) = W(i,j,k) + (rnd-0.5_WP)*amp*W(i,j,k)
           end do
        end do
     end do

     U(:,j_end(pp),:) = 0.0_WP
     W(:,j_end(pp),:) = 0.0_WP
  end do

  return
end subroutine pipe_data

