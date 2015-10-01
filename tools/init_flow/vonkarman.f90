module vonkarman
  use precision
  use param
  use fileio
  implicit none

  ! Length and Width of the domain
  real(WP) :: L,width,H
  ! Radius and position of the cylinder
  real(WP) :: x_pos,y_pos,R,M
  
  ! Grid Parameters
  real(WP) ::yf,ny1,nx1,nx2,ny2,nf

  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P

  ! Use Immersed Boundary Option
  integer :: i_use_ibm, io
 
end module vonkarman


! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine vonkarman_grid
  use vonkarman
  use parser
  implicit none
  
  integer :: i,j,k
  integer, dimension(3) :: id_cpus,n_ms,n_ps,n_cpus,nps_global,npe_global 
  real(WP) :: dy,dx,ay,ax1,ax2

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  
  call parser_read('Length x',L)
  call parser_read('Height y',H)
  call parser_read('Width z',width)
  
  call parser_read('Cylinder Radius',R)
  call parser_read('Cyl x-pos',x_pos)
  call parser_read('Cyl y-pos',y_pos)
  
  ! Read in Fine Mesh parameters
  call parser_read('fine grid points',nf)
  call parser_read('points before cylinder',nx1) 

  ! Invoke Immersed Boundary module
  call parser_read('Use IB mod',i_use_ibm)

  ! Set Periodicity
  xper = 0
  yper = 0
  zper = 1
  
  ! Cylindrical
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Size and spacing of the fine grid
  yf = 2.0_WP * R
  dy = 2.0_WP * yf / nf  

  ! Create the grid
  ! Define # of y points having coarse grid in each half plane  
  nx2 = nx - nf - nx1
  ny1 = nf/2
  ny2 = ny/2 - ny1
  
 
  ! Stretching in y
  ay = 1.1_WP
  do i=1,40
     ay = (((H/2 - yf)*(ay - 1.0_WP)/dy) + 1)**(1.0_WP/ny2)
  end do
  print*,'ay = ',ay
 
  ! ================= !
  ! Grid generation y !
  ! ================= !
  
  ! Set centerline to 0
  y(ny/2+1) = 0.0_WP
  dy = 2.0_WP * yf /nf
  ! Fine grid (top)
  do j=ny/2+2,ny/2+2+ny1
     y(j) = y(j-1) + dy
  end do
  ! Coarse grid (top)
  do j = ny/2 +ny1 +2,ny+1
     y(j) = y(j-1) + dy
     dy = ay*dy
  end do
  
  ! Use Symmetry in y for the bottom part
  do j= 1,ny/2
     y(j) = -1.0_WP * y(ny+2-j)
  end do  
  
  ! ================= !
  ! Grid generation x !
  ! ================= !

  ! Fine grid (center)
  dx = 2.0_WP * yf / nf  
  do i=nx1+1,nx1+nf+1
     x(i) = x_pos - (3.0_WP/4.0_WP)*yf + dx*(i-(nx1+1))*1.0_WP
  end do
  
  ! Stretching in x (before the cyl)
  dx  = 2.0_WP * yf / nf
  ax1 = 1.1_WP
  do i=1,40
     ax1 = (((x_pos - (3.0_WP/4.0_WP)*yf)*(ax1 - 1.0_WP)/dx) + 1)**(1.0_WP/nx1)
  end do
  print*,'ax1 = ',ax1
  do i=nx1,1,-1
     x(i) = x(i+1) - dx
     dx   = ax1 * dx
  end do
  
  ! Stretching in x (after the cyl)
  dx  = 2.0_WP * yf / nf
  ax2 = 1.1_WP
  do i=1,40
     ax2=(((L - x_pos - 1.25_WP*yf)*(ax2 - 1.0_WP)/dx) + 1)**(1.0_WP/nx2)
  end do
  print*,'ax2 = ',ax2
  do i=nx1+nf+2,nx+1
     x(i) = x(i-1) + dx
     dx   = ax2 * dx
  end do
  x = x - x_pos
  
  ! ================= !
  ! Grid generation z !
  ! ================= !
  do k=1,nz+1
     z(k) = (k-1)*width/nz - 0.5_WP*width
  end do
  
  ! Create the mid points 
  do i=1,nx
     xm(i)= 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j)= 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k)= 0.5_WP*(z(k)+z(k+1))
  end do


  !--------------------------------------------------------------------
  !  Routine to initialize the IB variables
  !
  !  Parameters
  !  - Input
  !    imode:
  !        1: single cpu preprocessor mode (id_cpus will be ignored)
  !        2: runtime mode with file output
  !        3: runtime mode without file output
  !    nps,npe: the start/end index of grid points in each direction
  !    xs,ys,zs: positions of grid points
  !    ipr_s,jpr_s,kpr_s: option for periodic b.c.
  !                       => 1: periodic b.c.
  !    i_staggered_s: option for staggered mesh configuration
  !                => 0: collocated 1: kang's code 2: chuck's code
  !    i_type_index_s: option for index ordering
  !                => 1: kang's code 2: chuck's code
  !    id_cpus: rank of this cpu in each direction
  !    n_cpus: # of cpus in each direction
  !    npad_ms: # of buffer grid cells on the lower side of domain to be written
  !    npad_ps: # of buffer grid cells on the upper side of domain to be written
  !    nps_global_s,npe_global_s: the start/end index of grid points
  !                               in the global domain in each direction
  !    nf: file id for input/output
  !    add_bnd,search_bnd: small margin to a real value in order to make
  !                        two very similar values recognized as the same
  !    i_tri_mod: option to lower dimension of space where interpolation of
  !               the velocity is made, in order to avoid some case that an
  !               interpolated velocity component is determined by other
  !               interpolated velocity components
  !               (normally set to 0)
  !    s_ratio: a value > 1. used when i_tri_mod=1
  !--------------------------------------------------------------------
  if (i_use_ibm .eq. 1) then
     mask = 0
     id_cpus(1)=1
     id_cpus(2)=1
     id_cpus(3)=1
     n_cpus(1)=1
     n_cpus(2)=1
     n_cpus(3)=1
     ! running single cpu so no padding cells
     n_ms(1)=0
     n_ms(2)=0
     n_ms(3)=0
     n_ps(1)=0
     n_ps(2)=0
     n_ps(3)=0
     nps_global(1)=1
     nps_global(2)=1
     nps_global(3)=1
     npe_global(1)=nx+1
     npe_global(2)=ny+1
     npe_global(3)=nz+1
     
     io = iopen()
     !call ibm_init_preprocess(1,1,nx+1,1,ny+1,1,nz+1,x,y,z,xper,yper,zper, &
     !     1,1,id_cpus,n_cpus,n_ms,n_ps,nps_global,npe_global,io,           &
     !     1.0e-12_WP,1.0e-12_WP,0,1.1_WP)

  else
     ! Use stair stepping
     ! Create the masks
     mask = 0
     do i=1,nx-1
        do j=1,ny
           if (sqrt((xm(i))**2+(ym(j))**2).le.R) then
              mask(i,j) = 1
           end if
        end do
     end do
  end if
  
  return
end subroutine vonkarman_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine vonkarman_data
  use vonkarman
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
end subroutine vonkarman_data
