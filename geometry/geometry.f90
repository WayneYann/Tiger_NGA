module geometry
  use string
  use precision
  use config
  implicit none

  ! imino ... imin .... imax ... imaxo
  ! <- nover --><-- nx --><-- nover ->
  ! <------------ nxo --------------->

  ! x(i) < xm(i) < x(i+1) < xm(i+1)
  ! <----- dx ----->
  !         <------- dxm ----->
  
  ! Overlap
  integer :: nover
  
  ! Cartesian  : x-direction
  ! Cylindrical: x-direction
  integer :: imino
  integer :: imin
  integer :: imaxo
  integer :: imax
  integer :: nx
  integer :: nxo

  ! Cartesian  : y-direction
  ! Cylindrical: r-direction
  integer :: jmino
  integer :: jmin
  integer :: jmaxo
  integer :: jmax
  integer :: ny
  integer :: nyo

  ! Cartesian  : z-direction
  ! Cylindrical: theta-direction
  integer :: kmino
  integer :: kmin
  integer :: kmaxo
  integer :: kmax
  integer :: nz
  integer :: nzo

  ! Node locations
  real(WP), dimension(:), pointer :: x
  real(WP), dimension(:), pointer :: y
  real(WP), dimension(:), pointer :: z

  ! Cell Center locations
  real(WP), dimension(:), pointer :: xm
  real(WP), dimension(:), pointer :: ym
  real(WP), dimension(:), pointer :: zm

  ! Length of cells
  real(WP), dimension(:), pointer :: dx
  real(WP), dimension(:), pointer :: dy
  real(WP), dimension(:), pointer :: dxm
  real(WP), dimension(:), pointer :: dym
  real(WP) :: dz
  
  ! Inverse of length of cells
  real(WP), dimension(:), pointer :: dxi
  real(WP), dimension(:), pointer :: dyi
  real(WP), dimension(:), pointer :: dxmi
  real(WP), dimension(:), pointer :: dymi
  real(WP) :: dzi
  
  ! For cylindrical cases : 1/r
  ! ymm(j) \approx y(j)
  real(WP), dimension(:), pointer :: yi
  real(WP), dimension(:), pointer :: ymi
  real(WP), dimension(:), pointer :: ymm
  real(WP), dimension(:), pointer :: ymmi
  real(WP), dimension(:), pointer :: dzi_u, dzi_v
  
  ! For surface integration
  real(WP), dimension(:), pointer :: dA
  real(WP), dimension(:), pointer :: dz_v
  
  ! Total length
  real(WP) :: xL,yL,zL
  
  ! Volume
  real(WP), dimension(:,:), pointer :: vol
  real(WP) :: vol_total
  
end module geometry


! ================================= !
! INITIALIZE the geometry           !
!                                   !
! -> Read config file               !
! -> Set up the masks               !
! -> Compute short hand notations   !
! -> Initialize the parallel module !
! -> Decompose the domain           !
! -> Initialize the blocks          !
! ================================= !
subroutine geometry_init
  use geometry
  use parser
  implicit none
  integer :: iunit

  ! Main geometry inits
  call parallel_init_filesystem
  call config_get_schemes(nover)
  call geometry_read_config(iunit)
  call masks_init(iunit)
  call geometry_notations
  call parallel_init_topology(xper,yper,zper)
  call partition_init
  call structure_init
  call borders_init
  
  ! Compute the metrics for the solvers
  call metric_generic_init
  call metric_velocity_conv_init
  call metric_velocity_visc_init
  
  ! Initialize the filter
  call filter_init
  
  ! Generate an equivalent unstructured grid
  call unstruct_init
  
  return
end subroutine geometry_init


! ===================== !
! Read the mesh & masks !
! ===================== !
subroutine geometry_read_config(iunit)
  use geometry
  use fileio
  use parser
  use parallel
  implicit none
  
  integer, intent(out) :: iunit
  integer :: ierr
  integer :: i,j,k
  character(len=str_medium) :: filename
  integer, dimension(MPI_STATUS_SIZE) :: status
  real(WP) :: delta_x,delta_y,delta_z

  ! Open the configuration file
  call parser_read('Configuration file',filename)
  filename = trim(mpiiofs) // ":" // trim(filename)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,mpi_info,iunit,ierr)
  if (ierr .ne. 0) stop "A configuration file is required"
  
  ! Read configuartion parameters
  call MPI_FILE_READ_ALL(iunit,simu_name,str_medium,MPI_CHARACTER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,icyl,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,xper,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,yper,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,zper,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,nx,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,ny,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(iunit,nz,1,MPI_INTEGER,status,ierr)
  
  ! Set the indices
  nxo = nx + 2*nover
  nyo = ny + 2*nover
  nzo = nz + 2*nover
  
  imino = 1
  imin  = imino + nover
  imax  = imin  + nx - 1
  imaxo = imax  + nover
  
  jmino = 1
  jmin  = jmino + nover
  jmax  = jmin  + ny - 1
  jmaxo = jmax  + nover
  
  kmino = 1
  kmin  = kmino + nover
  kmax  = kmin  + nz - 1
  kmaxo = kmax  + nover
  
  ! Allocate the arrays
  allocate(x(imino:imaxo+1))
  allocate(y(jmino:jmaxo+1))
  allocate(z(kmino:kmaxo+1))
  
  ! Read the data
  call MPI_FILE_READ_ALL(iunit,x(imin:imax+1),nx+1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(iunit,y(jmin:jmax+1),ny+1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(iunit,z(kmin:kmax+1),nz+1,MPI_REAL_WP,status,ierr)
  
  ! Compute total length
  xL = x(imax+1) - x(imin)
  yL = y(jmax+1) - y(jmin)
  zL = z(kmax+1) - z(kmin)
  
  ! Compute grid size in the overlap - X
  if (xper.ne.1) then
     delta_x = x(imin+1) - x(imin)
     do i=imin-1,imino,-1
        x(i) = x(i+1) - delta_x
     end do
     delta_x = x(imax+1) - x(imax)
     do i=imax+2,imaxo+1
        x(i) = x(i-1) + delta_x
     end do
  else
     do i=imin-1,imino,-1
        x(i) = x(i+nx) - xL
     end do
     do i=imax+2,imaxo+1
        x(i) = x(i-nx) + xL
     end do
  end if
  
  ! Compute grid size in the overlap - Y
  if (yper.ne.1) then
     if (icyl.eq.1) then
        do j=jmin-1,jmino,-1
           y(j) = 2.0_WP*y(jmin) - y(2*jmin-j)
        end do
        delta_y = y(jmax+1) - y(jmax)
        do j=jmax+2,jmaxo+1
           y(j) = y(j-1) + delta_y
        end do        
     else
        delta_y = y(jmin+1) - y(jmin)
        do j=jmin-1,jmino,-1
           y(j) = y(j+1) - delta_y
        end do
        delta_y = y(jmax+1) - y(jmax)
        do j=jmax+2,jmaxo+1
           y(j) = y(j-1) + delta_y
        end do
     end if
  else
     do j=jmin-1,jmino,-1
        y(j) = y(j+ny) - yL
     end do
     do j=jmax+2,jmaxo+1
        y(j) = y(j-ny) + yL
     end do
  end if
  
  ! Compute grid size in the overlap - Z
  if (zper.ne.1) then
     delta_z = z(kmin+1) - z(kmin)
     do k=kmin-1,kmino,-1
        z(k) = z(k+1) - delta_z
     end do
     delta_z = z(kmax+1) - z(kmax)
     do k=kmax+2,kmaxo+1
        z(k) = z(k-1) + delta_z
     end do
  else
     do k=kmin-1,kmino,-1
        z(k) = z(k+nz) - zL
     end do
     do k=kmax+2,kmaxo+1
        z(k) = z(k-nz) + zL
     end do
  end if
  
  ! Check consistency with cylindrical coordinates
  if (icyl.eq.1 .and. y(jmino).lt.0.0_WP .and. y(jmin).ne.0.0_WP) print*,'WARNING: y(jmin)=',y(jmin)
  
  return
end subroutine geometry_read_config


! ================================ !
! Compute the short hand notations !
! ================================ !
subroutine geometry_notations
  use geometry
  use math
  implicit none
  integer :: i,j,k
  real(WP) :: facz

  ! Allocate the arrays
  allocate(xm (imino:imaxo))
  allocate(ym (jmino:jmaxo))
  allocate(zm (kmino:kmaxo))
  allocate(ymm(jmino:jmaxo))

  allocate(dx(imino:imaxo))
  allocate(dy(jmino:jmaxo))

  allocate(dxi(imino:imaxo))
  allocate(dyi(jmino:jmaxo))

  allocate(dxm(imino:imaxo-1))
  allocate(dym(jmino:jmaxo-1))

  allocate(dxmi(imino:imaxo-1))
  allocate(dymi(jmino:jmaxo-1))
  
  allocate(yi(jmino:jmaxo))
  allocate(ymi(jmino:jmaxo))
  allocate(ymmi(jmino:jmaxo))
  
  allocate(vol(imino:imaxo,jmino:jmaxo)) ! Computed in metric_generic
  
  ! Compute location of cell centers
  do i=imino,imaxo
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  if (xper.ne.1) then
     xm(imaxo) = 2.0_WP*x(imaxo) - xm(imaxo-1)
  else
     xm(imaxo) = xm(imaxo-nx) + xL
  end if
  
  do j=jmino,jmaxo
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  if (yper.ne.1) then
     ym(jmaxo) = 2.0_WP*y(jmaxo) - ym(jmaxo-1)
  else
     ym(jmaxo) = ym(jmaxo-ny) + yL
  end if
  
  do k=kmino,kmaxo
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do
  if (zper.ne.1) then
     zm(kmaxo) = 2.0_WP*z(kmaxo) - zm(kmaxo-1)
  else
     zm(kmaxo) = zm(kmaxo-nz) + zL
  end if
  
  ! Mean staggered grid
  ! Very important
  do j=jmino+1,jmaxo-1
     ymm(j) = 0.5_WP*(y(j+1)+y(j-1))
  end do
  ymm(jmino) = y(jmino)
  ymm(jmaxo) = y(jmaxo)
  
  ! Compute short hand notations - x
  do i=imino,imaxo
     dx(i)   = x(i+1) - x(i)
     dxi(i)  = 1.0_WP/dx(i)
  end do
  do i=imino,imaxo-1
     dxm(i)  = xm(i+1) - xm(i)
     dxmi(i) = 1.0_WP/dxm(i)
  end do
  
  ! Compute short hand notations - y
  yi=0.0_WP
  ymi=0.0_WP
  ymmi=0.0_WP
  do j=jmino,jmaxo
     dy(j)   = y(j+1) - y(j)
     dyi(j)  = 1.0_WP/dy(j)
     if (y(j)  .ne.0.0_WP) yi(j)   = 1.0_WP/y(j)
     if (ym(j) .ne.0.0_WP) ymi(j)  = 1.0_WP/ym(j)
     if (ymm(j).ne.0.0_WP) ymmi(j) = 1.0_WP/ymm(j)
  end do
  do j=jmino,jmaxo-1
     dym(j)  = ym(j+1) - ym(j)
     dymi(j) = 1.0_WP/dym(j)
  end do

  ! THIS IS FOR MORINISHI'S AXIS TREATMENT
  ymmi(jmin) = ymi(jmin)

  ! Compute short hand notations - z
  dz = z(kmin+1) - z(kmin)
  dzi = 1.0_WP/dz
  
  ! Check if sector or full domain
  if (icyl.eq.1) then
     if (zL.lt.twoPi-dz) then
        isect = 1
     else
        isect = 0
     end if
  else
     isect = 0
  end if
  
  ! Additionnal notations
  if (icyl .eq. 1) then

     ! Allocate arrays
     allocate(dzi_u(jmino:jmaxo))
     allocate(dzi_v(jmino:jmaxo))
     allocate(dz_v(jmin:jmax+1))
     allocate(dA(jmin:jmax))

     facz  = (dz/2.0_WP)/sin(dz/2.0_WP)
     dzi_u = ymi*dzi
     dzi_v = dzi_u * facz

     ! Cell "surface"
     dz_v = y(jmin:jmax+1)*dz
     do j=jmin,jmax
        dA(j) = 0.5_WP*(y(j+1)**2-y(j)**2)*dz
     end do

  else

     ! Allocate the only array
     allocate(dzi_u(jmino:jmaxo))
     allocate(dz_v(jmin:jmax+1))
     allocate(dA(jmin:jmax))

     dzi_u = dzi
     dzi_v => dzi_u

     ! Cell "surface"
     dz_v = dz
     dA = dy(jmin:jmax)*dz

  end if
  
  return
end subroutine geometry_notations

