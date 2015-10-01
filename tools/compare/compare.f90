module compare
  use precision
  use string
  use fileio
  implicit none
  
  ! Sizes
  integer :: nx1,ny1,nz1
  integer :: nx2,ny2,nz2
  integer :: xper1,yper1,zper1
  integer :: xper2,yper2,zper2
  integer :: icyl1
  integer :: icyl2
  integer :: nvar
  real(WP) :: dt,time
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: config1
  character(len=str_medium) :: config2
  integer, dimension(:,:), pointer :: mask1
  integer, dimension(:,:), pointer :: mask2
  
  ! Data
  real(WP), dimension(:,:,:), pointer :: data1
  real(WP), dimension(:,:,:), pointer :: data2
  real(WP), dimension(:,:,:), pointer :: data1i
  
  ! Mesh
  real(WP), dimension(:), pointer :: x1,xm1
  real(WP), dimension(:), pointer :: y1,ym1
  real(WP), dimension(:), pointer :: z1,zm1
  real(WP), dimension(:), pointer :: x2,xm2
  real(WP), dimension(:), pointer :: y2,ym2
  real(WP), dimension(:), pointer :: z2,zm2
  
  ! Volume
  real(WP), dimension(:,:), pointer :: vol
  real(WP) :: vol_total
  
end module compare


program comparator
  use compare
  implicit none
  
  character(len=str_medium) :: fconfig1
  character(len=str_medium) :: fconfig2
  character(len=str_medium) :: fdata1
  character(len=str_medium) :: fdata2
  integer :: iunit,iunit1,iunit2,var
  integer :: i,j,k,ierr
  real(WP) :: err_1_U,err_1_V,err_1_W,err_1_G
  real(WP) :: err_2_U,err_2_V,err_2_W,err_2_G
  real(WP) :: err_1v_U,err_1v_V,err_1v_W,err_1v_G
  real(WP) :: err_2v_U,err_2v_V,err_2v_W,err_2v_G
  real(WP) :: err_i_U,err_i_V,err_i_W,err_i_G
  
  ! Get file names from user
  print*,'==============================='
  print*,'| ARTS - Data file comparator |'
  print*,'==============================='
  print*
  print "(a27,$)", " Enter first config file : "
  read "(a)", fconfig1
  print "(a28,$)", " Enter second config file : "
  read "(a)", fconfig2
  print "(a25,$)", " Enter first data file : "
  read "(a)", fdata1
  print "(a26,$)", " Enter second data file : "
  read "(a)", fdata2
  
  ! Read the source config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig1),"r",ierr)
  call BINARY_FILE_READ(iunit,config1,str_medium,kind(config1),ierr)
  call BINARY_FILE_READ(iunit,icyl1,1,kind(icyl1),ierr)
  call BINARY_FILE_READ(iunit,xper1,1,kind(xper1),ierr)
  call BINARY_FILE_READ(iunit,yper1,1,kind(yper1),ierr)
  call BINARY_FILE_READ(iunit,zper1,1,kind(zper1),ierr)
  call BINARY_FILE_READ(iunit,nx1,1,kind(nx1),ierr)
  call BINARY_FILE_READ(iunit,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit,nz1,1,kind(nz1),ierr)
  print*,'Source file'
  print*,'Config : ',config1
  print*,'Grid :',nx1,'x',ny1,'x',nz1
  allocate(x1(nx1+1),y1(ny1+1),z1(nz1+1))
  allocate(mask1(nx1,ny1))
  call BINARY_FILE_READ(iunit,x1,nx1+1,kind(x1),ierr)
  call BINARY_FILE_READ(iunit,y1,ny1+1,kind(y1),ierr)
  call BINARY_FILE_READ(iunit,z1,nz1+1,kind(z1),ierr)
  call BINARY_FILE_READ(iunit,mask1,nx1*ny1,kind(mask1),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Read the destination config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig2),"r",ierr)
  call BINARY_FILE_READ(iunit,config2,str_medium,kind(config2),ierr)
  call BINARY_FILE_READ(iunit,icyl2,1,kind(icyl2),ierr)
  call BINARY_FILE_READ(iunit,xper2,1,kind(xper2),ierr)
  call BINARY_FILE_READ(iunit,yper2,1,kind(yper2),ierr)
  call BINARY_FILE_READ(iunit,zper2,1,kind(zper2),ierr)
  call BINARY_FILE_READ(iunit,nx2,1,kind(nx2),ierr)
  call BINARY_FILE_READ(iunit,ny2,1,kind(ny2),ierr)
  call BINARY_FILE_READ(iunit,nz2,1,kind(nz2),ierr)
  print*,'Destination file'
  print*,'Config : ',config2
  print*,'Grid :',nx2,'x',ny2,'x',nz2
  allocate(x2(nx2+1),y2(ny2+1),z2(nz2+1))
  allocate(mask2(nx2,ny2))
  call BINARY_FILE_READ(iunit,x2,nx2+1,kind(x2),ierr)
  call BINARY_FILE_READ(iunit,y2,ny2+1,kind(y2),ierr)
  call BINARY_FILE_READ(iunit,z2,nz2+1,kind(z2),ierr)
  call BINARY_FILE_READ(iunit,mask2,nx2*ny2,kind(mask2),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Finish to create the meshes
  allocate(xm1(nx1),ym1(ny1),zm1(nz1))
  do i=1,nx1
     xm1(i) = 0.5_WP*(x1(i)+x1(i+1))
  end do
  do j=1,ny1
     ym1(j) = 0.5_WP*(y1(j)+y1(j+1))
  end do
  do k=1,nz1
     zm1(k) = 0.5_WP*(z1(k)+z1(k+1))
  end do
  allocate(xm2(nx2),ym2(ny2),zm2(nz2))
  do i=1,nx2
     xm2(i) = 0.5_WP*(x2(i)+x2(i+1))
  end do
  do j=1,ny2
     ym2(j) = 0.5_WP*(y2(j)+y2(j+1))
  end do
  do k=1,nz2
     zm2(k) = 0.5_WP*(z2(k)+z2(k+1))
  end do
  
  ! Compute volume
  allocate(vol(nx2,ny2))
  if (icyl2.eq.1) then
     do j=1,ny2
        do i=1,nx2
           if (mask2(i,j).eq.1) then
              vol(i,j) = 0.0_WP
           else
              vol(i,j) = 0.5_WP*(x2(i+1)-x2(i))*(y2(j+1)**2-y2(j)**2)*(z2(2)-z2(1))
           end if
        end do
     end do
  else
     do j=1,ny2
        do i=1,nx2
           if (mask2(i,j).eq.1) then
              vol(i,j) = 0.0_WP
           else
              vol(i,j) = (x2(i+1)-x2(i))*(y2(j+1)-y2(j))*(z2(2)-z2(1))
           end if
        end do
     end do
  end if
  vol_total = sum(vol)
  
  ! Read the source data file
  call BINARY_FILE_OPEN(iunit1,trim(fdata1),"r",ierr)
  call BINARY_FILE_READ(iunit1,nx1, 1,kind(nx1), ierr)
  call BINARY_FILE_READ(iunit1,ny1, 1,kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1,nz1, 1,kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_READ(iunit1,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables in source data file: ',names
  
  ! Read the new data file
  call BINARY_FILE_OPEN(iunit2,trim(fdata2),"r",ierr)
  call BINARY_FILE_READ(iunit2,nx2, 1,kind(nx2), ierr)
  call BINARY_FILE_READ(iunit2,ny2, 1,kind(ny2), ierr)
  call BINARY_FILE_READ(iunit2,nz2, 1,kind(nz2), ierr)
  call BINARY_FILE_READ(iunit2,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_READ(iunit2,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_READ(iunit2,time,1,kind(time),ierr)
  do var=1,nvar
     call BINARY_FILE_READ(iunit2,names(var),str_short,kind(names),ierr)
  end do
  
  ! Prepare for comparison
  allocate(data1(nx1,ny1,nz1))
  allocate(data2(nx2,ny2,nz2))
  allocate(data1i(nx2,ny2,nz2))
  err_1_U = 0.0_WP; err_2_U = 0.0_WP; err_1v_U = 0.0_WP; err_2v_U = 0.0_WP; err_i_U = 0.0_WP
  err_1_V = 0.0_WP; err_2_V = 0.0_WP; err_1v_V = 0.0_WP; err_2v_W = 0.0_WP; err_i_V = 0.0_WP
  err_1_W = 0.0_WP; err_2_W = 0.0_WP; err_1v_W = 0.0_WP; err_2v_V = 0.0_WP; err_i_W = 0.0_WP
  
  ! Interpolate
  do var=1,nvar
     ! Read variables
     call BINARY_FILE_READ(iunit1,data1,nx1*ny1*nz1,kind(data1),ierr)
     call BINARY_FILE_READ(iunit2,data2,nx2*ny2*nz2,kind(data2),ierr)
     
     ! Perform interpolation & compute error
     select case (trim(adjustl(names(var))))
     case ('U')
        call interp_data(data1,data1i,'U')
        do k=1,nz2
           do j=1,ny2
              do i=1,nx2
                 err_1_U = err_1_U + abs(data1i(i,j,k)-data2(i,j,k))
                 err_2_U = err_2_U + (data1i(i,j,k)-data2(i,j,k))**2
                 err_1v_U = err_1v_U + abs(data1i(i,j,k)-data2(i,j,k))*vol(i,j)
                 err_2v_U = err_2v_U + (data1i(i,j,k)-data2(i,j,k))**2*vol(i,j)
                 err_i_U = max(err_i_U,abs(data1i(i,j,k)-data2(i,j,k)))
              end do
           end do
        end do
     case ('V')
        call interp_data(data1,data1i,'V')
        do k=1,nz2
           do j=1,ny2
              do i=1,nx2
                 err_1_V = err_1_V + abs(data1i(i,j,k)-data2(i,j,k))
                 err_2_V = err_2_V + (data1i(i,j,k)-data2(i,j,k))**2
                 err_1v_V = err_1v_V + abs(data1i(i,j,k)-data2(i,j,k))*vol(i,j)
                 err_2v_V = err_2v_V + (data1i(i,j,k)-data2(i,j,k))**2*vol(i,j)
                 err_i_V = max(err_i_V,abs(data1i(i,j,k)-data2(i,j,k)))
              end do
           end do
        end do
     case ('W')
        call interp_data(data1,data1i,'W')
        do k=1,nz2
           do j=1,ny2
              do i=1,nx2
                 err_1_W = err_1_W + abs(data1i(i,j,k)-data2(i,j,k))
                 err_2_W = err_2_W + (data1i(i,j,k)-data2(i,j,k))**2
                 err_1v_W = err_1v_W + abs(data1i(i,j,k)-data2(i,j,k))*vol(i,j)
                 err_2v_W = err_2v_W + (data1i(i,j,k)-data2(i,j,k))**2*vol(i,j)
                 err_i_W = max(err_i_W,abs(data1i(i,j,k)-data2(i,j,k)))
              end do
           end do
        end do
     case ('G')
        call interp_data(data1,data1i,'SC')
        do k=1,nz2
           do j=1,ny2
              do i=1,nx2
                 err_1_G = err_1_G + abs(data1i(i,j,k)-data2(i,j,k))
                 err_2_G = err_2_G + (data1i(i,j,k)-data2(i,j,k))**2
                 err_1v_G = err_1v_G + abs(data1i(i,j,k)-data2(i,j,k))*vol(i,j)
                 err_2v_G = err_2v_G + (data1i(i,j,k)-data2(i,j,k))**2*vol(i,j)
                 err_i_G = max(err_i_G,abs(data1i(i,j,k)-data2(i,j,k)))
              end do
           end do
        end do
     end select
     
  end do
  
  ! Rescale
  err_1v_U = err_1v_U/vol_total
  err_1v_V = err_1v_V/vol_total
  err_1v_W = err_1v_W/vol_total
  err_1v_G = err_1v_G/vol_total
  err_2v_U = sqrt(err_2v_U/vol_total)
  err_2v_V = sqrt(err_2v_V/vol_total)
  err_2v_W = sqrt(err_2v_W/vol_total)
  err_2v_G = sqrt(err_2v_G/vol_total)
  err_1_U = err_1_U/real(nx2*ny2*nz2,WP)
  err_1_V = err_1_V/real(nx2*ny2*nz2,WP)
  err_1_W = err_1_W/real(nx2*ny2*nz2,WP)
  err_1_G = err_1_G/real(nx2*ny2*nz2,WP)
  err_2_U = sqrt(err_2_U/real(nx2*ny2*nz2,WP))
  err_2_V = sqrt(err_2_V/real(nx2*ny2*nz2,WP))
  err_2_W = sqrt(err_2_W/real(nx2*ny2*nz2,WP))
  err_2_G = sqrt(err_2_G/real(nx2*ny2*nz2,WP))
  
  ! Print results
  print*,"L1(U)",err_1_U
  print*,"L1(V)",err_1_V
  print*,"L1(W)",err_1_W
  print*,"L1(G)",err_1_G

  print*,"L2(U)",err_2_U
  print*,"L2(V)",err_2_V
  print*,"L2(W)",err_2_W
  print*,"L2(G)",err_2_G

  print*,"L1v(U)",err_1v_U
  print*,"L1v(V)",err_1v_V
  print*,"L1v(W)",err_1v_W
  print*,"L1v(G)",err_1v_G

  print*,"L2v(U)",err_2v_U
  print*,"L2v(V)",err_2v_V
  print*,"L2v(W)",err_2v_W
  print*,"L2v(G)",err_2v_G

  print*,"Linf(U)",err_i_U
  print*,"Linf(V)",err_i_V
  print*,"Linf(W)",err_i_W
  print*,"Linf(G)",err_i_G
  
  ! Close files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
end program comparator


! --------------------------------------- !
! Interpolate the data between two arrays !
! For a given staggering 'dir'            !
! --------------------------------------- !
subroutine interp_data(A1,A2,dir)
  use compare
  implicit none
  
  real(WP), dimension(nx1,ny1,nz1), intent(in)  :: A1
  real(WP), dimension(nx2,ny2,nz2), intent(out) :: A2
  character(len=*) :: dir

  real(WP), dimension(:), pointer :: x,y,z
  integer  :: nx,ny,nz
  integer  :: i2,j2,k2
  integer  :: i1,j1,k1
  integer  :: ni,nj,nk
  real(WP) :: xd,yd,zd
  
  real(WP), dimension(0:1,0:1,0:1) :: ww
  real(WP) :: wx1,wx2,wy1,wy2,wz1,wz2

  ! dir = U  => x  - ym - zm
  ! dir = V  => xm - y  - zm
  ! dir = W  => xm - ym - z
  ! dir = SC => xm - ym - zm
  
  do i2=1,nx2
     do j2=1,ny2
        do k2=1,nz2
           ! Find the position for interpolation
           select case (trim(adjustl(dir)))
           case ('U')
              xd = x2(i2)
              yd = ym2(j2)
              zd = zm2(k2)
              x => x1;  nx = nx1+1
              y => ym1; ny = ny1
              z => zm1; nz = nz1
           case ('V')
              xd = xm2(i2)
              yd = y2(j2)
              zd = zm2(k2)
              x => xm1; nx = nx1
              y => y1;  ny = ny1+1
              z => zm1; nz = nz1
           case ('W')
              xd = xm2(i2)
              yd = ym2(j2)
              zd = z2(k2)
              x => xm1; nx = nx1
              y => ym1; ny = ny1
              z => z1;  nz = nz1+1
           case ('SC')
              xd = xm2(i2)
              yd = ym2(j2)
              zd = zm2(k2)
              x => xm1; nx = nx1
              y => ym1; ny = ny1
              z => zm1; nz = nz1
           end select

           ! Find the nearest points in mesh1
           call bisection(xd,i1,x,nx)
           call bisection(yd,j1,y,ny)
           call bisection(zd,k1,z,nz)
           i1 = max(1,min(i1,nx1-1))
           j1 = max(1,min(j1,ny1-1))
           k1 = max(1,min(k1,nz1-1))
           
           ! Interpolate the point
           wx1 = 1.0_WP
           wy1 = 1.0_WP
           wz1 = 1.0_WP
           if (nx.ne.1) wx1 = (x(i1+1)-xd   )/(x(i1+1)-x(i1))
           if (ny.ne.1) wy1 = (y(j1+1)-yd   )/(y(j1+1)-y(j1))
           if (nz.ne.1) wz1 = (z(k1+1)-zd   )/(z(k1+1)-z(k1))
           wx2 = 1.0_WP-wx1
           wy2 = 1.0_WP-wy1
           wz2 = 1.0_WP-wz1
           
           ! Combine the interpolation coefficients to form a tri-linear interpolation
           ww(0,0,0)=wx1*wy1*wz1
           ww(1,0,0)=wx2*wy1*wz1
           ww(0,1,0)=wx1*wy2*wz1
           ww(1,1,0)=wx2*wy2*wz1
           ww(0,0,1)=wx1*wy1*wz2
           ww(1,0,1)=wx2*wy1*wz2
           ww(0,1,1)=wx1*wy2*wz2
           ww(1,1,1)=wx2*wy2*wz2
           
           ! Perform the actual interpolation on A
           A2(i2,j2,k2) = 0.0_WP
           do ni=0,min(i1+1,nx1)-i1
              do nj=0,min(j1+1,ny1)-j1
                 do nk=0,min(k1+1,nz1)-k1
                    A2(i2,j2,k2) = A2(i2,j2,k2) + ww(ni,nj,nk)*A1(i1+ni,j1+nj,k1+nk)
                 end do
              end do
           end do

        end do
     end do
  end do
  
  return
end subroutine interp_data



! -------------------------------------------------- !
! Bisection routine                                  !
! Gets an array and its size as well as a position x !
! Returns the index between 1 and nx                 !
! x(iloc)<=xloc<x(iloc+1)                            !
! Assuming xarray is monotonically increasing        !
! -------------------------------------------------- !
subroutine bisection(xloc,iloc,x,nx)
  use precision
  implicit none
  
  real(WP), intent(in)  :: xloc
  integer,  intent(out) :: iloc
  integer,  intent(in)  :: nx
  real(WP), dimension(1:nx), intent(in) :: x
  
  integer :: il,im,iu
  
  ! Take care of outside points
  if (xloc.lt.x(1)) then
     iloc = 1
  else if (xloc.ge.x(nx)) then
     iloc = nx-1
  else
     
     ! Initialize lower and upper limits
     il=1
     iu=nx
     
     ! While not done
     do while (iu-il.gt.1)
        ! Compute a mid-point
        im=(iu+il)/2
        ! Replace lower of upper limit as appropriate
        if (xloc.ge.x(im)) then
           il=im
        else
           iu=im
        end if
     end do
     
     ! Return
     iloc = il
  end if
  
  return
end subroutine bisection
