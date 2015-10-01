module channel2jet
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
  integer :: nvar1,nvar2
  real(WP) :: dt,time
  character(len=str_short), dimension(:), pointer :: names1,names2
  character(len=str_medium) :: config1
  character(len=str_medium) :: config2
  integer, dimension(:,:), pointer :: mask1
  integer, dimension(:,:), pointer :: mask2
  
  ! Data
  real(WP), dimension(:,:,:), pointer :: data1
  real(WP), dimension(:,:,:), pointer :: data1_gc
  real(WP), dimension(:,:,:), pointer :: data2
  
  ! Mesh
  real(WP), dimension(:), pointer :: x1,xm1
  real(WP), dimension(:), pointer :: xmo1
  real(WP), dimension(:), pointer :: y1,ym1
  real(WP), dimension(:), pointer :: z1,zm1
  real(WP), dimension(:), pointer :: zmo1
  real(WP), dimension(:), pointer :: x2,xm2
  real(WP), dimension(:), pointer :: y2,ym2
  real(WP), dimension(:), pointer :: z2,zm2
  
end module channel2jet


program c2j
  use channel2jet
  implicit none
  
  character(len=str_medium) :: fconfig1
  character(len=str_medium) :: fconfig2
  character(len=str_medium) :: fdata1
  character(len=str_medium) :: fdata2
  integer :: iunit,iunit1,iunit2,var
  integer :: i,j,k,ierr
  
  ! Get file names from user
  print*,'===================================='
  print*,'| ARTS - Channel to jet conversion |'
  print*,'===================================='
  print*
  print "(a28,$)", " Channel config file : "
  read "(a)", fconfig1
  print "(a33,$)", " Jet config file : "
  read "(a)", fconfig2
  print "(a26,$)", " Channel data file : "
  read "(a)", fdata1
  print "(a31,$)", " Jet data file : "
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
  print*,'Channel file'
  print*,'Config : ',config1
  print*,'Grid :',nx1,'x',ny1,'x',nz1
  allocate(x1(nx1+1),y1(ny1+1),z1(nz1+1))
  allocate(mask1(nx1,ny1))
  call BINARY_FILE_READ(iunit,x1,nx1+1,kind(x1),ierr)
  call BINARY_FILE_READ(iunit,y1,ny1+1,kind(y1),ierr)
  call BINARY_FILE_READ(iunit,z1,nz1+1,kind(z1),ierr)
  call BINARY_FILE_READ(iunit,mask1,nx1*ny1,kind(mask1),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  print*,'Domain x :',x1(1),'x',x1(nx1+1)
  print*,'Domain y :',y1(1),'x',y1(ny1+1)
  print*,'Domain z :',z1(1),'x',z1(nz1+1)
  
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
  print*,'Jet file'
  print*,'Config : ',config2
  print*,'Grid :',nx2,'x',ny2,'x',nz2
  allocate(x2(nx2+1),y2(ny2+1),z2(nz2+1))
  allocate(mask2(nx2,ny2))
  call BINARY_FILE_READ(iunit,x2,nx2+1,kind(x2),ierr)
  call BINARY_FILE_READ(iunit,y2,ny2+1,kind(y2),ierr)
  call BINARY_FILE_READ(iunit,z2,nz2+1,kind(z2),ierr)
  call BINARY_FILE_READ(iunit,mask2,nx2*ny2,kind(mask2),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  print*,'Domain x :',x2(1),'x',x2(nx2+1)
  print*,'Domain y :',y2(1),'x',y2(ny2+1)
  print*,'Domain z :',z2(1),'x',z2(nz2+1)
  
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
  
  ! Mesh with gc
  allocate(xmo1(0:nx1+1),zmo1(0:nz1+1))
  xmo1(1:nx1)=xm1;xmo1(0)=xmo1(1)-(xmo1(2)-xmo1(1));xmo1(nx1+1)=xmo1(nx1)+(xmo1(nx1)-xmo1(nx1-1))
  zmo1(1:nz1)=zm1;zmo1(0)=zmo1(1)-(zmo1(2)-zmo1(1));zmo1(nz1+1)=zmo1(nz1)+(zmo1(nz1)-zmo1(nz1-1))
  
  ! Read the source data file
  call BINARY_FILE_OPEN(iunit1,trim(fdata1),"r",ierr)
  call BINARY_FILE_READ(iunit1,nx1,  1,kind(nx1),  ierr)
  call BINARY_FILE_READ(iunit1,ny1,  1,kind(ny1),  ierr)
  call BINARY_FILE_READ(iunit1,nz1,  1,kind(nz1),  ierr)
  call BINARY_FILE_READ(iunit1,nvar1,1,kind(nvar1),ierr)
  call BINARY_FILE_READ(iunit1,dt,   1,kind(dt),   ierr)
  call BINARY_FILE_READ(iunit1,time, 1,kind(time), ierr)
  allocate(names1(nvar1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit1,names1(var),str_short,kind(names1),ierr)
  end do
  print*,'Variables in channel data file: ',names1
  
  ! Change the variables: replace P by G
  nvar2=nvar1
  allocate(names2(nvar2))
  do var=1,nvar2
     if (trim(names1(var)).eq.'P') then
        names2(var)='G'
     else
        names2(var)=names1(var)
     end if
  end do
  print*,'Variables in jet data file: ',names2
  
  ! Reset time
  dt=0.0_WP
  time=0.0_WP
  
  ! Dump the new data file
  call BINARY_FILE_OPEN(iunit2,trim(fdata2),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,nx2,  1,kind(nx2),  ierr)
  call BINARY_FILE_WRITE(iunit2,ny2,  1,kind(ny2),  ierr)
  call BINARY_FILE_WRITE(iunit2,nz2,  1,kind(nz2),  ierr)
  call BINARY_FILE_WRITE(iunit2,nvar2,1,kind(nvar2),ierr)
  call BINARY_FILE_WRITE(iunit2,dt,   1,kind(dt),   ierr)
  call BINARY_FILE_WRITE(iunit2,time, 1,kind(time), ierr)
  do var=1,nvar2
     call BINARY_FILE_WRITE(iunit2,names2(var),str_short,kind(names2),ierr)
  end do
  
  ! Prepare for interpolation
  allocate(data1(nx1,ny1,nz1))
  allocate(data1_gc(0:nx1+1,1:ny1,0:nz1+1))
  allocate(data2(nx2,ny2,nz2))
  
  ! Interpolate
  do var=1,nvar2
     ! Read variable
     call BINARY_FILE_READ(iunit1,data1,nx1*ny1*nz1,kind(data1),ierr)
     ! Communicate
     data1_gc(1:nx1,1:ny1,1:nz1)=data1
     data1_gc(nx1+1,:,:)=data1_gc( 1 ,:,:)
     data1_gc(  0  ,:,:)=data1_gc(nx1,:,:)
     data1_gc(:,:,nz1+1)=data1_gc(:,:, 1 )
     data1_gc(:,:,  0  )=data1_gc(:,:,nz1)
     print*,data1_gc(104,41,256),'merde'
     print*,data1(104,41,256),'merde'
     ! Perform interpolation
     select case (trim(adjustl(names2(var))))
     case ('U')
        call interp_data('U')
     case ('V')
        call interp_data('V')
     case ('W')
        call interp_data('W')
     case ('G')
        call create_Gfield(data2)
     end select
     ! Write variable
     call BINARY_FILE_WRITE(iunit2,data2,nx2*ny2*nz2,kind(data2),ierr)
     ! Dump status
     print*,trim(adjustl(names2(var))),' done'
  end do
  
  ! Close files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
end program c2j


! ------------------------------ !
! Add a G field to the data file !
! ------------------------------ !
subroutine create_Gfield(A)
  use channel2jet
  
  integer :: i,j,k
  real(WP), dimension(nx2,ny2,nz2), intent(inout) :: A
  
  ! Create planar liquid sheet
  do i=1,nx2
     do j=1,ny2
        do k=1,nz2
           A(i,j,k) = min(0.5_WP-ym2(j),0.5_WP+ym2(j))
        end do
     end do
  end do
  
  return
end subroutine create_Gfield


! --------------------------------------- !
! Interpolate the data between two arrays !
! For a given staggering 'dir'            !
! --------------------------------------- !
subroutine interp_data(dir)
  use channel2jet
  implicit none
  
  character(len=*) :: dir
  real(WP), dimension(:), pointer :: x,y,z
  integer  :: nx,ny,nz
  integer  :: i2,j2,k2
  integer  :: i1,j1,k1
  integer  :: ni,nj,nk
  real(WP) :: xd,yd,zd
  integer  :: imin,jmin,kmin
  integer  :: imax,jmax,kmax
  
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
              xd = x2 (i2)
              yd = ym2(j2)
              zd = zm2(k2)
              x => x1;   nx = nx1+1; imin=1; imax=nx1+1
              y => ym1;  ny = ny1  ; jmin=1; jmax=ny1
              z => zmo1; nz = nz1+2; kmin=0; kmax=nz1+1
           case ('V')
              xd = xm2(i2)
              yd = y2 (j2)
              zd = zm2(k2)
              x => xmo1; nx = nx1+2; imin=0; imax=nx1+1
              y => y1;   ny = ny1+1; jmin=1; jmax=ny1+1
              z => zmo1; nz = nz1+2; kmin=0; kmax=nz1+1
           case ('W')
              xd = xm2(i2)
              yd = ym2(j2)
              zd = z2 (k2)
              x => xmo1; nx = nx1+2; imin=0; imax=nx1+1
              y => ym1;  ny = ny1  ; jmin=1; jmax=ny1
              z => z1;   nz = nz1+1; kmin=1; kmax=nz1+1
           end select
           
           ! Treat walls
           if (yd.le.-0.5_WP.or.yd.ge.0.5_WP) then
              data2(i2,j2,k2)=0.0_WP
              cycle
           end if
           
           ! Find the nearest points in mesh1
           call bisection(xd,i1,x,imin,imax)
           call bisection(yd,j1,y,jmin,jmax)
           call bisection(zd,k1,z,kmin,kmax)
           
           ! Interpolate the point
           wx1 = 1.0_WP
           wy1 = 1.0_WP
           wz1 = 1.0_WP
           if (nx.ne.1) wx1 = (x(i1+1)-xd)/(x(i1+1)-x(i1))
           if (ny.ne.1) wy1 = (y(j1+1)-yd)/(y(j1+1)-y(j1))
           if (nz.ne.1) wz1 = (z(k1+1)-zd)/(z(k1+1)-z(k1))
           
           ! Compute the other coefficients
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
           data2(i2,j2,k2) = 0.0_WP
           do ni=0,1!min(i1+1,nx)-i1
              do nj=0,1!min(j1+1,ny)-j1
                 do nk=0,1!min(k1+1,nz)-k1
                    data2(i2,j2,k2) = data2(i2,j2,k2) + ww(ni,nj,nk)*data1_gc(i1+ni,j1+nj,k1+nk)
                 end do
              end do
           end do
           
           if (abs(data2(i2,j2,k2)).gt.100.0_WP) then
              print*,'problem',data2(i2,j2,k2)
              print*,'at',i2,j2,k2
              print*,'interp',i1,j1,k1
              print*,'w',wx1,wy1,wz1
              print*,wx2,wy2,wz2
              print*,'A1',data1_gc(i1-1,j1  ,k1  ),data1_gc(i1,j1  ,k1  ),data1_gc(i1+1,j1  ,k1  )
              print*,'A1',data1_gc(i1-1,j1-1,k1  ),data1_gc(i1,j1-1,k1  ),data1_gc(i1+1,j1-1,k1  )
              print*,'A1',data1_gc(i1-1,j1+1,k1  ),data1_gc(i1,j1+1,k1  ),data1_gc(i1+1,j1+1,k1  )
              print*,'A1',data1_gc(i1-1,j1  ,k1-1),data1_gc(i1,j1  ,k1-1),data1_gc(i1+1,j1  ,k1-1)
              print*,'A1',data1_gc(i1-1,j1-1,k1-1),data1_gc(i1,j1-1,k1-1),data1_gc(i1+1,j1-1,k1-1)
              print*,'A1',data1_gc(i1-1,j1+1,k1-1),data1_gc(i1,j1+1,k1-1),data1_gc(i1+1,j1+1,k1-1)
              print*,'A1',data1_gc(i1-1,j1  ,k1+1),data1_gc(i1,j1  ,k1+1),data1_gc(i1+1,j1  ,k1+1)
              print*,'A1',data1_gc(i1-1,j1-1,k1+1),data1_gc(i1,j1-1,k1+1),data1_gc(i1+1,j1-1,k1+1)
              print*,'A1',data1_gc(i1-1,j1+1,k1+1),data1_gc(i1,j1+1,k1+1),data1_gc(i1+1,j1+1,k1+1)
           end if
           
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
subroutine bisection(xloc,iloc,x,indmin,indmax)
  use precision
  implicit none
  
  real(WP), intent(in)  :: xloc
  integer,  intent(out) :: iloc
  integer,  intent(in)  :: indmin,indmax
  real(WP), dimension(indmin:indmax), intent(in) :: x
  
  integer :: il,im,iu
  
  ! Take care of outside points
  if (xloc.lt.x(indmin)) then
     iloc = indmin
  else if (xloc.ge.x(indmax)) then
     iloc = indmax-1
  else
     
     ! Initialize lower and upper limits
     il=indmin
     iu=indmax-1
     
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
