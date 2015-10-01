module filter
  use geometry
  use partition
  use parallel
  implicit none
  
  ! Temporary work arrays
  real(WP), dimension(:,:,:), pointer :: buffer
  
  ! Filter - 3D
  real(WP), dimension(:,:,:), pointer :: filterX_3Dn
  real(WP), dimension(:,:,:), pointer :: filterY_3Dn
  real(WP), dimension(:,:,:), pointer :: filterX_3Dd
  real(WP), dimension(:,:,:), pointer :: filterY_3Dd
  real(WP), dimension(:),     pointer :: filterZ_3D
  real(WP), dimension(:,:),   pointer :: delta_3D
  real(WP), dimension(:,:),   pointer :: ratio_3D
  
  ! Temporary arrays for filtering
  real(WP), dimension(-1:+1,-1:+1) :: tmp2D
  real(WP), dimension(-1:+1) :: tmp1D  

  !$OMP THREADPRIVATE(tmp1D,tmp2D)
  
end module filter


! ---------------------- !
! Filters initialization !
! ---------------------- !
subroutine filter_init
  use filter
  use masks
  implicit none
  
  real(WP) :: dx1,dx2,dy1,dy2
  integer :: i,j
  
  ! Allocate work array
  allocate(buffer(imino_:imaxo_,jmino_:max(jmaxo_,jmino_+5),kmino_:kmaxo_))
  

  ! Compute 3D filters
  ! ------------------

  ! Allocate filter arrays - 3D filter
  allocate(filterX_3Dd(imin_:imax_,jmin_-1:jmax_+1,-1:+1))
  allocate(filterY_3Dd(imin_:imax_,jmin_-1:jmax_+1,-1:+1))
  allocate(filterX_3Dn(imin_:imax_,jmin_-1:jmax_+1,-1:+1))
  allocate(filterY_3Dn(imin_:imax_,jmin_-1:jmax_+1,-1:+1))
  allocate(filterZ_3D (-1:+1))
  allocate(delta_3D(imino_:imaxo_,jmino_:jmaxo_))
  allocate(ratio_3D(imin_ :imax_, jmin_ :jmax_ ))
  
  ! Simpsons Rule for filter coefficients
  filterZ_3D(-1) = 1.0_WP/6.0_WP
  filterZ_3D( 0) = 2.0_WP/3.0_WP
  filterZ_3D(+1) = 1.0_WP/6.0_WP
  do i=imin_,imax_
     dx1 = xm(i+1) - xm(i)
     dx2 = xm(i)   - xm(i-1)
     filterX_3Dd(i,:,-1) = (-0.5_WP*(dx1-dx2)-(dx1-dx2)**2/(6.0_WP*dx2)+dx1/3.0_WP)/(dx1+dx2)
     filterX_3Dd(i,:, 0) = 2.0_WP/3.0_WP + (dx1-dx2)**2/(6.0_WP*dx1*dx2)
     filterX_3Dd(i,:,+1) = ( 0.5_WP*(dx1-dx2)-(dx1-dx2)**2/(6.0_WP*dx1)+dx2/3.0_WP)/(dx1+dx2)
     filterX_3Dn(i,:,-1) = (-0.5_WP*(dx1-dx2)-(dx1-dx2)**2/(6.0_WP*dx2)+dx1/3.0_WP)/(dx1+dx2)
     filterX_3Dn(i,:, 0) = 2.0_WP/3.0_WP + (dx1-dx2)**2/(6.0_WP*dx1*dx2)
     filterX_3Dn(i,:,+1) = ( 0.5_WP*(dx1-dx2)-(dx1-dx2)**2/(6.0_WP*dx1)+dx2/3.0_WP)/(dx1+dx2)
  end do
  do j=max(jmin_-1,jmin),min(jmax_+1,jmax)
     dy1 = ym(j+1) - ym(j)
     dy2 = ym(j)   - ym(j-1)
     filterY_3Dd(:,j,-1) = (-0.5_WP*(dy1-dy2)-(dy1-dy2)**2/(6.0_WP*dy2)+dy1/3.0_WP)/(dy1+dy2)
     filterY_3Dd(:,j, 0) = 2.0_WP/3.0_WP + (dy1-dy2)**2/(6.0_WP*dy1*dy2)
     filterY_3Dd(:,j,+1) = ( 0.5_WP*(dy1-dy2)-(dy1-dy2)**2/(6.0_WP*dy1)+dy2/3.0_WP)/(dy1+dy2)
     filterY_3Dn(:,j,-1) = (-0.5_WP*(dy1-dy2)-(dy1-dy2)**2/(6.0_WP*dy2)+dy1/3.0_WP)/(dy1+dy2)
     filterY_3Dn(:,j, 0) = 2.0_WP/3.0_WP + (dy1-dy2)**2/(6.0_WP*dy1*dy2)
     filterY_3Dn(:,j,+1) = ( 0.5_WP*(dy1-dy2)-(dy1-dy2)**2/(6.0_WP*dy1)+dy2/3.0_WP)/(dy1+dy2)
  end do
  ! Delta and filter ratio
  do j=jmino_,jmaxo_
     do i=imino_,imaxo_
        if (icyl.eq.0) then
           delta_3D(i,j) = (dz*dy(j)*dx(i))**(1.0_WP/3.0_WP)
        else
           delta_3D(i,j) = (dz*abs(ym(j))*dy(j)*dx(i))**(1.0_WP/3.0_WP)
        end if
     end do
  end do
  do j=jmin_,jmax_
     do i=imin_,imax_
        ratio_3D(i,j) = (2.0_WP*(xm(i+1)-xm(i-1))/(x(i+1)-x(i))*(ym(j+1)-ym(j-1))/(y(j+1)-y(j)))**(1.0_WP/3.0_WP)
     end do
  end do


  ! Account for walls
  ! -----------------

  ! Dirichlet
  do j = max(jmin_-1,jmin),min(jmax_+1,jmax)
     do i = imin_,imax_
        if (mask(i-1,j).eq.1) then
           filterX_3Dd(i,j, 0) = filterX_3Dd(i,j,0) - filterX_3Dd(i,j,-1)
           filterX_3Dd(i,j,-1) = 0.0_WP
        end if
        if (mask(i+1,j).eq.1) then
           filterX_3Dd(i,j, 0) = filterX_3Dd(i,j,0) - filterX_3Dd(i,j,+1)
           filterX_3Dd(i,j,+1) = 0.0_WP
        end if
        if (mask(i,j-1).eq.1) then
           filterY_3Dd(i,j, 0) = filterY_3Dd(i,j,0) - filterY_3Dd(i,j,-1)
           filterY_3Dd(i,j,-1) = 0.0_WP
        end if
        if (mask(i,j+1).eq.1) then
           filterY_3Dd(i,j, 0) = filterY_3Dd(i,j,0) - filterY_3Dd(i,j,+1)
           filterY_3Dd(i,j,+1) = 0.0_WP
        end if
        if (mask(i  ,j).eq.1) then
           filterX_3Dd(i,j,: ) = 0.0_WP
           filterX_3Dd(i,j,0 ) = 1.0_WP 
           filterY_3Dd(i,j,: ) = 0.0_WP
           filterY_3Dd(i,j,0 ) = 1.0_WP
        end if
     end do
  end do
  
  ! Neumann
  do j = max(jmin_-1,jmin),min(jmax_+1,jmax)
     do i = imin_,imax_
        if (mask(i-1,j).eq.1) then
           filterX_3Dn(i,j, 0) = filterX_3Dn(i,j,0) + filterX_3Dn(i,j,-1)
           filterX_3Dn(i,j,-1) = 0.0_WP
        end if
        if (mask(i+1,j).eq.1) then
           filterX_3Dn(i,j, 0) = filterX_3Dn(i,j,0) + filterX_3Dn(i,j,+1)
           filterX_3Dn(i,j,+1) = 0.0_WP
        end if
        if (mask(i,j-1).eq.1) then
           filterY_3Dn(i,j, 0) = filterY_3Dn(i,j,0) + filterY_3Dn(i,j,-1)
           filterY_3Dn(i,j,-1) = 0.0_WP
        end if
        if (mask(i,j+1).eq.1) then
           filterY_3Dn(i,j, 0) = filterY_3Dn(i,j,0) + filterY_3Dn(i,j,+1)
           filterY_3Dn(i,j,+1) = 0.0_WP
        end if
        if (mask(i  ,j).eq.1) then
           filterX_3Dn(i,j,: ) = 0.0_WP
           filterX_3Dn(i,j,0 ) = 1.0_WP 
           filterY_3Dn(i,j,: ) = 0.0_WP
           filterY_3Dn(i,j,0 ) = 1.0_WP
        end if
     end do
  end do
  
  ! Prolongate the arrays
  if (jproc.eq.1) then
     filterY_3Dn(:,jmin-1,:) = filterY_3Dn(:,jmin,:)
     filterY_3Dd(:,jmin-1,:) = filterY_3Dd(:,jmin,:)
  end if
  if (jproc.eq.npy) then
     filterY_3Dn(:,jmax+1,:) = filterY_3Dn(:,jmax,:)
     filterY_3Dd(:,jmax+1,:) = filterY_3Dd(:,jmax,:)
  end if
  
  
  ! Enforce the physical boundary conditions of the domain
  ! ------------------------------------------------------
  if (xper.ne.1) then
     ! Left boundary
     ! -> Neumann
     if (iproc.eq.1) then
        filterX_3Dn(imin,:, 0) = filterX_3Dn(imin,:,0) + filterX_3Dn(imin,:,-1)
        filterX_3Dn(imin,:,-1) = 0.0_WP
        filterX_3Dd(imin,:, 0) = filterX_3Dd(imin,:,0) + filterX_3Dd(imin,:,-1)
        filterX_3Dd(imin,:,-1) = 0.0_WP
     end if
     
     ! Right boundary
     ! -> Neumann
     if (iproc.eq.npx) then
        filterX_3Dn(imax,:, 0) = filterX_3Dn(imax,:,0) + filterX_3Dn(imax,:,+1)
        filterX_3Dn(imax,:,+1) = 0.0_WP
        filterX_3Dd(imax,:, 0) = filterX_3Dd(imax,:,0) + filterX_3Dd(imax,:,+1)
        filterX_3Dd(imax,:,+1) = 0.0_WP
     end if
  end if
     
  if (yper.ne.1) then
     ! Lower boundary
     ! -> Neumann
     if (jproc.eq.1) then
        filterY_3Dn(:,jmin, 0) = filterY_3Dn(:,jmin,0) + filterY_3Dn(:,jmin,-1)
        filterY_3Dn(:,jmin,-1) = 0.0_WP
        filterY_3Dd(:,jmin, 0) = filterY_3Dd(:,jmin,0) + filterY_3Dd(:,jmin,-1)
        filterY_3Dd(:,jmin,-1) = 0.0_WP
     end if
     ! Upper boundary
     ! -> Neumann
     if (jproc.eq.npy) then
        filterY_3Dn(:,jmax, 0) = filterY_3Dn(:,jmax,0) + filterY_3Dn(:,jmax,+1)
        filterY_3Dn(:,jmax,+1) = 0.0_WP
        filterY_3Dd(:,jmax, 0) = filterY_3Dd(:,jmax,0) + filterY_3Dd(:,jmax,+1)
        filterY_3Dd(:,jmax,+1) = 0.0_WP
     end if
  end if

  return
end subroutine filter_init


! -------------------------------------------- !
! Filter only locally given a box of 27 points !
! -------------------------------------------- !
subroutine filter_local_3D(A,B,i,j,bc)
  use filter
  implicit none
  
  real(WP), dimension(-1:+1,-1:+1,-1:+1), intent(in) :: A
  real(WP), intent(out) :: B  
  integer, intent(in) :: i,j
  integer :: nj,nk
  character(len=1) :: bc
  
  ! Do not thread this subroutine

  select case (bc)
  case ('n')
     do nk=-1,+1
        do nj=-1,+1
           tmp2D (nj,nk) = sum(filterX_3Dn(i,j+nj,:) * A(:,nj,nk))
        end do
     end do
     do nk=-1,+1
        tmp1D (nk) = sum(filterY_3Dn(i,j,:) * tmp2D(:,nk))
     end do
     B = sum(filterZ_3D * tmp1D)
  case ('d')
     do nk=-1,+1
        do nj=-1,+1
           tmp2D (nj,nk) = sum(filterX_3Dd(i,j+nj,:) * A(:,nj,nk))
        end do
     end do
     do nk=-1,+1
        tmp1D (nk) = sum(filterY_3Dd(i,j,:) * tmp2D(:,nk))
     end do
     B = sum(filterZ_3D * tmp1D)
  end select
  
  return
end subroutine filter_local_3D


! ---------------------------------------- !
! Filter globally using a 27 point stencil !
! ---------------------------------------- !
subroutine filter_global_3D(A,B,sym,bc)
  use filter
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: B
  integer :: i,j,k
  character(len=1), intent(in) :: sym
  character(len=1) :: bc
  
  ! Filter in x
  select case (bc)
  case ('n')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buffer(i,j,k) = filterX_3Dn(i,j,-1) * A(i-1,j,k) + &
                              filterX_3Dn(i,j, 0) * A(i  ,j,k) + &
                              filterX_3Dn(i,j,+1) * A(i+1,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  case ('d')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buffer(i,j,k) = filterX_3Dd(i,j,-1) * A(i-1,j,k) + &
                              filterX_3Dd(i,j, 0) * A(i  ,j,k) + &
                              filterX_3Dd(i,j,+1) * A(i+1,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end select

  ! Update ghost cells
  call boundary_update_border(buffer,sym,'ym')
  
  ! Non periodic directions
  call boundary_neumann(buffer,'+ym')
  call boundary_neumann(buffer,'-ym')
  
  ! Filter in y
  select case (bc)
  case ('n')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              B(i,j,k) = filterY_3Dn(i,j,-1) * buffer(i,j-1,k) + &
                         filterY_3Dn(i,j, 0) * buffer(i,j  ,k) + &
                         filterY_3Dn(i,j,+1) * buffer(i,j+1,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  case ('d')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              B(i,j,k) = filterY_3Dd(i,j,-1) * buffer(i,j-1,k) + &
                         filterY_3Dd(i,j, 0) * buffer(i,j  ,k) + &
                         filterY_3Dd(i,j,+1) * buffer(i,j+1,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end select

  ! Update ghost cells
  call boundary_update_border(B,sym,'ym')
  
  ! Non periodic directions
  call boundary_neumann(B,'+ym')
  call boundary_neumann(B,'-ym')

  !$OMP PARALLEL

  ! Filter in z
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buffer(i,j,k) = filterZ_3D(-1) * B(i,j,k-1) + &
                           filterZ_3D( 0) * B(i,j,k  ) + &
                           filterZ_3D(+1) * B(i,j,k+1)
        end do
     end do
  end do
  !$OMP END DO

  ! Put it back
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           B(i,j,k) = buffer(i,j,k)
        end do
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL
  
  ! Update ghost cells
  call boundary_update_border(B,sym,'ym')
  
  ! Non periodic directions
  call boundary_neumann(B,'+ym')
  call boundary_neumann(B,'-ym')

  ! Case of the inlet conditions
  if (iproc.eq.1 .and. xper.eq.0) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino,imin-1
              B(i,j,k) = A(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Case of the outlet conditions
  if (iproc.eq.npx .and. xper.eq.0) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imax+1,imaxo
              B(i,j,k) = A(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine filter_global_3D


