module fourier
  use pressure
  implicit none
  include 'fftw3.f'

  ! Use of fft
  logical :: use_fft

  ! Oddball
  logical :: oddball

  ! FFTW variables
  real(WP), dimension(:), pointer :: in_x,in_y,in_z
  real(WP), dimension(:), pointer :: out_x,out_y,out_z
  integer(KIND=8) :: fplan_x,fplan_y,fplan_z
  integer(KIND=8) :: bplan_x,bplan_y,bplan_z

  !$OMP THREADPRIVATE(in_x,in_y,in_z,out_x,out_y,out_z)
  !$OMP THREADPRIVATE(fplan_x,fplan_y,fplan_z)
  !$OMP THREADPRIVATE(bplan_x,bplan_y,bplan_z)

end module fourier


! ============================= !
! Initialize the fourier module !
! ============================= !
subroutine fourier_init
  use fourier
  use parser
  use parallel
  implicit none

  call parser_read("Pressure fft",use_fft,.false.)
  
  ! Nothing to do if no fourier
  if (.not.use_fft) return
  
  ! Check directions
  fft_x = .false.
  fft_y = .false.
  fft_z = .false.
  if (npx.eq.1 .and. xper.eq.1 .and. nx.ne.1) fft_x = .true.
  if (npy.eq.1 .and. yper.eq.1 .and. ny.ne.1) fft_y = .true.
  if (npz.eq.1 .and. zper.eq.1 .and. nz.ne.1) fft_z = .true.
  use_fft = fft_x .or. fft_y .or. fft_z
  
  ! Nothing to do if no fourier
  if (.not.use_fft) return

  ! Test if not even number of points
  if (fft_x .and. mod(nx,2).ne.0) call die('fourier_init: fft only with even number of points')
  if (fft_y .and. mod(ny,2).ne.0) call die('fourier_init: fft only with even number of points')
  if (fft_z .and. mod(nz,2).ne.0) call die('fourier_init: fft only with even number of points')
  
  !$OMP PARALLEL

  ! Create plan - X
  if (fft_x) then
     allocate(in_x (nx))
     allocate(out_x(nx))
     !$OMP CRITICAL
     call dfftw_plan_r2r_1d(fplan_x,nx,in_x,out_x,FFTW_R2HC,FFTW_MEASURE)
     call dfftw_plan_r2r_1d(bplan_x,nx,in_x,out_x,FFTW_HC2R,FFTW_MEASURE)
     !$OMP END CRITICAL
  end if
  
  ! Create plan - Y
  if (fft_y) then
     allocate(in_y (ny))
     allocate(out_y(ny))
     !$OMP CRITICAL
     call dfftw_plan_r2r_1d(fplan_y,ny,in_y,out_y,FFTW_R2HC,FFTW_MEASURE)
     call dfftw_plan_r2r_1d(bplan_y,ny,in_y,out_y,FFTW_HC2R,FFTW_MEASURE)
     !$OMP END CRITICAL
  end if
  
  ! Create plan - Z
  if (fft_z) then
     allocate(in_z (nz))
     allocate(out_z(nz))
     !$OMP CRITICAL
     call dfftw_plan_r2r_1d(fplan_z,nz,in_z,out_z,FFTW_R2HC,FFTW_MEASURE)
     call dfftw_plan_r2r_1d(bplan_z,nz,in_z,out_z,FFTW_HC2R,FFTW_MEASURE)
     !$OMP END CRITICAL
  end if

  !$OMP END PARALLEL
  
  ! Oddball
  oddball = .false.
  if ( (fft_x .or. nx.eq.1) .and. &
       (fft_y .or. ny.eq.1) .and. &
       (fft_z .or. nz.eq.1) ) oddball = .true.
  
  return
end subroutine fourier_init


! ==================== !
! Update the laplacian !
! ==================== !
subroutine fourier_operator
  use fourier
  use parallel
  use metric_velocity_conv
  use math
  use masks
  implicit none
  
  real(WP) :: k2eff,kk,dk
  integer  :: i,j,k,st,k0,n
  complex(WP), parameter :: ii = (0.0_WP,1.0_WP)

  ! Nothing to do if no fourier
  if (.not.use_fft) return
  
  !$OMP PARALLEL PRIVATE(dk,kk,k2eff,k0)

  ! Compute the keff^2 - X
  if (fft_x) then
     dk = twoPi/xL
     !$OMP DO
     do i=imin_,imax_
        kk = real(i-imin,WP)*dk
        if (i-imin.gt.(nx/2)) kk = real(nx-i+imin,WP)*dk
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              k2eff = real(sum(lap(i,j,k,1,:)*exp(ii*kk*xm(i-stcp:i+stcp)))/exp(ii*kk*xm(i)),WP)
              lap(i,j,k,1,:) = 0.0_WP
              lap(i,j,k,1,0) = k2eff
           end do
        end do
     end do
     !$OMP END DO
  end if

  ! Compute the keff^2 - Y
  if (fft_y) then
     dk = twoPi/yL
     !$OMP DO
     do j=jmin_,jmax_
        kk = real(j-jmin,WP)*dk
        if (j-jmin.gt.(ny/2)) kk = real(ny-j+jmin,WP)*dk
        do k=kmin_,kmax_
           do i=imin_,imax_
              k2eff = real(sum(lap(i,j,k,2,:)*exp(ii*kk*ym(j-stcp:j+stcp)))/exp(ii*kk*ym(j)),WP)
              lap(i,j,k,2,:) = 0.0_WP
              lap(i,j,k,2,0) = k2eff
           end do
        end do
     end do
     !$OMP END DO
  end if

  ! Compute the keff^2 - Z
  if (fft_z) then
     dk = twoPi/zL
     !$OMP DO
     do k=kmin_,kmax_
        kk = real(k-kmin,WP)*dk
        if (k-kmin.gt.(nz/2)) kk = real(nz-k+kmin,WP)*dk
        do j=jmin_,jmax_
           do i=imin_,imax_
              k2eff = real(sum(lap(i,j,k,3,:)*exp(ii*kk*zm(k-stcp:k+stcp)))/exp(ii*kk*zm(k)),WP)
              lap(i,j,k,3,:) = 0.0_WP
              lap(i,j,k,3,0) = k2eff
           end do
        end do
     end do
     !$OMP END DO
     ! Point on the other side of the axis
     if (icyl.eq.1 .and. jproc.eq.1) then
        !$OMP DO
        do k=kmin_,kmax_
           k0 = k-kmin
           if (k0.gt.(nz/2)) k0 = nz-k+kmin
           do i=imin_,imax_
              if (mask(i,jmin).ne.1) then
                 
                 do st=0,stcp-1
                    do n=-stcp,-st-1
                       lap(i,jmin+st,k,2,-n-2*st-1) = lap(i,jmin+st,k,2,-n-2*st-1) + lap(i,jmin+st,k,2,n)*(-1.0_WP)**k0
                       lap(i,jmin+st,k,2,n) = 0.0_WP
                    end do
                 end do
                 
              end if
           end do
        end do
        !$OMP END DO
     end if
  end if

  !$OMP END PARALLEL

  ! Necessary communications done in pressure_rescale
  
  ! Oddball
  if (oddball) lap(imin_,jmin_,kmin_,1,0) = 1.0_WP
  
  return
end subroutine fourier_operator


! ========================================= !
! PreProcess the DP array before the solver !
! ========================================= !
subroutine fourier_rhs
  use fourier
  implicit none

  integer :: i,j,k
  
  ! Nothing to do if no fourier
  if (.not.use_fft) return

  !$OMP PARALLEL

  ! Forward transform - X
  if (fft_x) then
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           in_x = RP(:,j,k)
           call dfftw_execute(fplan_x)
           RP(:,j,k) = out_x
        end do
     end do
     !$OMP END DO
  end if
  
  ! Forward transform - Y
  if (fft_y) then
     !$OMP DO
     do k=kmin_,kmax_
        do i=imin_,imax_
           in_y = RP(i,:,k)
           call dfftw_execute(fplan_y)
           RP(i,:,k) = out_y
        end do
     end do
     !$OMP END DO
  end if
  
  ! Forward transform - Z
  if (fft_z) then
     !$OMP DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           in_z = RP(i,j,:)
           call dfftw_execute(fplan_z)
           RP(i,j,:) = out_z
        end do
     end do
     !$OMP END DO
  end if

  !$OMP END PARALLEL

  ! Oddball
  if (oddball) RP(imin_,jmin_,kmin_) = 0.0_WP
  
  return
end subroutine fourier_rhs


! ========================================= !
! PostProcess the DP array after the solver !
! ========================================= !
subroutine fourier_dp
  use fourier
  implicit none

  integer :: i,j,k

  ! Nothing to do if no fourier
  if (.not.use_fft) return
  
  !$OMP PARALLEL

  ! Backward transform - X
  if (fft_x) then
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           in_x = DP(imin_:imax_,j,k)
           call dfftw_execute(bplan_x)
           DP(imin_:imax_,j,k) = out_x/real(nx,WP)
        end do
     end do
     !$OMP END DO
  end if
  
  ! Backward transform - Y
  if (fft_y) then
     !$OMP DO
     do k=kmin_,kmax_
        do i=imin_,imax_
           in_y = DP(i,jmin_:jmax_,k)
           call dfftw_execute(bplan_y)
           DP(i,jmin_:jmax_,k) = out_y/real(ny,WP)
        end do
     end do
     !$OMP END DO
  end if
  
  ! Backward transform - Z
  if (fft_z) then
     !$OMP DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           in_z = DP(i,j,kmin_:kmax_)
           call dfftw_execute(bplan_z)
           DP(i,j,kmin_:kmax_) = out_z/real(nz,WP)
        end do
     end do
     !$OMP END DO
  end if

  !$OMP END PARALLEL
  
  return
end subroutine fourier_dp

