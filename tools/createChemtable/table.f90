module table
  use flamelet
  use precision
  use string
  implicit none
  
  ! Table parameters
  integer :: nZVar, nZMean, n3
  real(WP), dimension(:), pointer :: ZVar, ZMean, Z3
  
  ! Beta pdf
  real(WP), dimension(:), pointer :: pdf
  
  ! List of Flamelet files
  integer :: nfiles
  character(len=str_long), dimension(:), pointer :: files
  
  ! Variables to be wrtitten into the table
  integer :: nvar_out
  character(len=str_medium), dimension(:), pointer :: output_name
  real(WP), dimension(:,:,:,:), pointer :: output_data
  integer, dimension(:,:,:), pointer :: mask
  
  ! Variable after convolution
  real(WP), dimension(:,:,:,:), pointer :: postconv
  
contains
  
  ! One third of the points between 0 and Zst : uniform mesh
  ! Two third between Zst and 1 : linear growth
  subroutine create_zmean
    use parser
    implicit none
    
    integer :: zm, zcut
    real(WP) :: m11,m12,m21,m22,r1,r2,delta
    real(WP) :: a,b,c
    real(WP) :: Zst,dz
    
    call parser_read('Stoichiometric mixture fraction',Zst)
    
    ! Create such mesh only if Zst sufficiently small
    if (Zst.gt.0.30_WP) then
       dz = 1.0_WP / real(nZMean-1,WP)
       do zm=1,nZMean
          ZMean(zm) = real(zm-1,WP) * dz
       end do
       return
    end if
    
    ! Linear mesh for [0,Zst]
    zcut = nZMean/3+1
    dz = Zst / real(zcut-1,WP)
    do zm=1,zcut
       ZMean(zm) = real(zm-1,WP) * dz
    end do
    
    ! Mesh with linear growth to reach Z=1
    m11 = real(nZMean**2-zcut**2,WP)
    m12 = real(nZMean-zcut,WP)
    m21 = real(2*zcut+1,WP)
    m22 = real(1,WP)
    r1 = 1.0_WP - Zst
    r2 = dz
    delta = m11*m22-m12*m21
    a = (+ m22*r1 - m12*r2 )/delta
    b = (- m21*r1 + m11*r2 )/delta
    c = Zst - a*zcut**2-b*zcut
    
    do zm=zcut+1,nZMean
       ZMean(zm) = a*real(zm,WP)**2 + b*real(zm,WP) + c
    end do
    
    return
  end subroutine create_zmean
  
  
  ! Uniform mesh between 0 and 0.25
  subroutine create_zvar
    implicit none
    
    integer :: zv
    
    do zv=1,nZVar
       !ZVar(zv) = 0.25_WP* real(zv-1,WP) / real(nZVar-1,WP)
       ZVar(zv) = 0.25_WP* (real(zv-1,WP) / real(nZVar-1,WP))**2
    end do
    
    return
  end subroutine create_zvar
  
  
  ! Precompute the PDF for a beta distribution with
  ! -> mean mixture fraction : zm
  ! -> variance of mixture fraction : zv
  subroutine create_beta_pdf(zm,zv)
    use math
    implicit none
    
    real(WP), intent(in) :: zm, zv
    real(WP) :: a,b,factor,tmp,dz,mean
    integer :: index1, n
    
    pdf = 0.0_WP
    
    ! Zero mean : delta at Z=0
    if (zm.le.1.0E-10_WP) then
       pdf(1) = 1.0_WP
       return
    end if
    
    ! Max mean : delta at Z=1
    if (zm.ge.1.0_WP-1.0E-10_WP) then
       pdf(nPoints) = 1.0_WP
       return
    end if
    
    ! Zero variance : delta at Z=zm
    if (zv.le.1.0E-10_WP) then
       index1 = 1
       do while (Z(index1).lt.zm) 
          index1 = index1+1
       end do
       pdf(index1-1) = (Z(index1)-zm)  /(Z(index1)-Z(index1-1))
       pdf(index1)   = (zm-Z(index1-1))/(Z(index1)-Z(index1-1))
       return
    end if
        
    ! Impossible cases => two delta at 0 and 1
    if (zv.ge.zm*(1.0_WP-zm)) then
       pdf(1) = 1.0_WP-zm
       pdf(nPoints) = zm
       return
    end if
    
    a = zm*(zm*(1.0_WP-zm)/zv - 1.0_WP)
    b = a/zm - a
    factor = gammaln(a+b) - gammaln(a) - gammaln(b)
    
    ! Left BC : explicit integration
    dz = 0.5_WP*(Z(2)-Z(1))
    tmp = a*log(dz) + factor
    pdf(1) = exp(tmp) / a
    ! Right BC : explicit integration
    dz = 0.5_WP*(Z(nPoints)-Z(nPoints-1))
    tmp = b*log(dz) + factor
    pdf(nPoints) = exp(tmp) / b
    ! Other Points
    do n=2,nPoints-1
       dz = 0.5_WP*(Z(n+1)-Z(n-1))
       tmp = (a-1.0_WP)*log(Z(n)) + (b-1.0_WP)*log(1.0_WP-Z(n))
       tmp = tmp + factor
       pdf(n) = exp(tmp) * dz
    end do
    
    ! Normalize the pdf
    pdf = pdf / sum(pdf)
    
    ! Check mean
    !mean = sum(pdf*Z)
    !pdf(nPoints) = pdf(nPoints) + (zm-mean)
    !pdf(1) = pdf(1) - (zm-mean)
    
    return
  end subroutine create_beta_pdf
  
end module table


subroutine table_init
  use table
  use parser
  implicit none
  
  ! Read the dimension of the final table
  call parser_read('Number of points for mean Z', nZMean)
  call parser_read('Number of points for variance of Z', nZVar)
  allocate(ZMean(nZMean))
  allocate(ZVar(nZVar))
  if (trim(combModel).eq.'Enthalpy Flamelet') then
     n3 = nfiles
  else
     call parser_read('Number of points for third direction', n3)
  end if
  allocate(Z3(n3))
  
  ! Create the first two directions of the table
  call create_zmean
  call create_zvar
  
  ! Allocate arrays
  allocate(postconv(nvar_in,nZMean,nZVar,nfiles))
  
  return
end subroutine table_init
  

! ================================================ !
! Convolute the data with pdfs of Mixture Fraction !
! ================================================ !
subroutine table_convolute(file)
  use table
  implicit none
  
  integer, intent(in) :: file
  integer :: izm, izv, k, var
  real(WP) :: meanVal
  
  ! Prepare the convolution
  allocate(pdf(nPoints))
  
  ! Convolutes
  do izv=1,nZVar
     do izm=1,nZMean
        call create_beta_pdf(ZMean(izm),ZVar(izv))
        
        do var=1,nvar_in
           meanVal  = 0.0_WP
           
           if (trim(input_name(var)).eq.'density') then
              do k=1,nPoints
                 meanVal = meanVal + pdf(k)/input_data(k,var)
              end do
              postconv(var,izm,izv,file)  = 1.0_WP / meanVal
           else
              do k=1,nPoints
                 meanVal = meanVal + pdf(k)*input_data(k,var)
              end do
              postconv(var,izm,izv,file)  = meanVal
           end if
        end do
     end do
  end do
  
  ! Add to the 3rd dimension
  if (trim(combModel).eq.'Enthalpy Flamelet') then
     Z3(file) = TOT_ENTHALPY
  end if

  ! Finish the convolution
  deallocate(pdf)
  nullify(pdf)
  
  return
end subroutine table_convolute


! ========================================================== !
! Convert the names from FlameMaster to user specified names !
! ========================================================== !
subroutine table_convert_names
  use table
  use parser
  implicit none
  
  character(len=str_medium), dimension(:), pointer :: conversion
  character(len=str_medium) :: varname
  integer :: i,n,var
  
  ! Default name as in FlameMaster
  nvar_out = nvar_in! - 1
  allocate(output_name(nvar_out))
  output_name = input_name(1:nvar_in)!-1)
  
  ! Get the number of name conversions
  call parser_getsize('Name conversion',n)
  if (mod(n,3).ne.0) stop "table_convert_names: Problem in the definition of conversion names"
  
  ! Allocate array and read
  allocate(conversion(n))
  call parser_read('Name conversion',conversion)
  
  ! Convert the names
  n = n / 3
  do i=1,n
     varname = trim(conversion((i-1)*3+3))
     loop1:do var=1,nvar_in
        if (trim(input_name(var)).eq.trim(varname)) exit loop1
     end do loop1
     if (var.eq.nvar_in+1) then
        print*, "table_convert_names: Unknown variable name : " // varname
        stop
     end if
     output_name(var) = trim(conversion((i-1)*3+1))
  end do
  
  return
end subroutine table_convert_names


! ============================================== !
! Setup the table by mapping the third direction !
! ============================================== !
subroutine table_setup
  use table
  use parser
  implicit none
  
  integer :: izm,izv,i3,var
  real(WP) :: min3, max3
  
  integer :: file, file_up, file_down
  real(WP) :: err, err_up, err_down
  real(WP) :: alpha_up,alpha_down
  
  real(WP), dimension(:), pointer :: tmp
  character(str_short) :: scale
  
  ! Convert to a chi=0 flamelet
  if (trim(combModel).eq.'Steady Flamelet') then
     file = minloc(maxval(maxval(postconv(nvar_in,:,:,:),dim=1),dim=1),dim=1)
     print*,''
     print*,'Flamelet #',file,'used as chi=0 flamelet'
     postconv(nvar_in,:,:,file) = 0.0_WP
  end if
  
  ! Allocate final table
  allocate(output_data(nZMean,nZVar,n3,nvar_out))
  allocate(mask(nZMean,nZVar,n3))
  mask=0
  
  ! Find min and max
  select case (trim(combModel))
  case ('Enthalpy Flamelet')
     max3 = maxval(Z3)
     min3 = minval(Z3)
  case default
     max3 = maxval(postconv(nvar_in,:,:,:))
     min3 = minval(postconv(nvar_in,:,:,:))
  end select
  
  ! Linear or log progression
  call parser_read('Scale for third direction',scale)
  select case(trim(scale))
  case ('lin')
     do i3=1,n3
        Z3(i3) = min3 + real(i3-1,WP)*(max3-min3)/real(n3-1,WP)
     end do
  case ('log')
     do i3=1,n3
        Z3(i3) = min3 * (max3/min3) ** (real(i3-1,WP)/real(n3-1,WP))
     end do
  case default
     stop "table_setup: Unknown Scale for third direction"
  end select
  
  ! Loop over the three mapping directions
  do i3=1,n3
     do izv=1,nZVar
        do izm=1,nZMean
           
           tmp => postconv(nvar_in,izm,izv,:)
           
           ! Find the two files right above and right below
           err_up = huge(1.0_WP)
           err_down = -huge(1.0_WP)
           
           file_up   = 0
           file_down = 0
           
           do file=1,nfiles
              err = tmp(file) - Z3(i3)
              
              if ((err.ge.0.0_WP) .and. (err.le.err_up)) then
                 file_up = file
                 err_up = err
              end if
              if ((err.le.0.0_WP) .and. (err.ge.err_down)) then
                 file_down = file
                 err_down = err
              end if
           end do
           
           ! Interpolate
           if (file_up.eq.0 .or. file_down.eq.0) then
              if (file_up.eq.0) then
                 alpha_up   = 0.0_WP
                 alpha_down = 1.0_WP
                 file_up = 1
              end if
              if (file_down.eq.0) then
                 alpha_up   = 1.0_WP
                 alpha_down = 0.0_WP
                 file_down = 1
              end if
              ! Mask it
              mask(izm,izv,i3) = 1
           else
              if (file_up.eq.file_down) then
                 alpha_up   = 1.0_WP
                 alpha_down = 0.0_WP
              else
                 alpha_up   = (Z3(i3)-tmp(file_down)) / (tmp(file_up)-tmp(file_down))
                 alpha_down = (tmp(file_up)-Z3(i3))   / (tmp(file_up)-tmp(file_down))
              end if
           end if
           
           do var=1,nvar_out
              if (trim(input_name(var)).eq.'density') then
                 output_data(izm,izv,i3,var) =  1.0_WP/( &
                      alpha_up  /postconv(var,izm,izv,file_up) + &
                      alpha_down/postconv(var,izm,izv,file_down) )
              else
                 output_data(izm,izv,i3,var) =  &
                      alpha_up  *postconv(var,izm,izv,file_up) + &
                      alpha_down*postconv(var,izm,izv,file_down)
              end if
           end do
           
        end do
     end do
  end do
  
  return
end subroutine table_setup


! ============================================================= !
! Extent the value of the chemtable outside the physical bounds !
! ============================================================= !
subroutine table_extent
  use table
  implicit none
  
  integer  :: izm,izv,i3
  integer  :: i,j,imin,jmin
  real(WP) :: d,dmin
  
  ! Loop over the three mapping directions
  do i3=1,n3
     do izv=1,nZVar
        do izm=1,nZMean
           
           ! If masked recompute value from nearest neighboor
           if (mask(izm,izv,i3).eq.1) then

              dmin = huge(WP)
              do i=1,nZMean
                 do j=1,nZVar
                    if (mask(i,j,i3).eq.0) then
                       d = sqrt((ZMean(izm)-ZMean(i))**2+(ZVar(izv)-ZVar(j))**2)
                       if (d.lt.dmin) then
                          imin = i
                          jmin = j
                          dmin = d
                       end if
                    end if
                 end do
              end do

              output_data(izm,izv,i3,:) = output_data(imin,jmin,i3,:) 

           end if
           
        end do
     end do
  end do
  
  return
end subroutine table_extent



! ===================== !
! Print some statistics !
! ===================== !
subroutine table_stats
  use table
  implicit none
  
  integer :: var
  
  print*,''
  
  ! Min and Max of the coordinates
  print*, '** Coordinates of the table **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  write(*,11) 'ZMEAN       ', minval(ZMean), maxval(ZMean)
  write(*,11) 'ZVAR        ', minval(ZVar), maxval(ZVar)
  if (trim(combModel).eq.'Enthalpy Flamelet') then
     write(*,11) 'TOT ENTHALPY', minval(Z3), maxval(Z3)
  else
     write(*,11) input_name(nvar_in), minval(Z3), maxval(Z3)
  end if
  print*,''
  
  ! Min and Max of all the mapped quantities
  print*, '** Mapped quantities **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  do var=1,nvar_out
     write(*,11) output_name(var),minval(output_data(:,:,:,var)),maxval(output_data(:,:,:,var))
  end do
  print*,''
  
10 format (A12,'  ',A12,'  ',A12)
11 format (A12,'  ',ES12.4,'  ',ES12.4)
  
  return
end subroutine table_stats


! ===================================== !
! Write the table back to a binary file !
! ===================================== !
subroutine table_write
  use table
  use parser
  implicit none
  
  integer :: ierr,var,iunit
  character(len=str_medium) :: filename
  real(WP) :: tmp_val
  character(len=str_medium) :: tmp_str
  integer :: i,j,k
  
  ! Open the data file
  call parser_read('Table filename', filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
  
!!$  ! Write VIDA header
!!$  tmp_val = 0.0_WP
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  tmp_str = 'FPVA'
!!$  call BINARY_FILE_WRITE(iunit,tmp_str,str_medium,kind(tmp_str),ierr)
!!$  tmp_str = 'SPECIES'
!!$  call BINARY_FILE_WRITE(iunit,tmp_str,str_medium,kind(tmp_str),ierr)
!!$  tmp_val = 1.0133e5_WP
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  tmp_val = 298.0_WP
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)

  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nZMean,1,kind(nZMean),ierr)
  call BINARY_FILE_WRITE(iunit,nZVar,1,kind(nZVar),ierr)
  call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  ! Write the axis coordinates
  call BINARY_FILE_WRITE(iunit,ZMean,nZMean,kind(ZMean),ierr)
  call BINARY_FILE_WRITE(iunit,ZVar,nZVar,kind(ZVar),ierr)
  call BINARY_FILE_WRITE(iunit,Z3,n3,kind(Z3),ierr)
!!$  ! Masks
!!$  call BINARY_FILE_WRITE(iunit,mask,nZMean*nZVar*n3,kind(mask),ierr)
  ! Write additional stuff
  call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
  ! Write variable names
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
  end do
  ! Write data field
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_data(:,:,:,var),nZMean*nZVar*n3,kind(output_data),ierr)
  end do
!!$  ! Backwards order for VIDA (written in C++)
!!$  do i=1,nZMean
!!$     do j=1,nZVar
!!$        do k=1,n3
!!$           do var=1,nvar_out
!!$              call BINARY_FILE_WRITE(iunit,output_data(i,j,k,var),1,kind(output_data),ierr)
!!$           end do
!!$        end do
!!$     end do
!!$  end do
  
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine table_write
