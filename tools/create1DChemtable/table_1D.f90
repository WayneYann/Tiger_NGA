module table_1D
  use flamelet_1D
  use precision
  use string
  implicit none
  
  ! Table parameters
  integer :: nZMean
  real(WP), dimension(:), pointer :: ZMean
  
  ! Delta pdf
  real(WP), dimension(:), pointer :: pdf
  
  ! Flamelet file
  character(len=str_long) :: flameletfile
  
  ! Variables to be wrtitten into the table
  integer :: nvar_out
  character(len=str_medium), dimension(:), pointer :: output_name
  real(WP), dimension(:,:), pointer :: output_data
  
  ! Variable after convolution
  real(WP), dimension(:,:), pointer :: postconv
  
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
  
  ! Precompute the PDF for a delta distribution with
  ! -> mean mixture fraction : zm
  subroutine create_delta_pdf(zm)
    implicit none 

    real(WP), intent(in) :: zm
    integer :: index1

    pdf = 0.0_WP

    ! Zero mean
    if (zm.le.1.0E-10_WP) then
       pdf(1) = 1.0_WP
       return
    end if
    
    ! Max mean
    if (zm.ge.1.0_WP-1.0E-10_WP) then
       pdf(nPoints) = 1.0_WP
       return
    end if

    ! Other mean : delta at Z=zm
    index1 = 1
    do while (Z(index1).lt.zm)
       index1 = index1+1
    end do
    pdf(index1-1) = (Z(index1)-zm)  /(Z(index1)-Z(index1-1))
    pdf(index1)   = (zm-Z(index1-1))/(Z(index1)-Z(index1-1))

    return
  end subroutine create_delta_pdf
  
end module table_1D


subroutine table_1D_init
  use table_1D
  use parser
  implicit none
  
  ! Read the dimension of the final table
  call parser_read('Number of points for mean Z', nZMean)
  allocate(ZMean(nZMean))
  
  ! Create the table dimension
  call create_zmean
  
  ! Allocate arrays
  allocate(postconv(nvar_in,nZMean))
  
  return
end subroutine table_1D_init
  

! ================================================ !
! Convolute the data with pdfs of Mixture Fraction !
! ================================================ !
subroutine table_1D_convolute
  use table_1D
  implicit none
  
  integer :: izm, izv, k, var
  real(WP) :: meanVal
  
  ! Prepare the convolution
  allocate(pdf(nPoints))
  
  ! Delta distribution for mixture fraction
  ! -> Interpolate flamelet onto table coordinate
  do izm=1,nZMean

     call create_delta_pdf(ZMean(izm))

     do var=1,nvar_in
        meanVal = 0.0_WP

        if (trim(input_name(var)).eq.'density') then
           do k=1,nPoints
              meanVal = meanVal + pdf(k)/input_data(k,var)
           end do
           postconv(var,izm)  = 1.0_WP / meanVal
        else
           do k=1,nPoints
              meanVal = meanVal + pdf(k)*input_data(k,var)
           end do
           postconv(var,izm)  = meanVal
        end if
     end do

  end do
  
  ! Finish the convolution
  deallocate(pdf)
  nullify(pdf)
  
  return
end subroutine table_1D_convolute


! ========================================================== !
! Convert the names from FlameMaster to user specified names !
! ========================================================== !
subroutine table_1D_convert_names
  use table_1D
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
  if (mod(n,3).ne.0) stop "table_1D_convert_names: Problem in the definition of conversion names"
  
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
        print*, "table_1D_convert_names: Unknown variable name : " // varname
        stop
     end if
     output_name(var) = trim(conversion((i-1)*3+1))
  end do
  
  return
end subroutine table_1D_convert_names


! ============================================== !
! Setup the table by mapping the third direction !
! ============================================== !
subroutine table_1D_setup
  use table_1D
  use parser
  implicit none
  
  integer :: izm
  
  ! Allocate final table
  allocate(output_data(nZMean,nvar_out))
  
  ! Loop over the mapping direction
  do izm=1,nZMean
     output_data(izm,:) = postconv(:,izm)
  end do
  
  return
end subroutine table_1D_setup



! ===================== !
! Print some statistics !
! ===================== !
subroutine table_1D_stats
  use table_1D
  implicit none
  
  integer :: var
  
  print*,''
  
  ! Min and Max of the coordinates
  print*, '** Coordinates of the table **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  write(*,11) 'ZMEAN       ', minval(ZMean), maxval(ZMean)
  print*,''
  
  ! Min and Max of all the mapped quantities
  print*, '** Mapped quantities **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  do var=1,nvar_out
     write(*,11) output_name(var),minval(output_data(:,var)),maxval(output_data(:,var))
  end do
  print*,''
  
10 format (A12,'  ',A12,'  ',A12)
11 format (A12,'  ',ES12.4,'  ',ES12.4)
  
  return
end subroutine table_1D_stats


! ===================================== !
! Write the table back to a binary file !
! ===================================== !
subroutine table_1D_write
  use table_1D
  use parser
  implicit none
  
  integer :: ierr,var,iunit
  character(len=str_medium) :: filename

  ! Open the data file
  call parser_read('Table filename', filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nZMean,1,kind(nZMean),ierr)
  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  ! Write the axis coordinate
  call BINARY_FILE_WRITE(iunit,ZMean,nZMean,kind(ZMean),ierr)
  ! Write additional stuff
  call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
  ! Write variable names
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
  end do
  ! Write data field
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_data(:,var),nZMean,kind(output_data),ierr)
  end do
  
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine table_1D_write
