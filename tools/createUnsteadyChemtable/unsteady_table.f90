module unsteady_table
  use unsteady_flamelet
  use precision
  use string
  implicit none
  
  ! Table parameters
  integer :: nZVar, nZMean, n3, n4
  real(WP), dimension(:), pointer :: ZVar, ZMean, Z3, Z4
  
  ! Beta pdf
  real(WP), dimension(:), pointer :: pdf
  
  ! List of Flamelet files
  integer :: nfiles
  character(len=str_long), dimension(:), pointer :: files
  
  ! Variables to be wrtitten into the table
  integer :: nvar_out
  character(len=str_medium), dimension(:), pointer :: output_name
  real(WP), dimension(:,:,:,:,:), pointer :: output_data
  integer, dimension(:,:,:,:), pointer :: mask
  
  ! Variable after convolution
  real(WP), dimension(:,:,:,:), pointer :: postconv

  ! VIDA or NGA?
  logical :: vida_table
  
contains
  
  ! One third of the points between 0 and Zst : uniform mesh
  ! Two third between Zst and 1 : linear growth
  subroutine create_zmean
    use parser
    implicit none
    
    integer :: zm, zcut
    real(WP) :: m11,m12,m21,m22,r1,r2,delta
    real(WP) :: a,b,c
    real(WP) :: Zst
    
    call parser_read('Stoichiometric mixture fraction',Zst)
    
    zcut = nZMean/3+1
    m11 = real(nZMean**2-zcut**2,WP)
    m12 = real(nZMean-zcut,WP)
    m21 = real(2*zcut+1,WP)
    m22 = real(1,WP)
    r1 = 1.0_WP - Zst
    r2 = Zst / real(zcut-1,WP)
    delta = m11*m22-m12*m21
    a = (+ m22*r1 - m12*r2 )/delta
    b = (- m21*r1 + m11*r2 )/delta
    c = Zst - a*zcut**2-b*zcut
    
    do zm=1,zcut
       ZMean(zm) = Zst * real(zm-1,WP) / real(zcut-1,WP)
    end do
    do zm=zcut+1,nZMean
       ZMean(zm) = a*real(zm,WP)**2 + b*real(zm,WP) + c
    end do
    
    return
  end subroutine create_zmean
  
  
  ! Uniform mesh between 0 and 0.25
  subroutine create_zvar
    implicit none
    
    integer :: zv
    
    if (nZvar.eq.1) then
       Zvar(1) = 0.0_WP
       return
    end if

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
  
end module unsteady_table


subroutine unsteady_table_init
  use unsteady_table
  use parser
  implicit none
  
  ! Read the dimension of the final table
  call parser_read('Number of points for mean Z', nZMean)
  call parser_read('Number of points for variance of Z', nZVar)
  allocate(ZMean(nZMean))
  allocate(ZVar(nZVar))
  call parser_read('Number of points for third direction', n3)
  call parser_read('Number of points for fourth direction', n4)
  allocate(Z3(n3))
  allocate(Z4(n4))
  
  ! Create the first two directions of the table
  call create_zmean
  call create_zvar
  
  ! Allocate arrays
  allocate(postconv(nvar_in,nZMean,nZVar,nfiles))

  ! Table for VIDA or NGA
  call parser_read('VIDA table',vida_table,.false.)
  
  return
end subroutine unsteady_table_init
  

! ================================================ !
! Convolute the data with pdfs of Mixture Fraction !
! ================================================ !
subroutine unsteady_table_convolute(file)
  use unsteady_table
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

  ! Finish the convolution
  deallocate(pdf)
  nullify(pdf)
  
  return
end subroutine unsteady_table_convolute


! ========================================================== !
! Convert the names from FlameMaster to user specified names !
! ========================================================== !
subroutine unsteady_table_convert_names
  use unsteady_table
  use parser
  implicit none
  
  character(len=str_medium), dimension(:), pointer :: conversion
  character(len=str_medium) :: varname
  integer :: i,n,var
  
  ! Default name as in FlameMaster
  nvar_out = nvar_in
  allocate(output_name(nvar_out))
  output_name = input_name(1:nvar_in)
  
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
end subroutine unsteady_table_convert_names


! ============================================== !
! Setup the table by mapping the third direction !
! ============================================== !
subroutine unsteady_table_setup
  use unsteady_table
  use parser
  use fileio
  implicit none
  
  integer :: izm,izv,i3,i4,var,i,ierr
  integer :: var_density,var_temp
  real(WP) :: min3, max3, min4, max4
  
  integer :: file, file1, file2, file3
  real(WP) :: err, err1, err2, err3
  real(WP) :: alpha_1,alpha_2,alpha_3,denom
  
  real(WP), dimension(:), pointer :: tmp3,tmp4,Z3n,Z4n
  character(str_short) :: scale
  
  if (trim(combModel).eq.'RFPVA') then
     file1=iopen()
     open(file1,file="c_vs_h_flmlets.txt",form="formatted",iostat=ierr,status="REPLACE")
     write(file1,'(3A16)') 'Z','C','H'
     do i=1,nfiles
        do izm=nZMean/3+1,nZMean/3+1
           write(file1,'(3ES16.6)') Zmean(izm),postconv(nvar_in-2,izm,1,i),postconv(nvar_in,izm,1,i)
        end do
     end do
     close(iclose(file1))
  end if

  ! Allocate final table
  allocate(output_data(nZMean,nZVar,n3,n4,nvar_out))
  allocate(mask(nZMean,nZVar,n3,n4))
  mask=0
  
  ! Find min and max
  max3 = maxval(postconv(nvar_in-2,:,:,:))
  min3 = minval(postconv(nvar_in-2,:,:,:))
  max4 = maxval(postconv(nvar_in,:,:,:))
  min4 = minval(postconv(nvar_in,:,:,:))
  
  ! Linear or log progression
  call parser_read('Scale for third direction',scale)
  select case(trim(scale))
  case ('lin')
     do i3=1,n3
        Z3(i3) = min3 + real(i3-1,WP)*(max3-min3)/real(n3-1,WP)
     end do
  case ('quad')
     do i3=1,n3
        Z3(i3) = min3 + (max3-min3) * (real(i3-1,WP)/real(n3-1,WP))**2
     end do
  case ('log')
     do i3=1,n3
        Z3(i3) = min3 * (max3/min3) ** (real(i3-1,WP)/real(n3-1,WP))
     end do
  case default
     stop "table_setup: Unknown Scale for third direction"
  end select

  ! Linear or log progression
  call parser_read('Scale for fourth direction',scale)
  select case(trim(scale))
  case ('lin')
     do i4=1,n4
        Z4(i4) = min4 + real(i4-1,WP)*(max4-min4)/real(n4-1,WP)
     end do
  case ('quad')
     do i4=1,n4
        !Z4(i4) = min4 + (max4-min4) * (real(i4-1,WP)/real(n4-1,WP))**2
        Z4(i4) = max4 - (max4-min4) * ((real(i4-1,WP)/real(n4-1,WP))-1.0_WP)**2
     end do
  case ('log')
     do i4=1,n4
        Z4(i4) = min4 * (max4/min4) ** (real(i4-1,WP)/real(n4-1,WP))
     end do
  case default
     stop "table_setup: Unknown Scale for fourth direction"
  end select

  allocate(tmp3(nfiles))
  allocate(tmp4(nfiles))
  allocate(Z3n(n3))
  allocate(Z4n(n4))

  Z3n = (Z3-min3)/(max3-min3)
  Z4n = (Z4-min4)/(max4-min4)

  call interpolate_unsteady_flamelets_delaunay

!!$  ! Loop over the four mapping directions
!!$  do i4=1,n4
!!$     do i3=1,n3
!!$        do izv=1,nZVar
!!$           do izm=1,nZMean
!!$           
!!$              tmp3 = (postconv(nvar_in-2,izm,izv,:)-min3)/(max3-min3)
!!$              tmp4 = (postconv(nvar_in,izm,izv,:)-min4)/(max4-min4)
!!$
!!$              !if (i4.eq.1 .and. i3.eq.1) then
!!$              !   print*, "Mean/Var", ZMean(izm), ZVar(izv)
!!$              !   do file=1,nfiles
!!$              !      print*, tmp3(file), tmp4(file)
!!$              !   end do
!!$              !end if
!!$           
!!$              ! Find the three closest points
!!$              err1 = huge(1.0_WP)
!!$              err2 = huge(1.0_WP)
!!$              err3 = huge(1.0_WP)
!!$
!!$              file1 = 0
!!$              file2 = 0
!!$              file3 = 0
!!$
!!$              do file=1,nfiles
!!$                 err = sqrt((tmp3(file)-Z3n(i3))**2+(tmp4(file)-Z4n(i4))**2)
!!$                 if (err.lt.err1) then
!!$                    err3 = err2
!!$                    err2 = err1
!!$                    file3 = file2
!!$                    file2 = file1
!!$                    err1 = err
!!$                    file1 = file
!!$                 elseif (err.lt.err2) then
!!$                    err3 = err2
!!$                    file3 = file2
!!$                    err2 = err
!!$                    file2 = file
!!$                 elseif (err.lt.err3) then
!!$                    err3 = err
!!$                    file3 = file
!!$                 end if
!!$              end do
!!$
!!$              ! Interpolate
!!$              denom = tmp3(file2)*tmp4(file3)-tmp4(file2)*tmp3(file3)-tmp4(file3)*tmp3(file1)+tmp4(file2)*tmp3(file1)-tmp3(file2)*tmp4(file1)+tmp3(file3)*tmp4(file1)
!!$              alpha_1 = tmp3(file2)*tmp4(file3)-tmp4(file2)*tmp3(file3)-tmp4(file3)*Z3n(i3)+tmp4(file2)*Z3n(i3)-tmp3(file2)*z4n(i4)+tmp3(file3)*z4n(i4)
!!$              alpha_2 = Z3n(i3)*tmp4(file3)-Z4n(i4)*tmp3(file3)-tmp4(file3)*tmp3(file1)+Z4n(i4)*tmp3(file1)-Z3n(i3)*tmp4(file1)+tmp3(file3)*tmp4(file1)
!!$              alpha_3 = tmp3(file2)*Z4n(i4)-tmp4(file2)*Z3n(i3)-Z4n(i4)*tmp3(file1)+tmp4(file2)*tmp3(file1)-tmp3(file2)*tmp4(file1)+Z3n(i3)*tmp4(file1)
!!$
!!$              if (denom.eq.0) then
!!$                 alpha_1 = 1.0_WP
!!$                 alpha_2 = 0.0_WP
!!$                 alpha_3 = 0.0_WP
!!$              else
!!$                 alpha_1 = alpha_1 / denom
!!$                 alpha_2 = alpha_2 / denom
!!$                 alpha_3 = alpha_3 / denom
!!$              end if
!!$
!!$              !if (alpha_1.gt.10.0_WP .or. alpha_2.gt.10.0_WP .or. alpha_3.gt.10.0_WP) then
!!$              !   print*, alpha_1, alpha_2, alpha_3
!!$              !end if
!!$
!!$              do var=1,nvar_out
!!$                 if (trim(input_name(var)).eq.'density') then
!!$                    output_data(izm,izv,i3,i4,var) =  1.0_WP/( &
!!$                         alpha_1/postconv(var,izm,izv,file1) + &
!!$                         alpha_2/postconv(var,izm,izv,file2) + &
!!$                         alpha_3/postconv(var,izm,izv,file3) )
!!$                 else
!!$                    output_data(izm,izv,i3,i4,var) =  &
!!$                         alpha_1*postconv(var,izm,izv,file1) + &
!!$                         alpha_2*postconv(var,izm,izv,file2) + &
!!$                         alpha_3*postconv(var,izm,izv,file3)
!!$                 end if
!!$              end do
!!$              
!!$           end do
!!$        end do
!!$     end do
!!$  end do

  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'density') then
        var_density = var
     end if
  end do

  ! Remultiply by density for density and PAH source terms
  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'RhoDot') output_data(:,:,:,:,var) = output_data(:,:,:,:,var) * output_data(:,:,:,:,var_density)
     if (trim(input_name(var)).eq.'Dimer_ProdRate') output_data(:,:,:,:,var) = output_data(:,:,:,:,var) * output_data(:,:,:,:,var_density)
     if (trim(input_name(var)).eq.'ProdRatePos-PAH') output_data(:,:,:,:,var) = output_data(:,:,:,:,var) * output_data(:,:,:,:,var_density)
     if (trim(input_name(var)).eq.'ProdRateNeg-PAH') output_data(:,:,:,:,var) = output_data(:,:,:,:,var) * output_data(:,:,:,:,var_density)
  end do

  ! Compute dT/dC
  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'temperature') var_temp = var
  end do
  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'dTdC') then
        do i4 = 1,n4
           do izv = 1,nZVar
              do izm = 1,nZMean
                 output_data(izm,izv,1,i4,var) = (output_data(izm,izv,2,i4,var_temp) - output_data(izm,izv,2,i4,var_temp)) / (Z3(2)-Z3(1))
                 do i3 = 2,n3-1
                    output_data(izm,izv,i3,i4,var) = (output_data(izm,izv,i3+1,i4,var_temp) - output_data(izm,izv,i3-1,i4,var_temp)) / (Z3(i3+1)-Z3(i3-1))
                 end do
                 output_data(izm,izv,n3,i4,var) = (output_data(izm,izv,n3,i4,var_temp) - output_data(izm,izv,n3-1,i4,var_temp)) / (Z3(n3)-Z3(n3-1))
                 do i3 = 2,n3
                    if (output_data(izm,izv,i3,i4,var).eq.0.0_WP) then
                       do i = 1,i3-1
                          if (output_data(izm,izv,i,i4,var).ne.0.0_WP) output_data(izm,izv,i3,i4,var) = output_data(izm,izv,i,i4,var)
                       end do
                    end if
                 end do
              end do
           end do
        end do
     end if
  end do

  return
end subroutine unsteady_table_setup


! ============================================================= !
! Extent the value of the chemtable outside the physical bounds !
! ============================================================= !
subroutine unsteady_table_extent
  use unsteady_table
  implicit none
  
  integer  :: izm,izv,i3,i4
  integer  :: i,j,imin,jmin
  real(WP) :: d,dmin
  
  ! Loop over the three mapping directions
!!$  do i4=1,n4
!!$     do i3=1,n3
!!$        do izv=1,nZVar
!!$           do izm=1,nZMean
!!$           
!!$              ! If masked recompute value from nearest neighboor
!!$              if (mask(izm,izv,i3,i4).eq.1) then
!!$
!!$                 dmin = huge(WP)
!!$                 do i=1,nZMean
!!$                    do j=1,nZVar
!!$                       if (mask(i,j,i3).eq.0) then
!!$                          d = sqrt((ZMean(izm)-ZMean(i))**2+(ZVar(izv)-ZVar(j))**2)
!!$                          if (d.lt.dmin) then
!!$                             imin = i
!!$                             jmin = j
!!$                             dmin = d
!!$                          end if
!!$                       end if
!!$                    end do
!!$                 end do
!!$                 
!!$                 output_data(izm,izv,i3,:) = output_data(imin,jmin,i3,:) 
!!$
!!$              end if
!!$           
!!$           end do
!!$        end do
!!$     end do
!!$  end do
  
  return
end subroutine unsteady_table_extent



! ===================== !
! Print some statistics !
! ===================== !
subroutine unsteady_table_stats
  use unsteady_table
  implicit none
  
  integer :: var
  
  print*,''
  
  ! Min and Max of the coordinates
  print*, '** Coordinates of the table **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  write(*,11) 'ZMEAN       ', minval(ZMean), maxval(ZMean)
  write(*,11) 'ZVAR        ', minval(ZVar), maxval(ZVar)
  write(*,11) input_name(nvar_in-2), minval(Z3), maxval(Z3)
  write(*,11) input_name(nvar_in), minval(Z4), maxval(Z4)
  print*,''
  
  ! Min and Max of all the mapped quantities
  print*, '** Mapped quantities **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  do var=1,nvar_out
     write(*,11) output_name(var),minval(output_data(:,:,:,:,var)),maxval(output_data(:,:,:,:,var))
  end do
  print*,''
  
10 format (A12,'  ',A12,'  ',A12)
11 format (A12,'  ',ES12.4,'  ',ES12.4)
  
  return
end subroutine unsteady_table_stats


! ===================================== !
! Write the table back to a binary file !
! ===================================== !
subroutine unsteady_table_write
  use unsteady_table
  use parser
  implicit none
  
  integer :: ierr,var,iunit,i,j,k,l
  character(len=str_medium) :: filename,tmp_str
  real(WP) :: tmp_val
  
  ! Open the data file
  call parser_read('Table filename', filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

  ! Header for VIDA table
  if (vida_table) then
     tmp_val = 0.0_WP
     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
     call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
     tmp_str = 'SPECIES'
     call BINARY_FILE_WRITE(iunit,tmp_str,str_medium,kind(tmp_str),ierr)
     call parser_read('Reference pressure',tmp_val)
     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
     call parser_read('Oxidizer temperature',tmp_val)
     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
     call parser_read('Fuel temperature',tmp_val)
     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
  end if
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nZMean,1,kind(nZMean),ierr)
  call BINARY_FILE_WRITE(iunit,nZVar,1,kind(nZVar),ierr)
  call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
  call BINARY_FILE_WRITE(iunit,n4,1,kind(n4),ierr)
  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  ! Write the axis coordinates
  call BINARY_FILE_WRITE(iunit,ZMean,nZMean,kind(ZMean),ierr)
  call BINARY_FILE_WRITE(iunit,ZVar,nZVar,kind(ZVar),ierr)
  call BINARY_FILE_WRITE(iunit,Z3,n3,kind(Z3),ierr)
  call BINARY_FILE_WRITE(iunit,Z4,n4,kind(Z4),ierr)
  if (.not.vida_table) then
     ! Write additional stuff
     call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
  end if
  ! Write variable names
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
  end do
  if (.not.vida_table) then
     ! Write data field: Fortran -> column-major
     do var=1,nvar_out
        call BINARY_FILE_WRITE(iunit,output_data(:,:,:,:,var),nZMean*nZVar*n3*n4,kind(output_data),ierr)
     end do
  else
     ! Write data field: C++ -> row-major
     do i=1,nZMean
        do j=1,nZVar
           do k=1,n3
              do l=1,n4
                 do var=1,nvar_out
                    call BINARY_FILE_WRITE(iunit,output_data(i,j,k,l,var),1,kind(output_data),ierr)
                 end do
              end do
           end do
        end do
     end do
  end if
  
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Write a 3D slice -- Zero heat loss
  call BINARY_FILE_OPEN(iunit,'chemtable.3Dslice',"w",ierr)
  call BINARY_FILE_WRITE(iunit,nZMean,1,kind(nZMean),ierr)
  call BINARY_FILE_WRITE(iunit,nZVar,1,kind(nZVar),ierr)
  call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  call BINARY_FILE_WRITE(iunit,ZMean,nZMean,kind(ZMean),ierr)
  call BINARY_FILE_WRITE(iunit,ZVar,nZVar,kind(ZVar),ierr)
  call BINARY_FILE_WRITE(iunit,Z3,n3,kind(Z3),ierr)
  call BINARY_FILE_WRITE(iunit,mask(:,:,:,n4),nZMean*nZVar*n3,kind(mask),ierr)
  call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
  end do
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_data(:,:,:,n4,var),nZMean*nZVar*n3,kind(output_data),ierr)
  end do
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine unsteady_table_write


! ====================================================================================================== !
! Delaunay triangulation to interpolate the unstructured flamelet data onto a structured Cartesian table !
! ====================================================================================================== !
subroutine interpolate_unsteady_flamelets_delaunay
  use unsteady_flamelet_tools
  use unsteady_table
  use fileio
  implicit none
  
  integer :: var
  integer :: i
  integer :: izm, izv, i3, i4
  integer :: i2_2, i3_2, i4_2
  integer :: t_indx
  integer :: file1,file2,file3
  integer :: ierr
  integer :: npts_tri
  integer :: pts
  integer :: ierror
  integer, dimension(nfiles) :: file_list, file_ptr
  real(WP), dimension(2) :: a,b,c,p 
  real(WP), dimension(3) :: alpha
  real(WP) :: min_3,max_3
  real(WP) :: min_4,max_4
  real(WP) :: alpha_up, alpha_dn
  real(WP) :: Z3_file1,Z4_file1
  real(WP) :: Z3_file2,Z4_file2
  real(WP), dimension(2,nfiles) :: tri_vertices
  integer, dimension(3,5*nfiles) :: tri_indices
  integer, dimension(3,5*nfiles) :: tri_nghbrs
  integer, dimension(:), pointer :: stack
  logical :: flamelet_found
  logical :: point_inside
  logical, external :: chk_pt_in_triangle
  character(len=str_medium) :: buffer

  integer :: weights_counter
  real(WP), dimension(1200) :: weights
  real(WP) :: weights_sum, local_weight
  real(WP) :: dist_Z2, dist_Z3, dist_Z4, dist
  real(WP) :: power
  real(WP) :: Z3_target, Z4_target

  real(WP) :: err, err1
  real(WP) :: file
  real(WP), dimension(nfiles) :: tmp3,tmp4

  write(*,*) ' '
  write(*,*) 'Interpolating flamelets onto chemtable mesh ...'
  write(*,*) ' '
  
  do izv=1,nZVar
     write(*,'(A10,F6.1,A33)') 'Done with ', &
        100.0_WP*real(izv-1,WP)/real(nZVar,WP), &
        '% of triangular interpolation ...'
     do izm=1,nZMean
        
        min_3 = minval(postconv(nvar_in-2,izm,izv,:))
        min_4 = minval(postconv(nvar_in  ,izm,izv,:))
        max_3 = maxval(postconv(nvar_in-2,izm,izv,:))
        max_4 = maxval(postconv(nvar_in  ,izm,izv,:))
        max_3 = max(max_3,min_3+0.01_WP)
        max_4 = max(max_4,min_4+100.0_WP)

        npts_tri = 1
        file1 = 1
        file_list(npts_tri) = 1
        file_ptr(npts_tri) = file1
        Z3_file1 = (postconv(nvar_in-2,izm,izv,file1)-min_3)/(max_3-min_3)
        Z4_file1 = (postconv(nvar_in  ,izm,izv,file1)-min_4)/(max_4-min_4)
        tri_vertices(1,npts_tri) = Z3_file1
        tri_vertices(2,npts_tri) = Z4_file1

        do while (file1.lt.nfiles)
           file1 = file1 + 1
           Z3_file1 = (postconv(nvar_in-2,izm,izv,file1)-min_3)/(max_3-min_3)
           Z4_file1 = (postconv(nvar_in ,izm,izv,file1)-min_4)/(max_4-min_4)
           
           loop_files: do pts=1,npts_tri
              !file2 = file_list(pts)
              file2 = file_ptr(pts)
              Z3_file2 = (postconv(nvar_in-2,izm,izv,file2)-min_3)/(max_3-min_3)
              Z4_file2 = (postconv(nvar_in  ,izm,izv,file2)-min_4)/(max_4-min_4)
              
              if (abs(Z3_file1-Z3_file2).le.1.0E-3_WP .and. &
                  abs(Z4_file1-Z4_file2).le.1.0E-3_WP) then
                 exit loop_files
              end if 

              if (pts.eq.npts_tri) then
                 npts_tri = npts_tri + 1
                 file_list(npts_tri) = npts_tri 
                 file_ptr(npts_tri) = file1
                 tri_vertices(1,npts_tri) = Z3_file1
                 tri_vertices(2,npts_tri) = Z4_file1
              end if
           end do loop_files 
        end do

        pts = max(npts_tri,10)
        if (associated(stack)) deallocate(stack); nullify(stack)
        allocate(stack(1:pts))
        if (npts_tri.ge.3) then
           call dtris2( npts_tri, pts, tri_vertices(1:2,1:npts_tri), & 
                        file_list(1:npts_tri), n_triangles, & 
                        tri_indices(1:3,1:n_triangles), & 
                        tri_nghbrs(1:3,1:n_triangles), stack(1:pts), ierror )
        else
           n_triangles = 0
        end if

        if (izv.eq.1 .and. izm.eq.nZMean/3+1) then
           file1=iopen()
           open(file1,file="c_vs_h_tris.txt",form="formatted",iostat=ierr,status="REPLACE")
           write(file1,'(5A16)') '#              C','H','C_FILE','H_FILE','TEMP'
           do i=1,n_triangles
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(1,i)), &
                                       tri_vertices(2,tri_indices(1,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(1,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(1,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(1,i)))
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(2,i)), &
                                       tri_vertices(2,tri_indices(2,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(2,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(2,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(2,i)))
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(3,i)), &
                                       tri_vertices(2,tri_indices(3,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(3,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(3,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(3,i)))
           end do
           close(iclose(file1))

           do i=1,min(9999,n_triangles)
              file1=iopen()
              buffer = "triangles/c_vs_h_tri_XXXX.txt"
              write(buffer(22:25),'(i4.4)') i  
              open(file1,file=trim(buffer),form="formatted",iostat=ierr,status="REPLACE")
              write(file1,'(5A16)') '#              C','H','C_FILE','H_FILE','TEMP'
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(1,i)), &
                                       tri_vertices(2,tri_indices(1,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(1,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(1,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(1,i)))
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(2,i)), &
                                       tri_vertices(2,tri_indices(2,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(2,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(2,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(2,i)))
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(3,i)), &
                                       tri_vertices(2,tri_indices(3,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(3,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(3,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(3,i)))
              write(file1,'(5ES16.6)') tri_vertices(1,tri_indices(1,i)), &
                                       tri_vertices(2,tri_indices(1,i)), &
                                       postconv(nvar_in-2,izm,izv,file_ptr(tri_indices(1,i))), & 
                                       postconv(nvar_in  ,izm,izv,file_ptr(tri_indices(1,i))), &
                                       postconv(4        ,izm,izv,file_ptr(tri_indices(1,i)))
              close(iclose(file1))
           end do
        end if

        do i3=1,n3
           do i4=1,n4
              p(1) = (Z3(i3)-min_3) / (max_3-min_3)
              p(2) = (Z4(i4)-min_4) / (max_4-min_4)
              t_indx=0
              ! Find the triangle in which this point exists
              seektri_loop: do i=1,n_triangles
                 a(1)=tri_vertices(1,tri_indices(1,i)); a(2)=tri_vertices(2,tri_indices(1,i))
                 b(1)=tri_vertices(1,tri_indices(2,i)); b(2)=tri_vertices(2,tri_indices(2,i))
                 c(1)=tri_vertices(1,tri_indices(3,i)); c(2)=tri_vertices(2,tri_indices(3,i))
                 point_inside = chk_pt_in_triangle(a,b,c,p,alpha,1)
                 if (point_inside) then
                    t_indx = i
                    !write(*,*) 't_indx = ',t_indx
                    if (     alpha(1).gt. 1.01_WP .or. alpha(2).gt. 1.01_WP .or. alpha(3).gt. 1.01_WP & 
                        .or. alpha(1).lt.-0.01_WP .or. alpha(2).lt.-0.01_WP .or. alpha(3).lt.-0.01_WP) then
                       write(*,'(A11,3F9.4)') 'Error: alpha = ',alpha(1),alpha(2),alpha(3)
                       write(*,*) '      a = ',a
                       write(*,*) '      b = ',b
                       write(*,*) '      c = ',c
                       write(*,*) '      p = ',p
                    end if
                    exit seektri_loop
                 end if
              end do seektri_loop

              ! Interpolate based on the surrounding points
              if (t_indx.ne.0) then
                 do var=1,nvar_out
                    if (trim(input_name(var)).eq.'density') then
                       output_data(izm,izv,i3,i4,var) = 1.0_WP / ( &
                          alpha(1)/postconv(var,izm,izv,file_ptr(tri_indices(1,t_indx))) + &
                          alpha(2)/postconv(var,izm,izv,file_ptr(tri_indices(2,t_indx))) + &
                          alpha(3)/postconv(var,izm,izv,file_ptr(tri_indices(3,t_indx))) )
                    else
                       output_data(izm,izv,i3,i4,var) = &
                          alpha(1)*postconv(var,izm,izv,file_ptr(tri_indices(1,t_indx))) + &
                          alpha(2)*postconv(var,izm,izv,file_ptr(tri_indices(2,t_indx))) + &
                          alpha(3)*postconv(var,izm,izv,file_ptr(tri_indices(3,t_indx))) 
                    end if
                 end do
              else
                 file1 = 0
                 err1 = huge(1.0_WP)

                 tmp3 = (postconv(nvar_in-2,izm,izv,:)-min_3) / (max_3-min_3)
                 tmp4 = (postconv(nvar_in  ,izm,izv,:)-min_4) / (max_4-min_4)

                 do file=1,nfiles
                    err = sqrt((tmp3(file)-p(1))**2+(tmp4(file)-p(2))**2)
                    if (err.lt.err1) then
                       err1 = err
                       file1 = file
                    end if
                 end do
                 output_data(izm,izv,i3,i4,:) = postconv(:,izm,izv,file1)
                 !output_data(izm,izv,i3,i4,:) = 0.0_WP
                 !mask(izm,izv,i3,i4) = 1
              end if

           end do
        end do

     end do
  end do
        
!!$  write(*,*) '   ' 
!!$  write(*,*) 'Doing extrapolation of data ...   ' 
!!$  write(*,*) '   ' 
!!$
!!$  ! Use inverse-distance-weighted interpolation
!!$  min_3 = minval(postconv(nvar_in  ,:,:,:))
!!$  max_3 = maxval(postconv(nvar_in  ,:,:,:))
!!$  max_3 = max(max_3,min_3+0.01_WP)
!!$  min_4 = minval(postconv(nvar_in-2,:,:,:))
!!$  max_4 = maxval(postconv(nvar_in-2,:,:,:))
!!$  max_4 = max(max_4,min_4+100.0_WP)
!!$  power = 4.0_WP
!!$  do izv=1,nZVar
!!$     write(*,'(A10,F6.1,A22)') 'Done with ', &
!!$        100.0_WP*real(izv-1,WP)/real(nZVar,WP),'% of extrapolation ...'
!!$     do izm=1,nZMean
!!$        do i4=1,n4
!!$           do i3=1,n3
!!$          
!!$              if (mask_4D(izm,izv,i3,i4).eq.1) then
!!$                 pts = 0
!!$                 weights(:) = 0.0_WP
!!$                 weights_counter = 0
!!$              
!!$                 Z3_target = (Z3(i3)-min_3)/(max_3-min_3)
!!$                 Z4_target = (Z4(i4)-min_4)/(max_4-min_4)
!!$                 
!!$                 ! Set up the interpolation weights
!!$                 do i4_2=max(1,i4-4),min(i4+4,n4)
!!$                    do i3_2=max(1,i3-4),min(i3+4,n3)
!!$                       do i2_2=max(1,izm-3),min(izm+3,nZMean)
!!$                          if (mask_4D(i2_2,izv,i3_2,i4_2).eq.0) then
!!$                             Z3_file1 = (output_data_4D(i2_2,izv,i3_2,i4_2,nvar_in  )-min_3)/(max_3-min_3)
!!$                             Z4_file1 = (output_data_4D(i2_2,izv,i3_2,i4_2,nvar_in-2)-min_4)/(max_4-min_4)
!!$                      
!!$                             dist_Z2 = 250.0_WP*(ZMean(i2_2) - ZMean(izm))
!!$                             dist_Z3 = 600.0_WP*(Z3_file1 - Z3_target)
!!$                             dist_Z4 = 250.0_WP*(Z4_file1 - Z4_target)
!!$                             dist = sqrt(dist_Z2**2 + dist_Z3**2 + dist_Z4**2)
!!$
!!$                             weights_counter = weights_counter + 1
!!$                             weights(weights_counter) = 1.0_WP / (1.0E-9_WP + dist**power)
!!$                             pts = pts+1
!!$                          end if
!!$                       end do
!!$                    end do
!!$                 end do
!!$                 weights_sum = sum(weights(:))
!!$              
!!$                 if (pts.eq.0) then
!!$                    do var=1,nvar_out
!!$                       output_data_4D(izm,izv,i3,i4,var) = postconv(var,izm,izv,2)
!!$                       output_data_4D(izm,izv,i3,i4,var) = 0.0_WP
!!$                    end do
!!$                 else
!!$                    do var=1,nvar_out
!!$                       output_data_4D(izm,izv,i3,i4,var) = 0.0_WP
!!$                       weights_counter = 0
!!$                       do i4_2=max(1,i4-4),min(i4+4,n4)
!!$                          do i3_2=max(1,i3-4),min(i3+4,n3)
!!$                             do i2_2=max(1,izm-3),min(izm+3,nZMean)
!!$
!!$                                if (mask_4D(i2_2,izv,i3_2,i4_2).eq.0) then
!!$                                   weights_counter = weights_counter + 1
!!$                                   local_weight = weights(weights_counter) / weights_sum
!!$   
!!$                                   output_data_4D(izm,izv,i3,i4,var) = output_data_4D(izm,izv,i3,i4,var) & 
!!$                                                + output_data_4D(i2_2,izv,i3_2,i4_2,var) * local_weight
!!$                                end if
!!$
!!$                             end do
!!$                          end do
!!$                       end do
!!$
!!$                    end do
!!$
!!$                 end if
!!$
!!$              end if
!!$
!!$           end do
!!$        end do
!!$     end do
!!$  end do

  write(*,*) 'done with interpolate_unsteady_flamelets_delaunay'

  return
end subroutine interpolate_unsteady_flamelets_delaunay
