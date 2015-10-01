module flamelet
  use precision
  use string
  implicit none
  
  ! Combustiopn model
  character(len=str_medium) :: combModel
  
  ! Number of points in the flamelet
  integer :: nPoints
  
  ! Coordinate in mixture fraction space
  real(WP), dimension(:), pointer :: Z
  
  ! List of names/variables to get from FlameMaster files
  integer :: nvar_in
  character(len=str_medium), dimension(:), pointer :: input_name
  real(WP), dimension(:,:), pointer :: input_data
  logical, dimension(:), pointer :: found
  
  ! FPVA constituant variables
  integer :: nFPVA
  character(len=str_short), dimension(:), pointer :: FPVA_name
  
  ! Particular arrays for steady flamelets / FPVA
  real(WP), dimension(:), pointer :: PROG
  real(WP), dimension(:), pointer :: SRC_PROG

  ! Particular variable for premixed flamelets
  real(WP) :: TOT_ENTHALPY
  
end module flamelet


! ============================================== !
! FLAMELET INITIALIZATION                        !
! Convert the names to be those from FlameMaster !
! ============================================== !
subroutine flamelet_init
  use flamelet
  use parser
  implicit none
  
  ! Combustion model dependent parameters
  call parser_read('Combustion model',combModel)
  
  ! Get number of additionnal variables to store in the table
  call parser_getsize("FlameMaster variables", nvar_in)
  
  ! Treat the case the model is FPVA
  select case(trim(combModel))
  case ('FPVA')
     call parser_getsize('FPVA variables',nFPVA)
     allocate(FPVA_name(nFPVA))
     call parser_read('FPVA variables',FPVA_name)
     
     allocate(input_name(nvar_in+2))
     call parser_read("FlameMaster variables", input_name(1:nvar_in))
     
     nvar_in = nvar_in + 2
     input_name(nvar_in-1) = 'SRC_PROG'
     input_name(nvar_in)   = 'PROG'
     
  case ('Steady Flamelet')
     allocate(input_name(nvar_in+1))
     call parser_read("FlameMaster variables", input_name(1:nvar_in))
     
     nvar_in = nvar_in + 1
     input_name(nvar_in)   = 'chi'

  case ('Enthalpy Flamelet')
     allocate(input_name(nvar_in+1))
     call parser_read("FlameMaster variables", input_name(1:nvar_in))
     
     nvar_in = nvar_in + 1
     input_name(nvar_in)   = 'H'
     
  case default
     stop "flamelet_init: Unknown combustion model"
  end select
  
  ! Allocate array to specify wether the variables have been found
  allocate(found(nvar_in))
  
  return
end subroutine flamelet_init



! ================================================ !
! Read a flamelet file and store its value in data !
! ================================================ !
subroutine flamelet_readfile(filename)
  use flamelet
  use fileio
  implicit none
  
  character(len=*), intent(in) :: filename
  integer :: iunit, ierr, var, n, nlines, index1
  character(len=str_long) :: buffer
  character(len=4*str_long) :: line
  character(len=str_medium) :: varname, tmpname
  real(WP), dimension(:), pointer :: tmp
  
  ! Open the file
  iunit = iopen()
  open(iunit,file=trim(filename),form='formatted',status='old',iostat=ierr)
  if (ierr.ne.0) then
     print*,"flamelet_readfile: Error opening the file : " // filename
     stop
  end if
  
  nlines = 0
  found = .false.
  ierr = 0
  buffer = ''
  do while(index(buffer,'body').eq.0)
     
     ! First get some parameters
     read(iunit,'(a)',iostat=ierr) buffer
     
     ! Get nPoints and allocate arrays
     if (index(buffer,'gridPoints').ne.0) then
        read(buffer(index(buffer,'=')+1:),*) nPoints
        nlines = ceiling(nPoints / 5.0_WP)
        
        allocate(input_data(nPoints,nvar_in))
        allocate(tmp(nPoints))
        allocate(Z(nPoints))
        
        if (trim(combModel).eq.'FPVA') then
           SRC_PROG => input_data(:,nvar_in-1)
           PROG => input_data(:,nvar_in)
           PROG = 0.0_WP
           SRC_PROG = 0.0_WP
        end if
     end if

     ! If enthalpy is used as 3rd dimension, read that quantity
     if (trim(combModel).eq.'Enthalpy Flamelet') then
        if (index(buffer,'TotEnthalpymin').ne.0) then
           read(buffer(index(buffer,'=')+1:),*) TOT_ENTHALPY
        end if
     end if
     
  end do
  
  ! Test
  if (nlines.eq.0) stop "flamelet_readfile: missing gridPoints in flamemet file"
  
  ! Preset diffusivity to 1
  loop0:do var=1,nvar_in
     if (trim(input_name(var)).eq.'diffusivity') exit loop0
  end do loop0
  if (var.le.nvar_in) input_data(:,var) = 1.0_WP
  
  ! Set Total enthalpy
  if (trim(combModel).eq.'Enthalpy Flamelet') then
     input_data(:,nvar_in) = TOT_ENTHALPY
  end if
  
  loop1:do while (ierr .eq. 0)
     
     ! Read name of variable
     read(iunit,'(a)',iostat=ierr) buffer
     if (trim(buffer).eq.'trailer') exit loop1
     
     ! Read name of variable
     read(buffer,'(a)',iostat=ierr) varname
     index1 = index(varname,' ')
     if (index1.ne.0) varname(index1:) = ''
     
     ! Read the array
     line = ''
     do n=1,nlines
        read(iunit,'(a)',iostat=ierr) buffer
        line = trim(line) // adjustl(trim(buffer))
     end do
     
     ! Is it the coordinate Z?
     if (trim(varname).eq.'Z') then
        read(line,*) Z
     end if
     
     ! Is it part of the diffusivity?
     if (trim(varname).eq.'lambda') then
        read(line,*) tmp
        loop4:do var=1,nvar_in
           if (trim(input_name(var)).eq.'diffusivity') exit loop4
        end do loop4
        if (var.le.nvar_in) then
           input_data(:,var) = input_data(:,var) * tmp
           found(var) = .true.
        end if
     end if
     if (trim(varname).eq.'cp') then
        read(line,*) tmp
        loop5:do var=1,nvar_in
           if (trim(input_name(var)).eq.'diffusivity') exit loop5
        end do loop5
        if (var.le.nvar_in) then
           input_data(:,var) = input_data(:,var) / tmp
           found(var) = .true.
        end if
     end if
     
     if (trim(varname).eq.'lambdaOverCp') then
        ! WARNING - this is a non-dimensionalized diffusivity
        read(line,*) tmp
        loop9:do var=1,nvar_in
           if (trim(input_name(var)).eq.'diffusivity') exit loop9
        end do loop9
        if (var.le.nvar_in) then
           input_data(:,var) = input_data(:,var) / tmp
           found(var) = .true.
        end if
     end if
     
     ! Is it a variable for FPVA ?
     if (trim(combModel).eq.'FPVA') then
        loop2:do var=1,nFPVA
           tmpname = 'massfraction-' // trim(FPVA_name(var))
           if (trim(varname).eq.trim(tmpname)) then
              read(line,*) tmp
              PROG = PROG + tmp
              found(nvar_in) = .true.
              exit loop2
           end if
           tmpname = 'ProdRate' // trim(FPVA_name(var))
           if (trim(varname).eq.trim(tmpname)) then
              read(line,*) tmp
              SRC_PROG = SRC_PROG + tmp
              found(nvar_in-1) = .true.
              exit loop2
           end if
        end do loop2
     end if
     
     ! Do we want that variable?
     loop3:do var=1,nvar_in
        if (trim(input_name(var)).eq.varname) then
           read(line,*) input_data(:,var)
           found(var) = .true.
           exit loop3
        end if
     end do loop3
     
  end do loop1

  if (trim(combModel).eq.'Enthalpy Flamelet') then
     do var=1,nvar_in-1
        if (.not.found(var)) then
           print*,"Variable",trim(input_name(var))," not found in flamelet file"
           stop
        end if
     end do
  else
     do var=1,nvar_in
        if (.not.found(var)) then
           print*,"Variable ",trim(input_name(var))," not found in flamelet file"
           stop
        end if
     end do
  end if
  
  ! For FPVA, divide source term by density
  if (trim(combModel).eq.'FPVA') then
     loop6:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop6
     end do loop6
     if (var.le.nvar_in) SRC_PROG = SRC_PROG / input_data(:,var)
  end if
  
  ! Force 0 at Z=1 for chi
  if (trim(combModel).eq.'Steady Flamelet') then
     input_data(nPoints,nvar_in) = 0.0_WP
  end if
  
  ! Deallocate
  deallocate(tmp)
  nullify(tmp)
  close(iclose(iunit))
  
  return
end subroutine flamelet_readfile


! ========================================= !
! Deallocate the data array for new reading !
! ========================================= !
subroutine flamelet_cleanup
  use flamelet
  implicit none
  
  deallocate(input_data)
  deallocate(Z)
  nullify(input_data)
  nullify(Z)
  
  return
end subroutine flamelet_cleanup
