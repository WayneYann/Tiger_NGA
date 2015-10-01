module flamelet_1D
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
  
end module flamelet_1D


! ============================================== !
! FLAMELET INITIALIZATION                        !
! Convert the names to be those from FlameMaster !
! ============================================== !
subroutine flamelet_1D_init
  use flamelet_1D
  use parser
  implicit none
  
  ! Combustion model dependent parameters
  call parser_read('Combustion model',combModel)
  
  ! Get number of additionnal variables to store in the table
  call parser_getsize("FlameMaster variables", nvar_in)
  
  ! Get the names of FlameMaster variables to read
  select case(trim(combModel))
  case ('Single Flamelet')
     allocate(input_name(nvar_in))
     call parser_read("FlameMaster variables",input_name)
  case default
     stop "flamelet_1D_init: Only 'Single Flamelet' for 1D tables"
  end select
  
  ! Allocate array to specify wether the variables have been found
  allocate(found(nvar_in))
  
  return
end subroutine flamelet_1D_init



! ================================================ !
! Read a flamelet file and store its value in data !
! ================================================ !
subroutine flamelet_1D_readfile(filename)
  use flamelet_1D
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
     print*,"flamelet_1D_readfile: Error opening the file : " // trim(filename), ierr
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
     end if

  end do
  
  ! Test
  if (nlines.eq.0) stop "flamelet_1D_readfile: missing gridPoints in flamemet file"
  
  ! Preset diffusivity to 1
  loop0:do var=1,nvar_in
     if (trim(input_name(var)).eq.'diffusivity') exit loop0
  end do loop0
  if (var.le.nvar_in) input_data(:,var) = 1.0_WP
  
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
     
     ! Do we want that variable?
     loop3:do var=1,nvar_in
        if (trim(input_name(var)).eq.varname) then
           read(line,*) input_data(:,var)
           found(var) = .true.
           exit loop3
        end if
     end do loop3
     
  end do loop1

  do var=1,nvar_in
     if (.not.found(var)) then
        print*,"Variable ",trim(input_name(var))," not found in flamelet file"
        stop
     end if
  end do
  
  ! Deallocate
  deallocate(tmp)
  nullify(tmp)
  close(iclose(iunit))
  
  return
end subroutine flamelet_1D_readfile


! ========================================= !
! Deallocate the data array for new reading !
! ========================================= !
subroutine flamelet_1D_cleanup
  use flamelet_1D
  implicit none
  
  deallocate(input_data)
  deallocate(Z)
  nullify(input_data)
  nullify(Z)
  
  return
end subroutine flamelet_1D_cleanup
