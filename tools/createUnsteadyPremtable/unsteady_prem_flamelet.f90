module unsteady_prem_flamelet
  use precision
  use string
  implicit none
  
  ! Combustiopn model
  character(len=str_medium) :: combModel
  
  ! Number of points in the flamelet
  integer :: nPoints
  
  ! Coordinate in mixture fraction space
  real(WP), dimension(:), pointer :: C
  
  ! List of names/variables to get from FlameMaster files
  integer :: nvar_in
  character(len=str_medium), dimension(:), pointer :: input_name
  real(WP), dimension(:,:), pointer :: input_data
  logical, dimension(:), pointer :: found
  
  ! FPVA constituant variables
  integer :: nFPVA
  character(len=str_short), dimension(:), pointer :: FPVA_name

  ! Modified progress variable definition
  logical :: modprog

  ! Mixture fraction variance transport equation
  logical :: vartran

  ! Soot radiation
  logical :: sootrad
  
  ! Particular arrays for steady flamelets / FPVA
  !real(WP), dimension(:), pointer :: PROG
  real(WP), dimension(:), pointer :: SRC_PROG
  real(WP), dimension(:), pointer :: ENTHALPY
  real(WP), dimension(:), pointer :: SRC_RAD

  ! Particular arrays for variance transport equation
  real(WP), dimension(:), pointer :: Z2RHODOT
  real(WP), dimension(:), pointer :: ZSRC_ZMIX

  ! Particular variables for FPVA
  real(WP) :: PHI
  real(WP) :: Z
  real(WP) :: Zst

end module unsteady_prem_flamelet


! ============================================== !
! FLAMELET INITIALIZATION                        !
! Convert the names to be those from FlameMaster !
! ============================================== !
subroutine unsteady_prem_flamelet_init
  use unsteady_prem_flamelet
  use parser
  implicit none
  
  ! Combustion model dependent parameters
  call parser_read('Combustion model',combModel)
  
  ! Get number of additionnal variables to store in the table
  call parser_getsize("FlameMaster variables", nvar_in)

  ! Treat the case the model is FPVA
  select case(trim(combModel))
  case ('RPFPVA')
     call parser_getsize('FPVA variables',nFPVA)
     allocate(FPVA_name(nFPVA))
     call parser_read('FPVA variables',FPVA_name)

     call parser_read('Modified progress variable',modprog,.false.)

     call parser_read('Soot radiation',sootrad,.false.)
     
     call parser_read('Solve variance transport equation',vartran,.false.)

     if (vartran) then
        allocate(input_name(nvar_in+5))
     else
        allocate(input_name(nvar_in+3))
     end if
     call parser_read("FlameMaster variables", input_name(1:nvar_in))
     
     if (vartran) then
        nvar_in = nvar_in + 5
        input_name(nvar_in-4) = 'Z2RHODOT'
        input_name(nvar_in-3) = 'ZSRC_ZMIX'
     else
        nvar_in = nvar_in + 3
     end if
     input_name(nvar_in-2) = 'SRC_PROG'
     input_name(nvar_in-1) = 'SRC_RAD'
     input_name(nvar_in)   = 'ENTHALPY'
     
  case default
     stop "flamelet_init: Unknown combustion model"
  end select
  
  ! Allocate array to specify wether the variables have been found
  allocate(found(nvar_in))
  
  return
end subroutine unsteady_prem_flamelet_init



! ================================================ !
! Read a flamelet file and store its value in data !
! ================================================ !
subroutine unsteady_prem_flamelet_readfile(filename)
  use unsteady_prem_flamelet
  use fileio
  use parser
  implicit none
  
  character(len=*), intent(in) :: filename
  integer :: iunit, ierr, var, n, nlines, index1, index_zmixsrc, index_dimer
  integer :: index_rhodot, index_enthdot, index_pahsrc_pos, index_pahsrc_neg, index_sootrad
  character(len=str_long) :: buffer
  character(len=8*str_long) :: line
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
        allocate(C(nPoints))

        if (trim(combModel).eq.'RPFPVA') then
           if (vartran) then
              Z2RHODOT  => input_data(:,nvar_in-4)
              ZSRC_ZMIX => input_data(:,nvar_in-3)
           end if
           SRC_PROG => input_data(:,nvar_in-2)
           SRC_RAD  => input_data(:,nvar_in-1)
           ENTHALPY => input_data(:,nvar_in)
           SRC_PROG = 0.0_WP
           SRC_RAD  = 0.0_WP
           ENTHALPY = 0.0_WP
        end if
     end if

     ! If prem flame, read in phi
     if (trim(combModel).eq.'RPFPVA') then
        if (index(buffer,'fuel-air-equivalence-ratio').ne.0) then
           read(buffer(index(buffer,'=')+1:),*) PHI
           call parser_read('Stoichiometric mixture fraction',Zst)
           Z = PHI/((1-Zst)/Zst + PHI)
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
           input_data(:,var) = input_data(:,var) * tmp
           found(var) = .true.
        end if
     end if
     
     ! Is it a variable for FPVA ?
     if (trim(combModel).eq.'RPFPVA') then
        if (.not.modprog) then
           loop2:do var=1,nFPVA
              tmpname = 'massfraction-' // trim(FPVA_name(var))
              if (trim(varname).eq.trim(tmpname)) then
                 read(line,*) tmp
                 C = C + tmp
                 exit loop2
              end if
              tmpname = 'ProdRate' // trim(FPVA_name(var))
              if (trim(varname).eq.trim(tmpname)) then
                 read(line,*) tmp
                 SRC_PROG = SRC_PROG + tmp
                 found(nvar_in-2) = .true.
                 exit loop2
              end if
           end do loop2
        else
           if (trim(varname).eq.'ProgRat') then
              read(line,*) tmp
              C = tmp
           end if
           if (trim(varname).eq.'ProgSrc') then
              read(line,*) tmp
              SRC_PROG = tmp
              found(nvar_in-2) = .true.
           end if
        end if
     end if

     ! Is it the enthalpy?
     if (trim(combModel).eq.'RPFPVA') then
!!$        if (trim(varname).eq.'TotalEnthalpy') then
!!$           read(line,*) tmp
!!$           ENTHALPY = tmp
!!$           found(nvar_in) = .true.
!!$        end if
!!$        if (trim(varname).eq.'TotalEnthalpy2') then
!!$           read(line,*) tmp
!!$           ENTHALPY = tmp
!!$           found(nvar_in) = .true.
!!$        end if
        if (trim(varname).eq.'EnthLoss') then
           read(line,*) tmp
           ENTHALPY = tmp
           found(nvar_in) = .true.
        end if
     end if

     ! Is it the radiation source term?
     if (trim(combModel).eq.'RPFPVA') then
        if (trim(varname).eq.'GasRadiation') then
           read(line,*) tmp
           SRC_RAD = tmp
           found(nvar_in-1) = .true.
        end if
        if (trim(varname).eq.'RadiationSource') then
           read(line,*) tmp
           SRC_RAD = -tmp
           found(nvar_in-1) = .true.
        end if
     end if

     ! Do we want that variable?
     loop3:do var=1,nvar_in
        if (trim(input_name(var)).eq.varname) then
           read(line,*) input_data(:,var)
           found(var) = .true.

           if (trim(varname).eq.'ZBilgerSrc') then
              index_zmixsrc = var
           end if

           if (trim(varname).eq.'Dimer_ProdRate') then
              index_dimer = var
           end if

           if (trim(varname).eq.'RhoDot') then
              index_rhodot = var
           end if

!!$           if (trim(varname).eq.'EnthDot') then
!!$              index_enthdot = var
!!$           end if

           if (trim(varname).eq.'ProdRatePos-PAH') then
              index_pahsrc_pos = var
           end if

           if (trim(varname).eq.'ProdRateNeg-PAH') then
              index_pahsrc_neg = var
           end if

           if (trim(varname).eq.'SootRadCoeff') then
              index_sootrad = var
           end if

           exit loop3
        elseif (trim(input_name(var)).eq.'cp' .and. varname.eq.'Cp') then
           ! Specific heat has different format depending on flamelet type
           read(line,*) input_data(:,var)
           found(var) = .true.
        end if

     end do loop3

     if (vartran) then
        found(nvar_in-4) = .true.
        found(nvar_in-3) = .true.
     end if
     
  end do loop1

!!$  ! Progress variable derivative squared
!!$  do var=1,nvar_in
!!$     if (trim(input_name(var)).eq.'ProgDer') then
!!$        found(var) = .true.
!!$        input_data(1,var) = (C(2)-C(1)) / (Z(2)-Z(1))
!!$        input_data(2:nPoints-1,var) = (PROG(3:nPoints)-PROG(1:nPoints-2)) / (Z(3:nPoints)-Z(1:nPoints-2))
!!$        input_data(nPoints,var) = (PROG(nPoints)-PROG(nPoints-1)) / (Z(nPoints)-Z(nPoints-1))
!!$        input_data(:,var) = input_data(:,var)*input_data(:,var)
!!$     end if
!!$  end do

  ! Heat loss parameter source term source term
  do var=1,nvar_in
     if (trim(input_name(var)).eq.'EnthDotDummy') then
        found(var) = .true.
        input_data(:,var) = input_data(:,index_rhodot)*ENTHALPY
        index_enthdot = var
     end if
  end do

  ! Temperature gradient in progress variable
  do var=1,nvar_in
     if (trim(input_name(var)).eq.'dTdC') then
        found(var) = .true.
        input_data(:,var) = 1.0_WP
     end if
  end do

  if (.not.found(nvar_in-1))then
     found(nvar_in-1) = .true.
  endif

  do var=1,nvar_in
     if (.not.found(var)) then
        print*,"Variable ",trim(input_name(var))," not found in flamelet file"
        stop
     end if
  end do
  
  ! For FPVA, divide source term by density
  if (trim(combModel).eq.'RPFPVA') then
     loop6:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop6
     end do loop6
     if (var.le.nvar_in) SRC_PROG = SRC_PROG / input_data(:,var)
  end if

  ! For mixture fraction source term, divide by density
  if (vartran) then
     if (found(index_zmixsrc)) then
        loop7:do var=1,nvar_in
           if (trim(input_name(var)).eq.'density') exit loop7
        end do loop7
        if (var.le.nvar_in) input_data(:,index_zmixsrc) = input_data(:,index_zmixsrc) / input_data(:,var)
     end if
  end if

  ! Mixture fraction source term for mixture fraction squared equation
  if (vartran) then
     if (found(index_zmixsrc)) then
        ZSRC_ZMIX = input_data(:,index_zmixsrc) * Z
     end if
  end if

  ! For density source term, divide by density, REMULTIPLY AFTER CONVOLUTION
  if (sootrad) then
     if (found(index_rhodot)) then
        loop8:do var=1,nvar_in
           if (trim(input_name(var)).eq.'density') exit loop8
        end do loop8
        if (var.le.nvar_in) input_data(:,index_rhodot) = input_data(:,index_rhodot) / input_data(:,var)
     end if
  end if

  ! Density source term for mixture fraction squared equation
  if (vartran) then
     if (found(index_rhodot)) then
        Z2RHODOT = input_data(:,index_rhodot) * Z**2
     end if
  end if

  ! For enthalpy source term, divide by density
  if (sootrad) then
     if (found(index_enthdot)) then
        loop12:do var=1,nvar_in
           if (trim(input_name(var)).eq.'density') exit loop12
        end do loop12
        if (var.le.nvar_in) input_data(:,index_enthdot) = input_data(:,index_enthdot) / input_data(:,var)
     end if
  end if

  ! For radiation source term, divide by density
  if (trim(combModel).eq.'RPFPVA') then
     loop13:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop13
     end do loop13
     if (var.le.nvar_in) SRC_RAD = SRC_RAD / input_data(:,var)
  end if

  ! For soot radiation source term coefficient, divide by density
  if (trim(combModel).eq.'RPFPVA' .and. sootrad) then
     loop14:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop14
     end do loop14
     if (var.le.nvar_in) input_data(:,index_sootrad) = input_data(:,index_sootrad) / input_data(:,var)
  end if

  ! For dimer production rate, divide by density, REMULTIPLY AFTER CONVOLUTION
  if (sootrad) then
     if (found(index_dimer)) then
        loop10:do var=1,nvar_in
           if (trim(input_name(var)).eq.'density') exit loop10
        end do loop10
        if (var.le.nvar_in) input_data(:,index_dimer) = input_data(:,index_dimer) / input_data(:,var)
     end if
  end if

  ! For PAH source terms, divide by density, REMULTIPLY AFTER CONVOLUTION
  if (sootrad) then
     if (found(index_pahsrc_pos) .and. found(index_pahsrc_neg)) then
        loop11:do var=1,nvar_in
           if (trim(input_name(var)).eq.'density') exit loop11
        end do loop11
        if (var.le.nvar_in) input_data(:,index_pahsrc_pos) = input_data(:,index_pahsrc_pos) / input_data(:,var)
        if (var.le.nvar_in) input_data(:,index_pahsrc_neg) = input_data(:,index_pahsrc_neg) / input_data(:,var)
     elseif (found(index_pahsrc_pos) .or. found(index_pahsrc_neg)) then
        write(*,*) 'Both positive and negative source terms needed for PAH transport equation'
     end if
  end if

  ! For dimer production rate, convert from kmol to mol
  if (sootrad) then
     if (found(index_dimer)) input_data(:,index_dimer) = input_data(:,index_dimer) * 1000.0_WP
  end if

  !write(*,'(a20,7es20.12)'), 'RawData', Z(700), input_data(700,16), input_data(700,18), Z2RHODOT(700), input_data(700,nvar_in-5), ZSRC_ZMIX(700), input_data(700,nvar_in-4)
  !minval(input_data(:,22))

  ! Deallocate
  deallocate(tmp)
  nullify(tmp)
  close(iclose(iunit))

  return
end subroutine unsteady_prem_flamelet_readfile


! ========================================= !
! Deallocate the data array for new reading !
! ========================================= !
subroutine unsteady_prem_flamelet_cleanup
  use unsteady_prem_flamelet
  implicit none
  
  deallocate(input_data)
  deallocate(C)
  nullify(input_data)
  nullify(C)
  
  return
end subroutine unsteady_prem_flamelet_cleanup
