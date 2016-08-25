program editData
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar,nvar2,n_del,iloc,i
  character(len=str_short), dimension(:), pointer :: names,names2
  character(len=str_short), dimension(10) :: names_del
  real(WP), dimension(:,:,:), pointer :: data
  integer :: iunit1,iunit2,ierr,var,choice,iunit3,var2
  character(len=str_medium) :: filename1,filename2,filename3
  character(len=str_short) :: varname,varname2
  real(WP) :: dt,time,value
  
  ! for new mech
  character(len=str_short), dimension(:), pointer :: newmech
  integer :: nnew,nold
  real(WP), dimension(:,:,:,:), pointer :: data_tmp

  ! Read file name from standard input
  print*,'======================'
  print*,'| ARTS - data Editor |'
  print*,'======================'
  print*
  print "(a28,$)", " data file before edition : "
  read "(a)", filename1
  print "(a27,$)", " data file after edition : "
  read "(a)", filename2

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
    
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  print*,'There are ',nvar,' variables.'

  ! Allocate arrays
  allocate(data(nx,ny,nz))

  ! ** Ask what to do **
  print*
  print*, "1. Print Min/Max of variable"
  print*, "2. Add variable"
  print*, "3. Delete variable"
  print*, "4. Reset time"
  print*, "5. Empty variable list"
  print*, "6. Append data"  
  print*, "7. Rename variable" 
  print*, "8. Chop and recover domain" 
  print*, "9. Switch mechanism" 
  print "(a9,$)", "Choice : "
  read "(i1)", choice
  
  ! Case dependent operation
  select case(choice)
     
  case(1) ! Print min/max of all variables
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        print*,"min: ",minval(data)," - max: ",maxval(data)
        print*,maxloc(data)
     end do
     
  case (2) ! Add variable
     print "(a16,$)", "Variable name : "
     read "(a)", varname
     print "(a16,$)", "Default value : "
     read(*,*) value
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar+1,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     call BINARY_FILE_WRITE(iunit2,varname,str_short,kind(varname),ierr)
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     data = value
     call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     call BINARY_FILE_CLOSE(iunit2,ierr)
     
  case (3) ! Delete variables 
     print*,'You can delete at most 10 variables, press q to exit'
     n_del = 0
     do while (n_del .le. 10)
        print "(a16,$)", "Variable name : "
        read "(a)", varname
        if (varname .ne. 'q') then
           n_del = n_del + 1
           names_del(n_del) = varname
        else
           exit   
        end if
     end do

     print*,'How many are deleted? ',n_del
     print*,'They are :',names_del(1:n_del)

     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar-n_del,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)

     do var=1,nvar
        if (.not. any(names_del.eq.(trim(adjustl(names(var)))))) &
             call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        if (.not. any(names_del.eq.(trim(adjustl(names(var)))))) &
             call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
          
  case (4) ! Reset time to zero
     time = 0.0_WP ! SD
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     do var=1,nvar
        call BINARY_FILE_READ (iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
     
  case (5) ! Empty variable list
     nvar = 0
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     call BINARY_FILE_CLOSE(iunit2,ierr)

  case (6) ! Append data 
     print "(a19,$)", " extra data file : "
     read "(a)", filename3
     ! ** Open the data file to read **
     call BINARY_FILE_OPEN(iunit3,trim(filename3),"r",ierr)

     ! Read sizes
     call BINARY_FILE_READ(iunit3,nx,1,kind(nx),ierr)
     call BINARY_FILE_READ(iunit3,ny,1,kind(ny),ierr)
     call BINARY_FILE_READ(iunit3,nz,1,kind(nz),ierr)
     call BINARY_FILE_READ(iunit3,nvar2,1,kind(nvar2),ierr)
     print*,'Grid :',nx,'x',ny,'x',nz

     ! Read additional stuff
     call BINARY_FILE_READ(iunit3,dt,1,kind(dt),ierr)
     call BINARY_FILE_READ(iunit3,time,1,kind(time),ierr)
     print*,'Data file at time :',time

     ! Read variable names
     allocate(names2(nvar2))
     do var=1,nvar2
        call BINARY_FILE_READ(iunit3,names2(var),str_short,kind(names2),ierr)
     end do
     print*,'Variables : ',names2
     print*,'There are ',nvar2,' variables.'
     print*,'Appending data...'

     ! Allocate arrays
     allocate(data(nx,ny,nz))

     ! Copy from the first file
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar+nvar2,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     ! Append the second file at the end of the first file
     do var=1,nvar2
        call BINARY_FILE_WRITE(iunit2,names2(var),str_short,kind(names2),ierr)
     end do
     ! Again from the first file
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     ! Again the second file
     do var=1,nvar2
        call BINARY_FILE_READ(iunit3,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
     call BINARY_FILE_CLOSE(iunit3,ierr)
     
  case (7)  ! Rename a variable 
     print "(a16,$)", "Variable name : "
     read "(a)", varname
     print "(a20,$)", "New variable name : "
     read "(a)", varname2
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        if (trim(adjustl(names(var))).ne.trim(adjustl(varname))) then
           call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
        else
           call BINARY_FILE_WRITE(iunit2,varname2,str_short,kind(names),ierr)
        end if
     end do
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)

  case (8) ! Chop and recover domain
     print "(a22,$)", "Location index in x : "
     read *, iloc
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        do i=iloc+1,nx
           data(i,:,:) = data(iloc,:,:)
        end do
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do

  case (9) ! Switch mechanism 
     print*,'You need to hard code the variable names.'
     nnew = 158! number of species in the new mech
     nold = 46! number of species in the old mech
     print*,'Number of species in the new mech: ',nnew
     print*,'Number of species in teh old mech: ',nold
     allocate(newmech(nnew))
     allocate(data_tmp(nx,ny,nz,nold))

     ! Species names in new mech
     newmech(1)='N2'
     newmech(2)='H'
     newmech(3)='O2'
     newmech(4)='O'
     newmech(5)='OH'
     newmech(6)='H2'
     newmech(7)='H2O'
     newmech(8)='CO2'
     newmech(9)='HO2'
     newmech(10)='H2O2'
     newmech(11)='CO'
     newmech(12)='HCO'
     newmech(13)='C'
     newmech(14)='CH'
     newmech(15)='TXCH2'
     newmech(16)='CH3'
     newmech(17)='CH2O'
     newmech(18)='HCCO'
     newmech(19)='C2H'
     newmech(20)='CH2CO'
     newmech(21)='C2H2'
     newmech(22)='SXCH2'
     newmech(23)='AR'
     newmech(24)='CH3OH'
     newmech(25)='CH2OH'
     newmech(26)='CH3O'
     newmech(27)='CH4'
     newmech(28)='CH3O2'
     newmech(29)='C2H3'
     newmech(30)='C2H4'
     newmech(31)='C2H5'
     newmech(32)='HCCOH'
     newmech(33)='CH2CHO'
     newmech(34)='CH3CHO'
     newmech(35)='H2C2'
     newmech(36)='C2H5O'
     newmech(37)='NXC3H7'
     newmech(38)='C2H6'
     newmech(39)='C3H8'
     newmech(40)='C3H6'
     newmech(41)='C3H3'
     newmech(42)='PXC3H4'
     newmech(43)='AXC3H4'
     newmech(44)='SXC3H5'
     newmech(45)='NXC4H3'
     newmech(46)='C2H3CHO'
     newmech(47)='AXC3H5'
     newmech(48)='C2O'
     newmech(49)='C4H4'
     newmech(50)='C3H2'
     newmech(51)='C3H2O'
     newmech(52)='C4H2'
     newmech(53)='IXC4H3'
     newmech(54)='TXC3H5'
     newmech(55)='C3H5O'
     newmech(56)='C4H'
     newmech(57)='C8H2'
     newmech(58)='C6H2'
     newmech(59)='C4H6'
     newmech(60)='NXC4H5'
     newmech(61)='IXC4H5'
     newmech(62)='A1XC6H6'
     newmech(63)='NXC7H16'
     newmech(64)='C5H11'
     newmech(65)='PXC4H9'
     newmech(66)='C7H15'
     newmech(67)='PXC4H8'
     newmech(68)='C5H10'
     newmech(69)='C7H14'
     newmech(70)='C7H15O'
     newmech(71)='C3H7CHO'
     newmech(72)='C4H7'
     newmech(73)='C7H13'
     newmech(74)='C5H9'
     newmech(75)='C4H7O'
     newmech(76)='NXC3H7O'
     newmech(77)='IXC8H18'
     newmech(78)='YXC7H15'
     newmech(79)='IXC4H8'
     newmech(80)='IXC3H7'
     newmech(81)='TXC4H9'
     newmech(82)='CXC8H17'
     newmech(83)='YXC7H14'
     newmech(84)='DXC8H17O'
     newmech(85)='CH3COCH3'
     newmech(86)='IXC4H7'
     newmech(87)='XXC7H13'
     newmech(88)='IXC3H5CH'
     newmech(89)='TXC4H9O'
     newmech(90)='IXC4H7O'
     newmech(91)='C5H4CH2'
     newmech(92)='A1XXC6H5'
     newmech(93)='A1C2H2XC'
     newmech(94)='A1C2H3XC'
     newmech(95)='A1C2HXC8'
     newmech(96)='A1C2HYXC'
     newmech(97)='A1C2H3YX'
     newmech(98)='A2XXC10H'
     newmech(99)='A2XC10H8'
     newmech(100)='A2YXC10H'
     newmech(101)='A2C2H2AX'
     newmech(102)='A2C2H2BX'
     newmech(103)='A2C2HAXC'
     newmech(104)='A2C2HBXC'
     newmech(105)='A2C2HAYX'
     newmech(106)='A2C2HBYX'
     newmech(107)='A2R5XC12'
     newmech(108)='A2R5XXC1'
     newmech(109)='A2R5C2H2'
     newmech(110)='A2R5C2HX'
     newmech(111)='A2R5C2HY'
     newmech(112)='P2XC12H1'
     newmech(113)='P2XXC12H'
     newmech(114)='A3XXC14H'
     newmech(115)='A3XC14H1'
     newmech(116)='A3YXC14H'
     newmech(117)='A3R5XXC1'
     newmech(118)='A3R5XC16'
     newmech(119)='A4XC16H1'
     newmech(120)='A4XXC16H'
     newmech(121)='A4R5XC18'
     newmech(122)='FLTNXC16'
     newmech(123)='C5H6'
     newmech(124)='C5H5'
     newmech(125)='TXC5H5O'
     newmech(126)='C5H4O'
     newmech(127)='SXC5H5O'
     newmech(128)='C9H8'
     newmech(129)='C9H7'
     newmech(130)='A1CH2XC7'
     newmech(131)='C9H6O'
     newmech(132)='OXC6H4'
     newmech(133)='A1CH3XC7'
     newmech(134)='A1OHXC6H'
     newmech(135)='HOA1CH3X'
     newmech(136)='OA1CH3XC'
     newmech(137)='A1CH2OXC'
     newmech(138)='A1CH2OHX'
     newmech(139)='A1CHOXC7'
     newmech(140)='A1OXC6H5'
     newmech(141)='A1CH3YXC'
     newmech(142)='A1C2H4XC'
     newmech(143)='A1C2H5XC'
     newmech(144)='C8H9O2'
     newmech(145)='C8H8OOH'
     newmech(146)='OC8H7OOH'
     newmech(147)='A1CH3CH3'
     newmech(148)='A1CH3CH2'
     newmech(149)='A1CH3CHO'
     newmech(150)='A2CH3XC1'
     newmech(151)='A1CHOCH2'
     newmech(152)='A1CHOCHO'
     newmech(153)='A2OHXC10'
     newmech(154)='A2CH2XC1'
     newmech(155)='A2CH2OXC'
     newmech(156)='A2CHOXC1'
     newmech(157)='A2OXC10H'
     newmech(158)='OC6H4O'

     ! Go pass the header   
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar-nold+nnew,1,kind(nvar),ierr) ! modified
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     ! U V W P RHO dRHO   
     do var=1,6
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     ! Deal with species
     do var=1,nnew
        call BINARY_FILE_WRITE(iunit2,newmech(var),str_short,kind(names),ierr)
     end do
     do var=6+nold+1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     ! Deal with data for U V W P RHO dRHO 
     do var=1,6
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     ! Store data with the old mech
     do var=1,nold
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        data_tmp(:,:,:,var) = data
     end do
     ! Write data with the new mech
     do var=1,nnew
        loop1:do var2=7,6+nold
           if (trim(adjustl(names(var2))).eq.trim(adjustl(newmech(var)))) then
              print*,'Shared species:',names(var2)
              data = data_tmp(:,:,:,var2-6)
              exit loop1
           else
              data = 0.0_WP
           end if       
        end do loop1
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     ! Deal with data not with finitechem
     do var=6+nold+1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
      
 case default
     stop "Unknown choice"
  end select
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  
end program editData
