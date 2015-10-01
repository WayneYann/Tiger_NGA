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
     time = 0.0025_WP ! SD
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
      
 case default
     stop "Unknown choice"
  end select
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  
end program editData
