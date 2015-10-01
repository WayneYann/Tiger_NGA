module param
  use precision
  use string
  implicit none

  ! --- CONFIG FILE ---

  ! Cylindrical or cartesian
  integer :: icyl
  ! Number of grid points
  integer :: nx,ny,nz
  ! Periodicity
  integer :: xper,yper,zper

  ! Grid/mesh variables for data file
  real(WP), dimension(:), pointer :: x
  real(WP), dimension(:), pointer :: y
  real(WP), dimension(:), pointer :: z
  real(WP), dimension(:), pointer :: xm
  real(WP), dimension(:), pointer :: ym
  real(WP), dimension(:), pointer :: zm
  integer, dimension(:,:), pointer :: mask

  ! --- DATA FILE ---

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nvar
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: data
  
  ! --- OPTDATA FILE ---
  
  ! Data array with all the optional variables
  integer :: nod
  character(len=str_short), dimension(:), pointer :: OD_names
  real(WP), dimension(:,:,:,:), pointer :: OD
  
  ! --- INFLOW FILE ---

  ! Inflow grid related parameters
  integer :: ntime
  real(WP), dimension(:), pointer :: t_inflow
  
  ! Number of time steps
  real(WP) :: dt_inflow,time_inflow

  ! Inflow array with all variables
  integer :: nvar_inflow
  character(len=str_short), dimension(:), pointer :: names_inflow
  real(WP), dimension(:,:,:,:), pointer :: inflow
  
  ! --- CHEMTABLE ---

  ! Coordinates of the chemtable
  integer :: n1,n2,n3
  real(WP), dimension(:), pointer :: x1,x2,x3
  ! Mask of the chemtable
  integer, dimension(:,:,:), pointer :: chem_mask

  ! Names in the chemtable
  integer :: nvar_chem
  character(len=str_medium), dimension(:), pointer :: names_chem

  ! Chemtable model
  character(len=str_medium) :: combModel

  ! Table of variables
  real(WP), dimension(:,:,:,:), pointer :: table
  

contains

  
  ! =============================== !
  ! Write the mesh file to the disk !
  ! =============================== !
  subroutine param_write_config(simulation)
    use parser
    implicit none

    character(len=str_medium) :: simulation
    integer :: ierr,iunit
    character(len=str_medium) :: filename

    ! Open the mesh file
    call parser_read('Init config file',filename)
    call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

    ! Write the mesh
    call BINARY_FILE_WRITE(iunit,simulation,str_medium,kind(simulation),ierr)
    call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
    call BINARY_FILE_WRITE(iunit,xper,1,kind(xper),ierr)
    call BINARY_FILE_WRITE(iunit,yper,1,kind(yper),ierr)
    call BINARY_FILE_WRITE(iunit,zper,1,kind(zper),ierr)
    call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
    call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit,x,nx+1,kind(x),ierr)
    call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
    call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)
    call BINARY_FILE_WRITE(iunit,mask,nx*ny,kind(mask),ierr)

    ! Close the file
    call BINARY_FILE_CLOSE(iunit,ierr)

    return
  end subroutine param_write_config


  ! =============================== !
  ! Write the data file to the disk !
  ! =============================== !
  subroutine param_write_data
    use parser
    implicit none

    integer :: ierr,iunit
    integer :: var
    character(len=str_medium) :: filename
    real(WP) :: dt,time

    ! Open the data file
    call parser_read('Init data file',filename)
    call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

    ! Write sizes
    call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
    call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
    ! Write additional stuff
    dt = 0.0_WP
    time = 0.0_WP
    call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
    call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
    ! Write variable names
    do var=1,nvar
       call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
    end do
    ! Write data field
    do var=1,nvar
       call BINARY_FILE_WRITE(iunit,data(:,:,:,var),nx*ny*nz,kind(data),ierr)
    end do
    
    call BINARY_FILE_CLOSE(iunit,ierr)

    return
  end subroutine param_write_data
  
  
  ! =============================== !
  ! Write the data file to the disk !
  ! =============================== !
  subroutine param_write_optdata
    use parser
    implicit none
    
    integer :: ierr,iunit
    integer :: var
    character(len=str_medium) :: filename
    real(WP) :: dt,time
    
    ! If optdata then continue
    if (.not.associated(OD)) return
    
    ! Open the data file
    call parser_read('Init optional data',filename)
    call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
    
    ! Write sizes
    call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
    call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit,nod,1,kind(nod),ierr)
    ! Write additional stuff
    dt = 0.0_WP
    time = 0.0_WP
    call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
    call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
    ! Write variable names
    do var=1,nod
       call BINARY_FILE_WRITE(iunit,OD_names(var),str_short,kind(OD_names),ierr)
    end do
    ! Write data field
    do var=1,nod
       call BINARY_FILE_WRITE(iunit,OD(:,:,:,var),nx*ny*nz,kind(OD),ierr)
    end do
    
    call BINARY_FILE_CLOSE(iunit,ierr)
    
    return
  end subroutine param_write_optdata
  
  
  ! ================================= !
  ! Write the inflow file to the disk !
  ! ================================= !
  subroutine param_write_inflow
     use parser
     implicit none
     
     integer :: ierr,iunit,var,n
     character(len=str_medium) :: filename
     
     ! If inflow then continue
     if (.not.associated(inflow)) return
          
     ! Read filename and frequency
     call parser_read('Init inflow file',filename)
     call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
     
     ! Write sizes
     call BINARY_FILE_WRITE(iunit,ntime,1,kind(ntime),ierr)
     call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit,nvar_inflow,1,kind(nvar_inflow),ierr)

     ! Write additional stuff
     call BINARY_FILE_WRITE(iunit,dt_inflow,1,kind(dt_inflow),ierr)
     call BINARY_FILE_WRITE(iunit,time_inflow,1,kind(time_inflow),ierr)
     ! Write variable names
     do var=1,nvar_inflow
        call BINARY_FILE_WRITE(iunit,names_inflow(var),str_short,kind(names_inflow),ierr)
     end do
     ! Write the grid
     call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
     call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
     call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)
     ! Write data field
     do n=1,ntime
        do var=1,nvar_inflow
           call BINARY_FILE_WRITE(iunit,inflow(n,:,:,var),ny*nz,kind(inflow),ierr)
        end do
     end do
    
     call BINARY_FILE_CLOSE(iunit,ierr)

     return
  end subroutine param_write_inflow 


  ! ==================================== !
  ! Write the chemtable file to the disk !
  ! ==================================== !
  subroutine param_write_chemtable
    use parser
    implicit none

      integer :: ierr,var,iunit
      character(len=str_medium) :: filename
      
      ! If chemtable then continue
      if (.not.associated(table)) return

      ! Open the data file
      call parser_read('Init chemtable file', filename)
      call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
      
      ! Write sizes
      call BINARY_FILE_WRITE(iunit,n1,1,kind(n1),ierr)
      call BINARY_FILE_WRITE(iunit,n2,1,kind(n2),ierr)
      call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
      call BINARY_FILE_WRITE(iunit,nvar_chem,1,kind(nvar_chem),ierr)
      ! Write the axis coordinates
      call BINARY_FILE_WRITE(iunit,x1,n1,kind(x1),ierr)
      call BINARY_FILE_WRITE(iunit,x2,n2,kind(x2),ierr)
      call BINARY_FILE_WRITE(iunit,x3,n3,kind(x3),ierr)
      ! Write the mask of the chemtable
      call BINARY_FILE_WRITE(iunit,chem_mask,n1*n2*n3,kind(chem_mask),ierr)
      ! Write additional stuff
      call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
      ! Write variable names
      do var=1,nvar_chem
         call BINARY_FILE_WRITE(iunit,names_chem(var),str_medium,kind(names_chem),ierr)
      end do
      ! Write data field
      do var=1,nvar_chem
         call BINARY_FILE_WRITE(iunit,table(:,:,:,var),n1*n2*n3,kind(table),ierr)
      end do
      
      call BINARY_FILE_CLOSE(iunit,ierr)
      
      return
    end subroutine param_write_chemtable
    
end module param
