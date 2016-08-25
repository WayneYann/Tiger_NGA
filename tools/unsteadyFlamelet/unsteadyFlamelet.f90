program unsteadyFlamelet
  use precision
  use string
  use fileio
  implicit none

  ! Arrays for needed variables (velocity, dissipation rate, scalar mean, scalar variance)
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: CHI
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: ZVAR
  real(WP), dimension(:,:,:), pointer :: TEMP
  real(WP), dimension(:,:,:), pointer :: PROG
  real(WP), dimension(:,:,:), pointer :: PAH
  real(WP), dimension(:,:,:), pointer :: PAH_S
  real(WP), dimension(:,:,:), pointer :: FV
  real(WP), dimension(:,:,:), pointer :: data

  integer :: nx,ny,nz,nvar
  integer :: xper,yper,zper
  integer :: icyl
  character(len=str_short), dimension(:), pointer :: names
  integer :: iunit,ierr,var
  character(len=str_medium) :: ufdatafilename,configfilename,directory,config,filename
  real(WP) :: dt,time
  real(WP), dimension(:), pointer :: x,y,z,xm,ym,zm
  integer, dimension(:,:), pointer :: mask

  real(WP) :: Zst

  real(WP), dimension(:), pointer :: ltime
  real(WP), dimension(:), pointer :: CHIst
  real(WP), dimension(:), pointer :: Zax
  real(WP), dimension(:), pointer :: ZVARax
  real(WP), dimension(:), pointer :: TEMPax
  real(WP), dimension(:), pointer :: PROGax
  real(WP), dimension(:), pointer :: PAHax
  real(WP), dimension(:), pointer :: PAH_Sax
  real(WP), dimension(:), pointer :: FVax

  integer, parameter :: NBins = 40

  real(WP), dimension(:), pointer :: Zpdf
  real(WP), dimension(:,:), pointer :: CHIpdf

  integer :: i,j,k

  integer :: index_vel,index_rho,index_chi,index_zmix,index_zvar,index_temp,index_prog,index_pah,index_pah_s,index_fv
  logical :: vel_present,rho_present,chi_present,zmix_present,zvar_present,temp_present,prog_present,pah_present,pah_s_present,fv_present

  ! Read file name from standard input
  print*, '============================'
  print*, '| ARTS - unsteady Flamelet |'
  print*, '============================'
  print*
  print "(a15,$)", " ufdata file : "
  read "(a)", ufdatafilename
  print "(a15,$)", " config file : "
  read "(a)", configfilename
  print "(a13,$)", " directory : "
  read "(a)", directory
  print "(a12,$)", " Z stoich : "
  read "(f)", Zst
  print*

  index_vel = -1
  index_vel = -1
  index_vel = -1
  index_vel = -1
  index_vel = -1
  index_vel = -1
  index_vel = -1
  index_pah = -1
  index_pah_s = -1
  index_fv = -1

  ! ** Open the ufdata file to read **
  call BINARY_FILE_OPEN(iunit,trim(ufdatafilename),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)

  ! Read additional stuff
  call BINARY_FILE_READ(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit,time,1,kind(dt),ierr)

  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
  end do

  ! Find the index of the needed variables
  do var=1,nvar
     select case(trim(names(var)))
     case ('U')
        print*, "axial velocity found"
        index_vel = var
        vel_present = .true.
     case ('RHO')
        print*, "density found"
        index_rho = var
        rho_present = .true.
     case ('CHI')
        print*, "dissipation rate found"
        index_chi = var
        chi_present = .true.
     case ('ZMIX')
        print*, "mixture fraction found"
        index_zmix = var
        zmix_present = .true.
     case ('ZVAR')
        print*, "mixture fraction variance found"
        index_zvar = var
        zvar_present = .true.
     case ('TEMP')
        print*, "temperature found"
        index_temp = var
        temp_present = .true.
     case ('Y_PROG')
        print*, "progress variable found"
        index_prog = var
        prog_present = .true.
     case ('Y_PAH')
        print*, "PAH mass fraction found"
        index_pah = var
        pah_present = .true.
     case ('Y_PAH_S')
        print*, "PAH mass fraction (steady model) found"
        index_pah_s = var
        pah_s_present = .true.
     case ('FV')
        print*, "volume fraction found"
        index_fv = var
        fv_present = .true.
     end select
  end do

  ! Allocate the variable
  allocate(data(nx,ny,nz))
  allocate(U(nx,ny,nz))
  allocate(RHO(nx,ny,nz))
  allocate(CHI(nx,ny,nz))
  allocate(ZMIX(nx,ny,nz))
  allocate(ZVAR(nx,ny,nz))
  allocate(TEMP(nx,ny,nz))
  allocate(PROG(nx,ny,nz))
  allocate(PAH(nx,ny,nz))
  allocate(PAH_S(nx,ny,nz))
  allocate(FV(nx,ny,nz))

  ! Read the needed data
  do var=1,nvar
     call BINARY_FILE_READ(iunit,data,nx*ny*nz,kind(data),ierr)
     if (var.eq.index_vel) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 U(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_chi) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 RHO(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_chi) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 CHI(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_zmix) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ZMIX(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_zvar) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ZVAR(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_temp) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 TEMP(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_prog) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 PROG(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_pah) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 PAH(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_pah_s) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 PAH_S(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
     if (var.eq.index_fv) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 FV(i,j,k) = data(i,j,k)
              end do
           end do
        end do
     end if
  end do

  ! Close the files
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Exit if not all needed variables found
  if (.not.(vel_present .and. chi_present .and. zmix_present .and. zvar_present)) then
     print*, 'Not all needed variables were found in the data file'
  end if

  ! ** Open the config file to read **
  call BINARY_FILE_OPEN(iunit,trim(configfilename),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)

  if (icyl.ne.1) print*, 'Only for cylindrical coordinates'
  
  ! Read grid field
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

  allocate(xm(nx),ym(ny),zm(nz))

  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do

  call system("mkdir -p "//trim(directory))

  allocate(ltime(nx))
  allocate(CHIst(nx))
  allocate(Zax(nx))
  allocate(ZVARax(nx))
  allocate(TEMPax(nx))
  allocate(PROGax(nx))
  allocate(PAHax(nx))
  allocate(PAH_Sax(nx))
  allocate(FVax(nx))

  call compute_ltime(xm,ZMIX,U,nx,ny,nz,Zst,ltime)
  call compute_chist(xm,ZMIX,CHI,nx,ny,nz,Zst,CHIst)
  call compute_ax(ZMIX(:,1,:),nx,nz,Zax)
  call compute_ax(ZVAR(:,1,:),nx,nz,ZVARax)
  if (temp_present) call compute_ax(TEMP(:,1,:),nx,nz,TEMPax)
  if (prog_present) call compute_ax(PROG(:,1,:),nx,nz,PROGax)
  if (pah_present) call compute_ax(PAH(:,1,:),nx,nz,PAHax)
  if (pah_s_present) call compute_ax(PAH_S(:,1,:),nx,nz,PAH_Sax)
  if (fv_present) call compute_ax(FV(:,1,:),nx,nz,FVax)

  ! Output the axial profile file
  filename = trim(directory) // '/CA.in'
  iunit = iopen()
  open(iunit,file=filename,form="formatted",iostat=ierr)
  write(iunit,'(10000a)') 'RPM = -1'
  write(iunit,'(10000a)') 'VarsIn = 9'
  write(iunit,'(10000a)') 'Time(s) ', 'x ', 'Pressure(Pa) ', 'TOx(K) ', 'TFuel(K) ', 'Sci(1/s) ', 'ZR ', 'ZMean ', 'ZVar '
  do i=1,nx
     if (xm(i).gt.0.0_WP) then
        write(iunit,'(10000ES20.12)') ltime(i), xm(i), 1.0133E+05, 298.0, 298.0, CHIst(i), 1.0, Zax(i), ZVARax(i)
     end if
  end do
  close(iclose(iunit))

  ! Output the dissipation rate pdf as a function of the Lagrangian time
  allocate(Zpdf(NBins+1))
  call compute_Zpdf(Zpdf,NBins)
  filename = 'zi.dat'
  iunit = iopen()
  open(iunit,file=filename,form="formatted",iostat=ierr)
  do i=1,NBins+1
     write(iunit,'(10000F6.4)') Zpdf(i)
  end do
  close(iclose(iunit))

  allocate(CHIpdf(nx,NBins+1))
  call compute_chiofZ(Zpdf,CHIpdf,CHI,ZMIX,nx,ny,nz,NBins)
  filename = 'chi.dat'
  iunit = iopen()
  open(iunit,file=filename,form="formatted",iostat=ierr)
  do i=1,nx
     if (xm(i).gt.0.0_WP) then
        do j=1,NBins+1
           write(iunit,'(10000ES20.12)') CHIpdf(i,j)
        end do
     end if
  end do
  close(iclose(iunit))

  ! Output the Lagrangian time
  filename = 'Time.dat'
  iunit = iopen()
  open(iunit,file=filename,form="formatted",iostat=ierr)
  do i=1,nx
     if (xm(i).gt.0.0_WP) then
        write(iunit,'(10000ES20.12)') ltime(i)
     end if
  end do
  close(iclose(iunit))

  ! Output the temperature as a function of the Lagrangian time
  if (temp_present) then
     filename = "Temp.dat"
     iunit = iopen()
     open(iunit,file=filename,form="formatted",iostat=ierr)
     do i=1,nx
        if (xm(i).gt.0.0_WP) then
           write(iunit,'(3ES20.12)') ltime(i), xm(i), TEMPax(i)
        end if
     end do
     close(iclose(iunit))
  end if

  ! Output the progress variable as a function of the Lagrangian time
  if (prog_present) then
     filename = "Prog.dat"
     iunit = iopen()
     open(iunit,file=filename,form="formatted",iostat=ierr)
     do i=1,nx
        if (xm(i).gt.0.0_WP) then
           write(iunit,'(3ES20.12)') ltime(i), xm(i), PROGax(i)
        end if
     end do
     close(iclose(iunit))
  end if

  ! Output the aromatic mass fraction as a function of the Lagrangian time
  if (pah_present) then
     filename = "PAH.dat"
     iunit = iopen()
     open(iunit,file=filename,form="formatted",iostat=ierr)
     do i=1,nx
        if (xm(i).gt.0.0_WP) then
           write(iunit,'(3ES20.12)') ltime(i), xm(i), PAHax(i)
        end if
     end do
     close(iclose(iunit))
  end if

  ! Output the aromatic mass fraction (steady model) as a function of the Lagrangian time
  if (pah_s_present) then
     filename = "PAH_S.dat"
     iunit = iopen()
     open(iunit,file=filename,form="formatted",iostat=ierr)
     do i=1,nx
        if (xm(i).gt.0.0_WP) then
           write(iunit,'(3ES20.12)') ltime(i), xm(i), PAH_Sax(i)
        end if
     end do
     close(iclose(iunit))
  end if

  ! Output the soot volume fraction as a function of the Lagrangian time
  if (fv_present) then
     filename = "FV.dat"
     iunit = iopen()
     open(iunit,file=filename,form="formatted",iostat=ierr)
     do i=1,nx
        if (xm(i).gt.0.0_WP) then
           write(iunit,'(3ES20.12)') ltime(i), xm(i), FVax(i)
        end if
     end do
     close(iclose(iunit))
  end if

end program unsteadyFlamelet


subroutine compute_ax(SC,nx,nz,SCax)
  use precision
  implicit none

  integer, intent(in) :: nx,nz
  real(WP), dimension(nx,nz), intent(in) :: SC
  real(WP), dimension(nx), intent(out) :: SCax
  integer :: i,k

  SCax = 0.0_WP

  do i=1,nx
     do k=1,nz
        SCax(i) = SCax(i) + SC(i,k)
     end do
  end do

  SCax = SCax / nz

  return
end subroutine compute_ax


subroutine compute_ltime(xm,ZMIX,U,nx,ny,nz,Zst,ltime)
  use precision
  implicit none

  integer, intent(in) :: nx,ny,nz
  real(WP), dimension(nx,ny,nz), intent(in) :: ZMIX,U
  real(WP), dimension(nx), intent(in) :: xm
  real(WP), intent(in) :: Zst
  real(WP), dimension(nx), intent(out) :: ltime
  integer :: i,j,k

  real(WP), dimension(nx) :: Ust

  ltime = 0.0_WP
  Ust = 0.0_WP

  do i=1,nx
     do k=1,nz
        if (ZMIX(i,1,k).gt.Zst) then
           loop1: do j=1,ny
              if (ZMIX(i,j,k).lt.Zst) then
                 Ust(i) = Ust(i) + U(i,j,k) - (U(i,j,k)-U(i,j-1,k))*((ZMIX(i,j,k)-Zst)/(ZMIX(i,j,k)-ZMIX(i,j-1,k)))
                 !if (i.eq.2) print*, U(i,j,k), Zst, ZMIX(i,j,k)
                 exit loop1
              end if
           end do loop1
        else
           Ust(i) = Ust(i) + U(i,1,k)
        end if
     end do
  end do
  Ust = Ust / nz

  do i=1,nx
     if (xm(i).lt.0.0_WP) then
        ltime(i) = 0.0_WP
     elseif (xm(i-1).lt.0.0_WP) then
        ltime(i) = 0.5_WP * (1.0_WP/Ust(i))*(xm(i))
     else
        ltime(i) = ltime(i-1) + 0.5_WP*(1.0_WP/Ust(i)+1.0_WP/Ust(i-1))*(xm(i)-xm(i-1))
     end if
  end do

  return
end subroutine compute_ltime


subroutine compute_Zpdf(Zpdf,NBins)
  use precision
  implicit none
  
  integer, intent(in) :: NBins
  real(WP), dimension(NBins+1), intent(out) :: Zpdf
  integer :: i

  do i=1,NBins+1
     Zpdf(i) = (i-1)*(1.0_WP/real(NBins,WP))
  end do

  return 
end subroutine compute_Zpdf


subroutine compute_chiofZ(Zpdf,CHIpdf,CHI,Z,nx,ny,nz,NBins)
  use precision
  implicit none

  integer, intent(in) :: nx,ny,nz,NBins

  real(WP), dimension(NBins+1), intent(in) :: Zpdf
  real(WP), dimension(nx,NBins+1), intent(out) :: CHIpdf

  real(WP), dimension(nx,ny,nz), intent(in) :: CHI
  real(WP), dimension(nx,ny,nz), intent(in) :: Z
  real(WP), dimension(nx,NBins+1) :: NUMpdf

  integer :: i,j,k,l,m,n

  real(WP) :: chi_lower, chi_upper, z_lower, z_upper

  CHIpdf = 0.0_WP
  NUMpdf = 0.0_WP

  do i=1,nx
     do j=1,ny
        loop1:do k=1,nz
           loop2:do l=2,NBins+1
              if (Z(i,j,k).lt.(0.5_WP*(Zpdf(l)+Zpdf(l-1)))) then
                 CHIpdf(i,l-1) = CHIpdf(i,l-1) + CHI(i,j,k)
                 !CHIpdf(i,l-1) = CHIpdf(i,l-1) + log(CHI(i,j,k))
                 NUMpdf(i,l-1) = NUMpdf(i,l-1) + 1.0_WP
                 cycle loop1
              end if
           end do loop2
           CHIpdf(i,NBins+1) = CHIpdf(i,NBins+1) + CHI(i,j,k)
           !CHIpdf(i,NBins+1) = CHIpdf(i,NBins+1) + log(CHI(i,j,k))
           NUMpdf(i,NBins+1) = NUMpdf(i,NBins+1) + 1.0_WP
        end do loop1
     end do
  end do

  do i=1,nx
     CHIpdf(i,1) = 0.0_WP
     do l=2,NBins
        if (NUMpdf(i,l).gt.0.0_WP) CHIpdf(i,l) = CHIpdf(i,l) / NUMpdf(i,l)
        !if (NUMpdf(i,l).gt.0.0_WP) CHIpdf(i,l) = exp(CHIpdf(i,l) / NUMpdf(i,l))
     end do
     CHIpdf(i,NBins+1) = 0.0_WP
     z_lower = -1.0_WP
     z_upper = 2.0_WP
     do l=2,NBins
        if (CHIpdf(i,l).eq.0.0_WP) then
           loop3:do m=l-1,1,-1
              if (CHIpdf(i,m).ne.0.0_WP) then
                 chi_lower = CHIpdf(i,m)
                 z_lower = Zpdf(m)
                 exit loop3
              end if
           end do loop3
           loop4:do n=l+1,NBins+1
              if (CHIpdf(i,n).ne.0.0_WP) then
                 chi_upper = CHIpdf(i,n)
                 z_upper = Zpdf(n)
                 exit loop4
              end if
           end do loop4
           if (z_lower.gt.-0.9_WP .and. z_upper.lt.1.9_WP) then
              CHIpdf(i,l) = chi_lower + (chi_upper-chi_lower)*(Zpdf(l)-z_lower)/(z_upper-z_lower)
              !CHIpdf(i,l) = chi_lower*(chi_upper/chi_lower)**((Zpdf(l)-z_lower)/(z_upper-z_lower))
           end if
        end if
     end do
  end do

  return
end subroutine compute_chiofZ


subroutine compute_chist(xm,ZMIX,CHI,nx,ny,nz,Zst,CHIst)
  use precision
  implicit none

  integer, intent(in) :: nx,ny,nz
  real(WP), dimension(nx,ny,nz), intent(in) :: ZMIX,CHI
  real(WP), dimension(nx), intent(in) :: xm
  real(WP), intent(in) :: Zst
  real(WP), dimension(nx), intent(out) :: CHIst
  integer :: i,j,k

  CHIst = 0.0_WP

  do i=1,nx
     do k=1,nz
        if (ZMIX(i,1,k).gt.Zst) then
           loop1: do j=1,ny
              if (ZMIX(i,j,k).lt.Zst) then
                 CHIst(i) = CHIst(i) + CHI(i,j,k) - (CHI(i,j,k)-CHI(i,j-1,k))*((ZMIX(i,j,k)-Zst)/(ZMIX(i,j,k)-ZMIX(i,j-1,k)))
                 !CHIst(i) = CHIst(i) + log(CHI(i,j,k)) - (log(CHI(i,j,k))-log(CHI(i,j-1,k)))*((ZMIX(i,j,k)-Zst)/(ZMIX(i,j,k)-ZMIX(i,j-1,k)))
                 !if (i.eq.2) print*, CHI(i,j,k), ZMIX(i,j,k)
                 exit loop1
              end if
           end do loop1
        else
           CHIst(i) = CHIst(i) + CHI(i,1,k)
           !CHIst(i) = CHIst(i) + log(CHI(i,1,k))
        end if
     end do
  end do
  CHIst = CHIst / nz
  !CHIst = exp(CHIst / nz)

  return
end subroutine compute_chist
