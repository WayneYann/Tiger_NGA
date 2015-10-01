module partition
  implicit none
  
  ! Cartesian  : x-direction
  ! Cylindrical: x-direction
  integer :: imino_
  integer :: imin_
  integer :: imaxo_
  integer :: imax_
  integer :: nx_
  integer :: nxo_
  
  ! Cartesian  : y-direction
  ! Cylindrical: r-direction
  integer :: jmino_
  integer :: jmin_
  integer :: jmaxo_
  integer :: jmax_
  integer :: ny_
  integer :: nyo_
  
  ! Cartesian  : z-direction
  ! Cylindrical: theta-direction
  integer :: kmino_
  integer :: kmin_
  integer :: kmaxo_
  integer :: kmax_
  integer :: nz_
  integer :: nzo_
  
end module partition

subroutine partition_init()
  use partition
  use geometry
  use parallel
  implicit none
  integer :: q,r

  ! Set initial partitions along x
  if (npx.gt.nx) call die('partition_init: nx has to be greater or equal to npx')
  q = nx/npx
  r = mod(nx,npx)
  if (iproc<=r) then
     nx_   = q+1
     imin_ = imin + (iproc-1)*(q+1)
  else
     nx_   = q
     imin_ = imin + r*(q+1) + (iproc-r-1)*q
  end if
  nxo_   = nx_ + 2*nover
  imax_  = imin_ + nx_ - 1
  imino_ = imin_ - nover
  imaxo_ = imax_ + nover

  ! Set initial partitions along y
  if (npy.gt.ny) call die('partition_init: ny has to be greater or equal to npy')
  q = ny/npy
  r = mod(ny,npy)
  if (jproc<=r) then
     ny_   = q+1
     jmin_ = jmin + (jproc-1)*(q+1)
  else
     ny_   = q
     jmin_ = jmin + r*(q+1) + (jproc-r-1)*q
  end if
  nyo_   = ny_ + 2*nover
  jmax_  = jmin_ + ny_ - 1
  jmino_ = jmin_ - nover
  jmaxo_ = jmax_ + nover
  
  ! Set initial partitions along z
  if (npz.gt.nz) call die('partition_init: nz has to be greater or equal to npz')
  if (icyl.eq.1 .and. npz.ne.1) call die('partition_init: npz has to be equal to 1 for cylindrical configuration')
  q = nz/npz
  r = mod(nz,npz)
  if (kproc<=r) then
     nz_   = q+1
     kmin_ = kmin + (kproc-1)*(q+1)
  else
     nz_   = q
     kmin_ = kmin + r*(q+1) + (kproc-r-1)*q
  end if
  nzo_   = nz_ + 2*nover
  kmax_  = kmin_ + nz_ - 1
  kmino_ = kmin_ - nover
  kmaxo_ = kmax_ + nover
  
  ! Initialize the memory module
  call memory_init(imino_,imaxo_,jmino_,jmaxo_,kmino_,kmaxo_,nover)
  
  return
end subroutine partition_init
