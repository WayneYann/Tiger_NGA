module memory
  use precision
  implicit none
  
  ! Arrays to be used for computation in solver:
  ! -> velocity
  ! -> scalar_quick
  ! -> SGS_lagrangian
  ! -> combustion
  real(WP), dimension(:,:,:), pointer :: tmp1
  real(WP), dimension(:,:,:), pointer :: tmp2
  real(WP), dimension(:,:,:), pointer :: tmp3
  real(WP), dimension(:,:,:), pointer :: tmp4
  real(WP), dimension(:,:,:), pointer :: tmp5
  real(WP), dimension(:,:,:), pointer :: tmp6
  real(WP), dimension(:,:,:), pointer :: tmp7
  real(WP), dimension(:,:,:), pointer :: tmp8
  real(WP), dimension(:,:,:), pointer :: tmp9
  real(WP), dimension(:,:,:), pointer :: tmp10
  real(WP), dimension(:,:,:), pointer :: tmp11
  real(WP), dimension(:,:,:), pointer :: tmp12
  real(WP), dimension(:,:,:), pointer :: tmp13
  real(WP), dimension(:,:,:), pointer :: tmp14
  real(WP), dimension(:,:,:), pointer :: tmp15
  
  ! The following arrays are just pointer to the previous arrays
  real(WP), dimension(:,:,:), pointer :: FX
  real(WP), dimension(:,:,:), pointer :: FY
  real(WP), dimension(:,:,:), pointer :: FZ
  real(WP), dimension(:,:,:), pointer :: rhoUi
  real(WP), dimension(:,:,:), pointer :: rhoVi
  real(WP), dimension(:,:,:), pointer :: rhoWi
  real(WP), dimension(:,:,:), pointer :: Fcyl
  real(WP), dimension(:,:,:), pointer :: Fcylv
  
  ! For the ADI
  real(WP), dimension(:,:,:),   pointer :: Rx,Ry,Rz
  real(WP), dimension(:,:,:,:), pointer :: Ax,Ay,Az
  
  ! Work array for multidiagonal solvers
  real(WP), dimension(:,:), pointer :: stackmem
  
end module memory


subroutine memory_init(i1,i2,j1,j2,k1,k2,n)
  use memory
  implicit none
  integer, intent(in) :: i1,i2,j1,j2,k1,k2,n
  
  ! Allocate the main temporary arrays
  allocate(tmp1 (i1:i2,j1:j2,k1:k2))
  allocate(tmp2 (i1:i2,j1:j2,k1:k2))
  allocate(tmp3 (i1:i2,j1:j2,k1:k2))
  allocate(tmp4 (i1:i2,j1:j2,k1:k2))
  allocate(tmp5 (i1:i2,j1:j2,k1:k2))
  allocate(tmp6 (i1:i2,j1:j2,k1:k2))
  allocate(tmp7 (i1:i2,j1:j2,k1:k2))
  allocate(tmp8 (i1:i2,j1:j2,k1:k2))
  allocate(tmp9 (i1:i2,j1:j2,k1:k2))
  allocate(tmp10(i1:i2,j1:j2,k1:k2))
  allocate(tmp11(i1:i2,j1:j2,k1:k2))
  allocate(tmp12(i1:i2,j1:j2,k1:k2))
  allocate(tmp13(i1:i2,j1:j2,k1:k2))
  allocate(tmp14(i1:i2,j1:j2,k1:k2))
  allocate(tmp15(i1:i2,j1:j2,k1:k2))
  
  ! Allocate arrays for tri/penta-diagonal solvers
  allocate(stackmem((i2-i1)*(j2-j1)*(k2-k1),2*n+1))
  
  ! For the approximate factorization
  allocate(Ax(j1+n:j2-n,k1+n:k2-n,i1+n:i2-n,-n:+n))
  allocate(Ay(i1+n:i2-n,k1+n:k2-n,j1+n:j2-n,-n:+n))
  allocate(Az(i1+n:i2-n,j1+n:j2-n,k1+n:k2-n,-n:+n))
  
  allocate(Rx(j1+n:j2-n,k1+n:k2-n,i1+n:i2-n))
  allocate(Ry(i1+n:i2-n,k1+n:k2-n,j1+n:j2-n))
  allocate(Rz(i1+n:i2-n,j1+n:j2-n,k1+n:k2-n))
  
  ! Create the link for some nicer named pointers
  rhoUi => tmp1
  rhoVi => tmp2
  rhoWi => tmp3
  FX    => tmp4
  FY    => tmp5
  FZ    => tmp6
  Fcyl  => tmp7
  Fcylv => tmp8
  
  return
end subroutine memory_init
