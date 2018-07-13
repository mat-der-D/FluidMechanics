!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 二次元拡散方程式
! u_t = D \delta u
! を ADI 法により解く.
! [境界条件]
!    -- z 方向: 固定境界(u(Lz)=1, u(0)=0)
!    -- x 方向: 周期境界

module settings

  implicit none
  ! mathematical const.
  real(8), parameter :: pi = 3.1415926535897932384626d0
  ! size parameters
  real(8), parameter :: Lx = 50d0, Lz = 50d0
  integer, parameter :: Nx = 100, Nz = 100
  real(8), parameter :: dx = Lx/Nx, dz = Lz/Nz
  integer :: iii, kkk
  real(8), parameter :: &
       & x(0:Nx-1) = (/ (dx*iii, iii=0, Nx-1) /), &
       & z(0:Nz) = (/ (dz*kkk, kkk=0, Nz) /)
  ! time parameters
  real(8), parameter :: dt = 1d-2
  integer, parameter :: Nt = 1000
  ! physical parameters
  real(8), parameter :: D = 1d0
  
end module settings

program diffusion2D_ADI

  use settings

  implicit none

  real(8) :: u(0:Nx-1, 0:Nz)

contains

  subroutine init_mat(A_x, A_z)

    real(8) :: A_x(


  end subroutine init_mat


end program diffusion2D_ADI
