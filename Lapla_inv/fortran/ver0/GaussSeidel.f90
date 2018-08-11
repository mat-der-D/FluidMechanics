module settings

  implicit none

  real(8), parameter :: pi = 3.1415926535897932384626d0
  real(8), parameter :: Lx = 50d0, Lz = 50d0
  integer, parameter :: Nx = 100, Nz = 100
  real(8), parameter :: dx = Lx/Nx, dz = Lz/Nz
  integer :: iii, kkk
  real(8), parameter :: &
       & x(0:Nx-1) = (/ (dx*iii, iii=0, Nx-1) /), &
       & z(0:Nz) = (/ (dz*kkk, kkk=0, Nz) /)
  
end module settings



program Lapla_inv_GS

  use settings
  
  implicit none

  real(8) :: u(0:Nx-1, 0:Nz), f(0:Nx-1, 0:Nz)
  real(8) :: u_exact(0:Nx-1, 0:Nz)
  real(8) :: xx(0:Nx-1, 0:Nz), zz(0:Nx-1, 0:Nz)
  integer :: i, k, it

  do k = 0, Nz
     xx(:,k) = x(:)
  end do

  do i = 0, Nx-1
     zz(i,:) = z(:)
  end do
  
  call initialize(u, f, u_exact)

  ! time development
  open(20, file="error.dat")
  
  do it = 1, 10000
     u = u_GaussSeidel_u(u, f)
     write(20, *) it, dist(f, u_Lapla_u(u))
  end do

  close(20)

  ! output
  open(10, file="test.dat")

  call output_u(f, u_Lapla_u(u), 10)
  
  close(10)
  
contains

  subroutine initialize(u, f, u_exact)

    real(8) :: u(0:Nx-1, 0:Nz), f(0:Nx-1, 0:Nz)
    real(8) :: u_exact(0:Nx-1, 0:Nz)

    ! f = sin(2*pi*zz/Lz)
    f = exp(zz/(0.3*Lz)) * cos(4*pi*xx/Lx)
    f(:,0) = 0d0
    f(:,Nz) = 0d0
    ! u_exact = -(Lz/(2*pi))**2*sin(2*pi*zz/Lz)
    u_exact = 0d0
    u = 0d0
    ! u = -63 * f
    
  end subroutine initialize

  function dist(u1, u2)

    real(8) :: dist
    real(8) :: u1(0:Nx-1, 0:Nz), u2(0:Nx-1, 0:Nz)
    integer :: i, k
    real(8) :: sum_tmp

    ! === Pattern I ===
    
    ! sum_tmp = 0d0
    ! do i = 0, Nx-1
    !    do k = 0, Nz
    !       sum_tmp = sum_tmp + (u1(i,k) - u2(i,k))**2
    !    end do
    ! end do

    ! write(*,*) sum_tmp
    ! dist = sqrt( sum_tmp / (Nx * (Nz + 1)) )

    ! === Pattern II ===
    dist = sqrt( sum((u1 - u2)**2) / (Nx * (Nz + 1)) )
    
  end function dist

  subroutine output_u(u, u_exact, unit)

    integer :: unit
    real(8) :: u(0:Nx-1, 0:Nz), u_exact(0:Nx-1, 0:Nz)
    integer, parameter :: igrid = 1, kgrid = 1
    integer :: i, k

    do i = 0, Nx-1, igrid
       do k = 0, Nz, kgrid
          ! write(unit, *) x(i), z(k), u(i,k)
          write(unit, *) x(i), z(k), u(i,k), u_exact(i,k)
       end do
    end do
        
  end subroutine output_u
  
  function next_u(i,k,u,f)

    ! output
    real(8) :: next_u
    ! input
    integer :: i, k
    real(8) :: u(0:Nx-1, 0:Nz), f(0:Nx-1, 0:Nz)

    if (k == 0 .or. k == Nz) then
       next_u = 0d0
    else
       next_u = - f(i,k)/( 2*(dx**(-2) + dz**(-2)) )&
            & + (u(inc_x(i), k) + u(dec_x(i), k))/&
            &   ( 2*(1 + (dx/dz)**2) )&
            & + (u(i, k+1) + u(i, k-1))/&
            &   ( 2*(1 + (dz/dx)**2) )
    end if
    
  end function next_u
  
  function u_GaussSeidel_u(u, f)

    ! output
    real(8) :: u_GaussSeidel_u(0:Nx-1, 0:Nz)
    ! input
    real(8) :: u(0:Nx-1, 0:Nz), f(0:Nx-1, 0:Nz)
    ! working variables
    integer :: i, k

    u_GaussSeidel_u = u
    
    do k = 0, Nz
       do i = 0, Nx-1
          u_GaussSeidel_u(i,k) = next_u(i,k,u_GaussSeidel_u,f)
       end do
    end do
    
  end function u_GaussSeidel_u
  
  function dec_x(num)
    integer :: dec_x, num

    dec_x = mod(num - 1 + Nx, Nx)
    
  end function dec_x

  function inc_x(num)
    integer :: inc_x, num

    inc_x = mod(num + 1, Nx)
    
  end function inc_x

  function u_Lapla_u(u)

    real(8) :: u_Lapla_u(0:Nx-1, 0:Nz)
    real(8) :: u(0:Nx-1, 0:Nz)
    real(8) :: u_xx(0:Nx-1, 0:Nz), u_zz(0:Nx-1, 0:Nz)

    u_xx = (cshift(u, +1, 1) - 2*u + cshift(u, -1, 1))/(dx**2)
    u_zz = (cshift(u, +1, 2) - 2*u + cshift(u, -1, 2))/(dz**2)
    u_Lapla_u = u_xx + u_zz
    u_Lapla_u(:,0) = 0d0
    u_Lapla_u(:,Nz) = 0d0
    
  end function u_Lapla_u


end program Lapla_inv_GS
