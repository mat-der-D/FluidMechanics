module settings

  implicit none
  ! mathematical const.
  real(8), parameter :: pi = 3.1415926535897932384626d0
  ! size parameters
  real(8), parameter :: Lx = 50d0, Lz = 50d0
  integer, parameter :: Nx = 200, Nz = 200
  integer, parameter :: dx_out = 10, dz_out = 10
  real(8), parameter :: dx = Lx/Nx, dz = Lz/Nz
  ! spectral parameters
  integer, parameter :: lm = 64 !cut-off wave number
  
  ! time parameters
  real(8), parameter :: dt = 1d-2
  integer, parameter :: Nt = 10000, dt_out = 100
  ! physical parameters
  real(8), parameter :: D = 1d0
  real(8), parameter :: lmb = D * dt / (2 * dz**2)
  
end module settings

program diffusion2D_spec_diff

  use settings
  use ae_module
  use lumatrix

  !==== Declarative Part ======================================
  !--- physical quantities ---
  real(8), dimension(0:Nz, 0:Nx-1) :: zx_u, zx_u_exact
  real(8), dimension(0:Nz, -lm:lm) :: ze_u, ze_u_exact

  real(8), dimension(0:Nx-1) :: x_uBCU, x_uBCL
  real(8), dimension(-lm:lm) :: e_uBCU, e_uBCL

  !--- linear equations ---
  real(8), dimension(0:Nz, 0:Nz) :: zz_Adiff
  integer :: kp(0:Nz)

  !--- basic variables ---
  real(8), dimension(0:Nz, 0:Nx-1) :: zx_x, zx_z

  integer :: i, k
  integer :: it
  
  !==== Main ==================================================

  !---- initialization ----
  call ae_initial(Nx, lm, 0d0, Lx)
  do i = 0, Nx-1
     do k = 0, Nz
        zx_x(k, i) = i*dx
        zx_z(k, i) = k*dz
     end do
  end do
  call BC_initial
  call value_initial
  call set_BC(zx_u, x_uBCL, x_uBCU)
  call set_BC(ze_u, e_uBCL, e_uBCU)
  call diff_mat_initial(zz_Adiff)

  call LUDecomp(zz_Adiff, kp)


  open(10, file="test.d")
  open(20, file="test_exact.d")
  
  call output_zx(zx_u, 10)
  call output_zx(zx_u, 20)
  do it = 1, Nt

     call ze_evolvediff_ze(ze_u, e_uBCL, e_uBCU, zz_Adiff, kp)

     if (mod(it, dt_out) == 0) then
        zx_u = ag_ae(ze_u)
        zx_u_exact = zx_exact(it * dt)
        call output_zx(zx_u, 10)
        call output_zx(zx_u_exact, 20)
        write(*,*) "t=", it * dt
     end if
  end do

  close(10)
  close(20)
  
contains
  !==== Functions & Subroutines ===============================
  subroutine value_initial

    ! zx_u = sin(2 * pi * zx_x / Lx) * sin(pi * zx_z / Lz) &
    !      & + sin((pi / 2) * zx_z / Lz)

    zx_u = sin(2 * pi * zx_x / Lx) * sin(pi * zx_z / Lz)

    !--- transform ---
    ze_u = ae_ag(zx_u)
    
  end subroutine value_initial

  function zx_exact(t)

    real(8) :: zx_exact(0:Nz, 0:Nx-1)
    real(8) :: t, dec_rate

    dec_rate = D*( (2*pi / Lx)**2 + (pi / Lz)**2 )
    zx_exact = sin(2 * pi * zx_x / Lx) * sin(pi * zx_z / Lz) &
         & * exp(- dec_rate * t)

  end function zx_exact
  
  subroutine diff_mat_initial(zz_A)

    real(8) :: zz_A(0:Nz, 0:Nz)
    real(8) :: zz_E(0:Nz, 0:Nz)

    integer :: ii
    
    zz_E = 0d0
    do ii = 0, Nz
       zz_E(ii, ii) = 1d0
    end do
    
    zz_A = (1 + 2*lmb)*zz_E &
         & - lmb * (cshift(zz_E, 1, 1) + cshift(zz_E, -1, 1))
    zz_A(0,:) = zz_E(0,:)
    zz_A(Nz,:) = zz_E(Nz,:)

  end subroutine diff_mat_initial
  
  subroutine BC_initial
    
    x_uBCU = 0d0
    x_uBCL = 0d0
    e_uBCU = e_g(x_uBCU)
    e_uBCL = e_g(x_uBCL)
    
  end subroutine BC_initial

  subroutine set_BC(func, funcBCL, funcBCU)

    real(8), dimension(:,:) :: func
    real(8), dimension(size(func, 2)) :: funcBCL, funcBCU

    if ( size(func, 1) /= Nz + 1 ) then

       write(*,*) "ERROR!: set_BC"
       stop
       
    end if
    
    func(1, :) = funcBCL(:)
    func(Nz+1,:) = funcBCU(:)
    
  end subroutine set_BC
  
  function mu_l(l)

    real(8) :: mu_l
    integer :: l
    real(8) :: alp_l

    alp_l = 2 * pi * l / Lx
    mu_l = exp(D * alp_l**2 * dt)

  end function mu_l

  subroutine ze_evolvediff_ze(ze, e_BCL, e_BCU, zz_LU, kp)

    real(8), dimension(0:Nz, -lm:lm) :: ze
    real(8), dimension(-lm:lm) :: e_BCL, e_BCU, e_BCL2, e_BCU2
    real(8), dimension(0:Nz, 0:Nz) :: zz_LU
    integer, dimension(0:Nz) :: kp
    
    real(8), dimension(-lm:lm, 0:Nz) :: ez_workrhs, ez_mid
    
    real(8), dimension(0:Nz, -lm:lm) :: ze_workrhs, ze_mid
    integer :: l

    do l = -lm, lm
       e_BCL2(l) = mu_l(l) * e_BCL(l)
       e_BCU2(l) = mu_l(l) * e_BCU(l)
    end do
    
    ze_workrhs = (1 - 2*lmb)*ze &
         & + lmb*(cshift(ze, 1, 1) + cshift(ze, -1, 1))
    call set_BC(ze_workrhs, e_BCL2, e_BCU2)

    ez_workrhs = transpose(ze_workrhs)

    ez_mid = LUSolve(zz_LU, kp, ez_workrhs)
    ze_mid = transpose(ez_mid)

    do l = -lm, lm
       ze(:,l) = ze_mid(:,l) / mu_l(l)
    end do
    
  end subroutine ze_evolvediff_ze

  subroutine output_zx(zx, unit)

    real(8), dimension(0:Nz, 0:Nx-1) :: zx
    integer :: unit
    integer :: i, k

    do i = 0, Nx-1, dx_out
       do k = 0, Nz, dz_out
          write(unit, *) &
               zx_x(k, i), zx_z(k, i), zx(k, i)
       end do
       write(unit, *)
    end do
    write(unit, *)
    
  end subroutine output_zx
  
end program diffusion2D_spec_diff
