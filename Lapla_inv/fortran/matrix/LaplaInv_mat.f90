module settings

  implicit none

  real(8), parameter :: pi = 3.1415926535897932d0
  real(8), parameter :: Lx = 50d0, Lz = 50d0
  integer, parameter :: Nx = 100, Nz = 100
  real(8), parameter :: dx = Lx/Nx, dz = Lz/Nz
  integer :: iii, kkk
  real(8), parameter :: &
       & x(0:Nx-1) = (/ (dx*iii, iii=0, Nx-1) /), &
       & z(0:Nz) = (/ (dz*kkk, kkk=0, Nz) /)
  
end module settings



program Lapla_inv_mat

  use lumatrix
  use dc_message, only : MessageNotify
  use settings
  
  implicit none

  real(8) :: test(0:Nx-1, 0:Nz)
  real(8) :: zx(0:Nx-1, 0:Nz)

  zx = 0d0
  test = zx_LaplaInv_zx(zx, 'LK')

contains

  function zx_LaplaInv_zx(zx, cond, new)

    real(8), intent(in) :: zx(0:Nx-1, 0:Nz)
    real(8) :: zx_LaplaInv_zx(0:Nx-1, 0:Nz)

    character(len=2), intent(in), optional :: cond
    ! cond = 'XY'
    ! X: 上端の条件, Y:下端の条件
    ! X,Y に入るもの:
    !   R: 粘着条件
    !   F: 応力無し条件
    !   K: 運動学的条件
    !   (default: 'RR')

    logical, intent(in), optional :: new
    
    real(8), allocatable :: alu(:,:)
    integer, allocatable :: kp(:)

    real(8) :: zxmix_work(Nx*(Nz+1), Nx*(Nz+1))

    character(len=1) :: cond_u = 'R'
    character(len=1) :: cond_l = 'R'
    logical :: first = .true.
    logical :: new_matrix = .false.
    integer :: n
    save    :: alu, kp, first

    if( present(cond) ) then
       cond_u = cond(1:1)
       cond_l = cond(2:2)
       write(*,*) cond_u
       write(*,*) cond_l
       if( (cond_u /= 'R' .and. cond_u /= 'F' .and. cond_u /= 'K') .or. &
           (cond_l /= 'R' .and. cond_l /= 'F' .and. cond_l /= 'K') ) then
          call MessageNotify('E','zx_LaplaInv_zx','B.C. not supported')
       end if
    else
       cond_u = 'R'
       cond_l = 'R'
    end if

    if( .not. present(new)) then
       new_matrix = .false.
    else
       new_matrix = new
    end if

    


    
    write(*,*) cond_u, cond_l

    zx_LaplaInv_zx = 0d0
    
  end function zx_LaplaInv_zx

  function zx_DDx_zx(zx)

    real(8), intent(in) :: zx(0:Nx-1, 0:Nz)
    real(8) :: zx_DDx_zx(0:Nx-1, 0:Nz)

    zx_DDx_zx = ( cshift(zx,+1,1) - 2*zx + cshift(zx,-1,1) )/( dx**2 )
    
  end function zx_DDx_zx

  function zx_DDz_zx(zx)
    ! 出力の z = 0,Nz の値は無意味(0に固定)
    
    real(8), intent(in) :: zx(0:Nx-1, 0:Nz)
    real(8) :: zx_DDz_zx(0:Nx-1, 0:Nz)

    zx_DDz_zx = ( cshift(zx,+1,2) - 2*zx + cshift(zx,-1,2) )/( dz**2 )
    zx_DDz_zx(:,0) = 0d0
    zx_DDz_zx(:,Nz) = 0d0

  end function zx_DDz_zx
  
end program Lapla_inv_mat
