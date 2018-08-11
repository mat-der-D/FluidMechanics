program test

  implicit none

  integer(8) :: a(3,3)

  a(1,:) = (/1, 2, 3/)
  a(2,:) = (/4, 5, 6/)
  a(3,:) = (/7, 8, 9/)

  a = cshift(a, +1, 2)

  write(*,*) a(1,:)
  write(*,*) a(2,:)
  write(*,*) a(3,:)

end program test
