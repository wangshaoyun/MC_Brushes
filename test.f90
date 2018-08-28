program test

  implicit none
  real*8 :: st, fn, x(3),a(3)
  integer (kind=8)  i
  complex(kind=8) :: b(3),c
  a=(/1.21D0,1.44D0,1.69D0/)
  b(1) = (1.,2)
  b(2) = (2.,1)
  b(3) = (1.,1)
  c=(1.232D0,2.345D0)
  i=1234567890
  call cpu_time(st)
  write(*,*) real(b),c,i,sqrt(a)
  call cpu_time(fn)


end program test
