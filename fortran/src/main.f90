program main
  use adjust_mod, only : adjust
  use utils,  only : logspace, linspace
  
  implicit none
  
  integer, parameter :: npz = 50
  real, dimension(npz+1) :: pedge
  real, dimension(npz) :: p, T, q, dp, Told
  logical, dimension(npz) :: mask
  integer :: k, file, iostat

  call logspace(1.,5.,pedge)
  call linspace(-600.,300., T)

  do k=1,npz
     dp(k) = pedge(k+1) - pedge(k)
     p(k) = dp(k)/ (log(pedge(k+1)/pedge(k)))
     if (T(k) .lt. 150) T(k) = 150.
     q(k) = 0.1
  enddo

  
  Told = T
  call adjust(p, dp, T, q, mask)

  file = 1
  open(1, file='output.txt')
  do k=1,npz
     write(1,*) p(k), Told(k), T(k), mask(k)
  enddo
  close(1)
 
end program main
