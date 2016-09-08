module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()
  !performs remapping of those fourier space components which will be carried out of the 
  !resolutable area of fourier space 
  !resolutable zone
  real(kind = rp)         ::ky_max  !maximum resolutable ky = k(:,ydim/2)
  ky_max = abs(aimag(state%iky%val(1,ydim/2)))

  do i =0,xdim-1
  do j =0,ydim-1
   if() 
  end do
  end do
end subroutine
end module
