module nabla
  !implements certain nabla operations like curl, div etc
  use sys_state
  use const
  use trafo
  implicit none

  contains

  function crossp(arr1,arr2)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            ::arr1 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            ::arr2 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            ::crossp
  crossp(:,:,1) = cmplx(real(arr1(:,:,1),rp)*real(arr2(:,:,2),rp) &
                       -real(arr1(:,:,2),rp)*real(arr2(:,:,1),rp),0.0_rp)
  crossp(:,:,2) = cmplx(real(arr1(:,:,2),rp)*real(arr2(:,:,1),rp) &
                       -real(arr1(:,:,1),rp)*real(arr2(:,:,2),rp),0.0_rp)
    IF(ANY(IsNaN(real(crossp))))  then
      write(*,*) 'func crossp(): NAN detected in output array'
      stop
    end if
  end function 

end module
