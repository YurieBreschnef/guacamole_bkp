module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()
  !performs remapping of those fourier space components which will be carried out of the 
  !resolutable area of fourier space 
 
  if(mod(state%step,remapping_rate)==0) then
    ! TODO: inefficient, since all transforms are conducted, only on dir necessary
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
    call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
    call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)
    sheartime = 0.0_rp
    call transform(state%u%val(:,:,1),state%u_f%val(:,:,1),+1,shearing,sheartime)
    call transform(state%u%val(:,:,2),state%u_f%val(:,:,2),+1,shearing,sheartime)
    call transform(state%temp%val,state%temp_f%val,+1,shearing,sheartime)
    call transform(state%chem%val,state%chem_f%val,+1,shearing,sheartime)
  end if

end subroutine

end module
