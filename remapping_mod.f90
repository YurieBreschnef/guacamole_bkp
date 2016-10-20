module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()

  if(sheartime>=T_rm) then
    call dealiase_all()
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
    sheartime = 0.0_rp !-T_rm/2.0_rp
    call transform(state%u%val(:,:,1),state%u_f%val(:,:,1), 1,shearing,sheartime)
    call transform(state%u%val(:,:,2),state%u_f%val(:,:,2), 1,shearing,sheartime)
    call dealiase_all()
  else
    call dealiase_all()
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
    call transform(state%u%val(:,:,1),state%u_f%val(:,:,1), 1,shearing,sheartime)
    call transform(state%u%val(:,:,2),state%u_f%val(:,:,2), 1,shearing,sheartime)
    call dealiase_all()
  end if
end subroutine

end module

    !call transform(state%u%val(:,:,1),state%u_f%val(:,:,1),+1,shearing,sheartime)
    !call transform(state%u%val(:,:,2),state%u_f%val(:,:,2),+1,shearing,sheartime)
    !call transform(state%temp%val,state%temp_f%val,+1,shearing,sheartime)
    !call transform(state%chem%val,state%chem_f%val,+1,shearing,sheartime)
