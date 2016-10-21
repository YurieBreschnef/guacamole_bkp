module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()
  !performs remapping of those fourier space components which will be carried out of the 
  !resolutable area of fourier space 
  call remap_brucker()
end subroutine




subroutine remap_brucker()
  ! TODO aliasing as in rogallo with zero padding 

    ! set all modes to zero in brucker space which can not be resolved on the real space grid
    do i =0,xdim-1 
      do j =0,ydim-1 
        if((aimag(state%iky_bar%val(i,j)) >= ky_max).OR.(aimag(state%iky_bar%val(i,j)) <=ky_min)) then
        state%u_f%val(i,j,:) = cmplx(0.0_rp,0.0_rp)
        end if
      end do
    end do

    ! set all modes in brucker space to zero which can  be resolved on the real space grid but not in brucker space
    do i =0,xdim-1 
      do j =0,ydim-1 
        if(     (aimag(state%iky%val(i,j)) >= maxval(aimag(state%iky_bar%val(i,:)))) &
            .OR.(aimag(state%iky%val(i,j)) <= minval(aimag(state%iky_bar%val(i,:))))) then
        state%u_f%val(i,j,:) = cmplx(0.0_rp,0.0_rp)
        end if
      end do
    end do

  if(sheartime>=T_rm/2.0) then
    
    write(*,*) '---sheartime reached---'
    ! TODO: inefficient, since all transforms are conducted, only one dir necessary
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,T_rm)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,T_rm)
    call transform(state%temp_f%val,state%temp%val,-1,shearing,T_rm)
    call transform(state%chem_f%val,state%chem%val,-1,shearing,T_rm)
    sheartime = -T_rm/2.0_rp
    call transform(state%u%val(:,:,1),state%u_f%val(:,:,1),+1,shearing,sheartime)
    call transform(state%u%val(:,:,2),state%u_f%val(:,:,2),+1,shearing,sheartime)
    call transform(state%temp%val,state%temp_f%val,+1,shearing,sheartime)
    call transform(state%chem%val,state%chem_f%val,+1,shearing,sheartime)
  end if

end subroutine

end module
