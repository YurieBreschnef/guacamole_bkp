module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()
  !performs remapping of those fourier space components which will be carried out of the 
  !resolutable area of fourier space 
    real(kind = rp)                              :: ky_max
  
  state_nm1%u_f     = state%u_f
  state_nm1%temp_f  = state%temp_f
  state_nm1%chem_f  = state%chem_f

   call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
   call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
   call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
   call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)
   sheartime = 0.0_rp
   call transform(state%u%val(:,:,1),state%u_f%val(:,:,1),1,shearing,sheartime)
   call transform(state%u%val(:,:,2),state%u_f%val(:,:,2),1,shearing,sheartime)
   call transform(state%temp%val,state%temp_f%val,1,shearing,sheartime)
   call transform(state%chem%val,state%chem_f%val,1,shearing,sheartime)

  ky_max = abs(aimag(state%iky%val(1,ydim/2)))
  do i=0,xdim/2
    do j=0,ydim-1 
     if(abs(real(aimag(state%iky%val(i,j)-shear*T_rm*state%ikx%val(i,j)),rp))>=ky_max) then
      write(*,*) 'remapping u_f(',i,',',j,') to u_f(',i,',',mod(j-i,ydim),')'
      state%u_f%val(i,j,1) = state_nm1%u_f%val(   i,mod(j-i,ydim),1) 
      state%u_f%val(i,j,2) = state_nm1%u_f%val(   i,mod(j-i,ydim),2) 
      state%temp_f%val(i,j)= state_nm1%temp_f%val(i,mod(j-i,ydim))
      state%chem_f%val(i,j)= state_nm1%chem_f%val(i,mod(j-i,ydim))
     end if
    end do
  end do

  ky_max = abs(aimag(state%iky%val(1,ydim/2+1)))
  do i=xdim/2+1,xdim-1
    do j=0,ydim-1 
     if(abs(real(aimag(state%iky%val(i,j)-shear*T_rm*state%ikx%val(i,j)),rp))>=ky_max) then
      write(*,*) 'remapping u_f(',i,',',j,') to u_f(',i,',',mod(j+i,ydim),')'
      state%u_f%val(i,j,1) = state_nm1%u_f%val(   i,mod(j+i,ydim),1) 
      state%u_f%val(i,j,2) = state_nm1%u_f%val(   i,mod(j+i,ydim),2) 
      state%temp_f%val(i,j)= state_nm1%temp_f%val(i,mod(j+i,ydim))
      state%chem_f%val(i,j)= state_nm1%chem_f%val(i,mod(j+i,ydim))
     end if
    end do
  end do

end subroutine
end module
