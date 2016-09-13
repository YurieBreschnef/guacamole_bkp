module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()
  !performs remapping of those fourier space components which will be carried out of the 
  !resolutable area of fourier space 
    real(kind = rp)                              :: ky_max
  if(state%t >=T_rm)   then
    state_nm1%u_f%val     = state%u_f%val
    state_nm1%temp_f%val  = state%temp_f%val
    state_nm1%chem_f%val  = state%chem_f%val

   call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
   call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
   call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
   call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)
   sheartime = 0.0_rp
   call transform(state%u%val(:,:,1),state%u_f%val(:,:,1),1,shearing,sheartime)
   call transform(state%u%val(:,:,2),state%u_f%val(:,:,2),1,shearing,sheartime)
   call transform(state%temp%val,state%temp_f%val,1,shearing,sheartime)
   call transform(state%chem%val,state%chem_f%val,1,shearing,sheartime)
  end if

open(unit=20,file="remap_vec.dat",status='replace',action='write',iostat=io_error) 
  ky_max = abs(aimag(state%iky%val(1,ydim/2)))
  do i=0,xdim/2
    do j=1,ydim-1 
     if(abs(real(aimag(state%iky%val(i,j)-shear*T_rm*state%ikx%val(i,j)),rp))>=ky_max) then
      state%u_f%val(i,j,1) =cmplx(0.0_rp,0.0_rp,rp) 
      state%u_f%val(i,j,2) =cmplx(0.0_rp,0.0_rp,rp) 
      state%temp_f%val(i,j)=cmplx(0.0_rp,0.0_rp,rp) 
      state%chem_f%val(i,j)=cmplx(0.0_rp,0.0_rp,rp) 
     end if
    end do
  end do

  ky_max = abs(aimag(state%iky%val(1,ydim/2+1)))
  do i=xdim/2+1,xdim-1
    do j=1,ydim-1 
     if(abs(real(aimag(state%iky%val(i,j)-shear*T_rm*state%ikx%val(i,j)),rp))>=ky_max) then
      state%u_f%val(i,j,1) =cmplx(0.0_rp,0.0_rp,rp) 
      state%u_f%val(i,j,2) =cmplx(0.0_rp,0.0_rp,rp) 
      state%temp_f%val(i,j)=cmplx(0.0_rp,0.0_rp,rp) 
      state%chem_f%val(i,j)=cmplx(0.0_rp,0.0_rp,rp) 
     end if
    end do
  end do
















!open(unit=20,file="remap_vec.dat",status='replace',action='write',iostat=io_error) 
!  ky_max = abs(aimag(state%iky%val(1,ydim/2)))
!  do i=0,xdim/2
!    do j=1,ydim-1 
!     if(abs(real(aimag(state%iky%val(i,j)-shear*T_rm*state%ikx%val(i,j)),rp))>=ky_max) then
!      !write(*,*) 'remapping u_f(',i,',',j,') to u_f(',i,',',mod(j-i,ydim),')'
!      write(20,*) i,j,i-i,mod(j-i,ydim)-j
!      state%u_f%val(i,j,1) = state_nm1%u_f%val(   i,mod(j-i,ydim),1) 
!      state%u_f%val(i,j,2) = state_nm1%u_f%val(   i,mod(j-i,ydim),2) 
!      state%temp_f%val(i,j)= state_nm1%temp_f%val(i,mod(j-i,ydim))
!      state%chem_f%val(i,j)= state_nm1%chem_f%val(i,mod(j-i,ydim))
!     end if
!    end do
!  end do
!
!  ky_max = abs(aimag(state%iky%val(1,ydim/2+1)))
!  do i=xdim/2+1,xdim-1
!    do j=1,ydim-1 
!     if(abs(real(aimag(state%iky%val(i,j)-shear*T_rm*state%ikx%val(i,j)),rp))>=ky_max) then
!      !write(*,*) 'remapping u_f(',i,',',j,') to u_f(',i,',',mod(j+ydim-i,ydim),')'
!      write(20,*) i,j,i-i,mod(j+ydim-i,ydim)-j
!      state%u_f%val(i,j,1) = state_nm1%u_f%val(   i,mod((j+ydim-i),ydim),1) 
!      state%u_f%val(i,j,2) = state_nm1%u_f%val(   i,mod((j+ydim-i),ydim),2) 
!      state%temp_f%val(i,j)= state_nm1%temp_f%val(i,mod((j+ydim-i),ydim))
!      state%chem_f%val(i,j)= state_nm1%chem_f%val(i,mod((j+ydim-i),ydim))
!     end if
!    end do
!  end do

 ! write(*,*) '------remapped------'
    close(20)
end subroutine

























































end module
