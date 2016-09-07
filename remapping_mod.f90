module remap
use sys_state
use const

contains


subroutine remap_stepwise()
  !performs remapping of those fourier space components which have been carried out of the 
  !resolutable zone
  integer                           :: x_rm_width = xdim/2-xdim/3
  integer                           :: y_rm_width = ydim/2-ydim/3

  state%dummy_f =state%u_f
  state%s_dummy_f =state%temp_f
  state%c_dummy_f =state%chem_f

  if(shearing == 0 .OR. remapping ==0) then
    write(*,*) 'ERROR: bad config. remapping has been called although shearing or remap is off.'
  end if
  

    !do i=xdim/3,2*(xdim/3)
    !    state%u_f%val(i,:,1) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
    !    state%u_f%val(i,:,2) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
    !    state%temp_f%val(i,:) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
    !    state%chem_f%val(i,:) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
	  !end do	
    !do j=ydim/3,2*(ydim/3)
    !    state%u_f%val(:,j,1) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
    !    state%u_f%val(:,j,2) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
    !    state%temp_f%val(:,j) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
    !    state%chem_f%val(:,j) = cmplx(real((rand()-0.5_rp),rp),real((rand()-0.5_rp),rp),rp)
	  !end do	

    !do i=0,x_rm_width
    !       state%u_f%val(xdim/3+i,:,1)  =    state%dummy_f%val(2*xdim/3-i+1,:,1)
    !       state%u_f%val(xdim/3+i,:,2)  =    state%dummy_f%val(2*xdim/3-i+1,:,2)
    !    state%temp_f%val(xdim/3+i,:)    =    state%s_dummy_f%val(2*xdim/3-i+1,:)  
    !    state%chem_f%val(xdim/3+i,:)    =    state%c_dummy_f%val(2*xdim/3-i+1,:)  
	  !end do	
    !do i=0,x_rm_width
    !    state%u_f%val(xdim/2+i+1,:,1)  =    state%dummy_f%val(xdim/2-i,:,1)
    !    state%u_f%val(xdim/2+i+1,:,2)  =    state%dummy_f%val(xdim/2-i,:,2)
    !    state%temp_f%val(xdim/2+i+1,:)    = state%s_dummy_f%val(xdim/2-i,:)  
    !    state%chem_f%val(xdim/2+i+1,:)    = state%c_dummy_f%val(xdim/2-i,:)  
	  !end do	
    do j=0,y_rm_width
           state%u_f%val(:,ydim/3+j,1)  =    state%dummy_f%val(:,2*ydim/3-j+1,1)
           state%u_f%val(:,ydim/3+j,2)  =    state%dummy_f%val(:,2*ydim/3-j+1,2)
        state%temp_f%val(:,ydim/3+j)    =    state%s_dummy_f%val(:,2*ydim/3-j+1)  
        state%chem_f%val(:,ydim/3+j)    =    state%c_dummy_f%val(:,2*ydim/3-j+1)  
	  end do	
    do j=0,y_rm_width
        state%u_f%val(:,ydim/2+j+1,1)  =    state%dummy_f%val(:,ydim/2-j,1)
        state%u_f%val(:,ydim/2+j+1,2)  =    state%dummy_f%val(:,ydim/2-j,2)
        state%temp_f%val(:,ydim/2+j+1)    =    state%s_dummy_f%val(:,ydim/2-j)  
        state%chem_f%val(:,ydim/2+j+1)    =    state%c_dummy_f%val(:,ydim/2-j)  
	  end do	

end subroutine
end module
