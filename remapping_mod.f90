module remap
use sys_state
use const
use trafo

contains


subroutine remap_stepwise()
  call set_ik_bar(sheartime)
  state%u_f%val(:,:,1) = remap_brucker(state%u_f%val(:,:,1))
  state%u_f%val(:,:,2) = remap_brucker(state%u_f%val(:,:,2))

  if(abs(T_rm/2.0-sheartime)<abs((T_rm/2.0_rp)-(sheartime+dt))) then
    sheartime = -T_rm/2.0_rp
    call set_ik_bar(sheartime)
  end if
end subroutine

function remap_brucker(in_arr)
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1)            :: in_arr,remap_brucker
    !remap_brucker = dealiase_field(in_arr)
    remap_brucker = in_arr
    if(abs(T_rm/2.0-sheartime)<abs((T_rm/2.0_rp)-(sheartime+dt))) then
    write(*,*) '---sheartime reached---'

    ! set regions to zero before they are carried out of resolvable domain by remeshing 
    do i =0,xdim-1 
      do j =0,ydim-1 
        if(     (aimag(state%iky%val(i,j)) >= maxval(aimag(state%iky_bar%val(i,:)))) &
            .OR.(aimag(state%iky%val(i,j)) <= minval(aimag(state%iky_bar%val(i,:))))) then
          remap_brucker(i,j) = cmplx(0.0_rp,0.0_rp)
        end if
      end do
    end do
    !----------------------------------------------------------------------------------------
    !STEP 1: Transform back to realspace in y-direction!
	  do i=0,xdim-1
	  	y_pen_f = remap_brucker(i,:)										
      	call dfftw_execute_dft(yf_y,y_pen_f, y_pen)						
	  	remap_brucker(i,:) = y_pen									
	  end do	
    !STEP 2: Phaseshift back to -T_rm/2.0_rp
	 	do j=1,ydim-1							
	         do i=0,xdim-1
	 		    	!multiply the fourier spectrum with corresponding inverse phase factor
	 		    	remap_brucker(i,j) = remap_brucker(i,j)*exp(-shear*(T_rm)*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
	 		    end do	
	 	end do	
    !STEP 3: transform to sheared Fourier-space on new base
	  do i=0,xdim-1
	  	y_pen = remap_brucker(i,:)									
      	call dfftw_execute_dft(y_yf,y_pen, y_pen_f)					
	  	remap_brucker(i,:) = y_pen_f/real(ydim,rp)			
	  end do	
    !----------------------------------------------------------------------------------------
    ! set regions to zero wich have been carried into the resolutable domain by remapping
    call set_ik_bar(-T_rm/2.0_rp)
    !do i =0,xdim-1 
    !  do j =0,ydim-1 
    !    if(     (aimag(state%iky%val(i,j)) >= maxval(aimag(state%iky_bar%val(i,:)))) &
    !        .OR.(aimag(state%iky%val(i,j)) <= minval(aimag(state%iky_bar%val(i,:))))) then
    !      remap_brucker(i,j) = cmplx(0.0_rp,0.0_rp)
    !    end if
    !  end do
    !end do
    !remap_brucker = dealiase_field(remap_brucker)
    else
      ! set all modes to zero in brucker space which can not be resolved on the real space grid
      remap_brucker = in_arr
      call surgery(remap_brucker)
    end if
end function


subroutine surgery(patient)
    ! delete those parts of spectrum wich can not be resolved
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(inout)            :: patient

      call set_ik_bar(sheartime)
      do i =0,xdim-1 
        do j =0,ydim-1 
          if((aimag(state%iky_bar%val(i,j)) > ky_max).OR.(aimag(state%iky_bar%val(i,j)) <ky_min)) then
            patient(i,j) = cmplx(0.0_rp,0.0_rp)
          end if
        end do
      end do
        ! set all modes in brucker space to zero which can  be resolved on the real space grid but not in brucker space
      !do i =0,xdim-1 
      !  do j =0,ydim-1 
      !    if(     (aimag(state%iky%val(i,j)) >= maxval(aimag(state%iky_bar%val(i,:)))) &
      !        .OR.(aimag(state%iky%val(i,j)) <= minval(aimag(state%iky_bar%val(i,:))))) then
      !      patient(i,j) = cmplx(0.0_rp,0.0_rp)
      !    end if
      !  end do
      !end do
end subroutine

end module
   ! ! TODO: inefficient, since all transforms are conducted, only one dir necessary
   ! call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,T_rm/2.0_rp)
   ! call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,T_rm/2.0_rp)
   ! call transform(state%temp_f%val,state%temp%val,-1,shearing,T_rm/2.0_rp)
   ! call transform(state%chem_f%val,state%chem%val,-1,shearing,T_rm/2.0_rp)
   ! sheartime = -T_rm/2.0_rp
   ! call transform(state%u%val(:,:,1),state%u_f%val(:,:,1),+1,shearing,sheartime)
   ! call transform(state%u%val(:,:,2),state%u_f%val(:,:,2),+1,shearing,sheartime)
   ! call transform(state%temp%val,state%temp_f%val,+1,shearing,sheartime)
   ! call transform(state%chem%val,state%chem_f%val,+1,shearing,sheartime)
