module trafo
  ! performs the fourier transforms
  use sys_state
  use plans
  use benchmark
  implicit none

  contains

  subroutine dealiase_all()
    !dealiase the complete system state
     state%u_f%val(:,:,1) = dealiase_field(state%u_f%val(:,:,1)) 
     state%u_f%val(:,:,2) = dealiase_field(state%u_f%val(:,:,2)) 
     state%temp_f%val     = dealiase_field(state%temp_f%val) 
     state%chem_f%val     = dealiase_field(state%chem_f%val) 
  end subroutine

  function dealiase_field(arr_f)
    !dealiase a scalar field by 2/3rd rule
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in) :: arr_f
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1)            :: dealiase_field
    dealiase_field = arr_f
    do i=xdim/3,2*(xdim/3)
		  dealiase_field(i,:) = cmplx(0.0_rp,0.0_rp,rp)
	  end do
    do j=ydim/3,2*(ydim/3)
		  dealiase_field(:,j) = cmplx(0.0_rp,0.0_rp,rp)
	  end do	
end function


subroutine transform(in_arr,out_arr,dir,shearing,time)
  ! transforms a given 2D-array into (or out of) fourier space
  ! dir= 1 means real-->fourier
  ! dir=-1 means fourier-->real
  ! shearing=1 means shearing factor will be multiplied
	complex(kind = rp),dimension(0:xdim-1,0:ydim-1),intent(in)		:: in_arr
	complex(kind = rp),dimension(0:xdim-1,0:ydim-1),intent(inout)	:: out_arr 
	integer,intent(in)		                     										:: shearing
	integer,intent(in)		                     										:: dir
  real(kind= rp),intent(in)                                     :: time
  if(benchmarking ==1) call cpu_time(bm_trafo_starttime)
	if(debuglevel.GE.3) write(*,*) 'starting transform..'

  if(dir==1) then
          ! trafo in forward dir:
          ! x-pencils____________________________________________________________
	        do j=0,ydim-1													
	        	if(debuglevel.GE.3) write(*,*) 'transforming x-pencil number!:', j	
	        	x_pen = in_arr(:,j)										
            	call dfftw_execute_dft(x_xf,x_pen, x_pen_f)		
	        	out_arr(:,j) = x_pen_f/real(xdim,rp)		
	        end do	

	        ! PHASE FACTOR____________________________________________________________
	        if(shearing.EQ.1) then
	        	do j=1,ydim-1			
	        		if(debuglevel.GE.3) write(*,*) 'multiplying phase factor'
	        		do i=1,xdim-1			
	        			!out_arr(i,j) = out_arr(i,j)*exp(imag*shear*t*kxd(i)*(real(j,rp)/real(ydim,rp))*Ly)
	        			out_arr(i,j) = out_arr(i,j)*exp(shear*time*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
	        			!multiply the fourier spectrum with corresponding phase factor
	        		end do
	        	end do
	        end if
	        ! PHASE FACTOR___________________________________________________________

          ! y-pencils:
	        do i=0,xdim-1													
	        	if(debuglevel.GE.3) write(*,*) 'transforming y-pencil number:', i	
	        	y_pen = out_arr(i,:)									
            	call dfftw_execute_dft(y_yf,y_pen, y_pen_f)					
	        	out_arr(i,:) = y_pen_f/real(ydim,rp)			
	        end do	
	        if(debuglevel.GE.3) write(*,*) 'transform done.'

    else if(dir==-1) then
          ! transform in backward direction
          ! y-pencils
	        do i=0,xdim-1													
	        	if(debuglevel.GE.3) write(*,*) 'transforming y-pencil number:', i
	        	y_pen_f = in_arr(i,:)										
            	call dfftw_execute_dft(yf_y,y_pen_f, y_pen)						
	        	out_arr(i,:) = y_pen									
            ! norm only if dir=1
	        end do	
	        
	        ! PHASE FACTOR____________________________________________________________
	        if(shearing.EQ.1) then
	        	do j=1,ydim-1							
	        		if(debuglevel.GE.3) write(*,*) 'multiplying inverse phase factor'
	        		do i=1,xdim-1			
	        			out_arr(i,j) = out_arr(i,j)*exp(-shear*time*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
	        			!multiply the fourier spectrum with corresponding inverse phase factor
	        		end do	
	        	end do	
	        end if
	        ! PHASE FACTOR___________________________________________________________________

          ! x-pencils:
	        do j=0,ydim-1
	        	if(debuglevel.GE.3) write(*,*) 'transforming x-pencil number:', j
	        	x_pen_f = out_arr(:,j)											! fill x-pencil 
           		call dfftw_execute_dft(xf_x,x_pen_f, x_pen)						! transform
	        	out_arr(:,j) = x_pen											! fill u with pencil
	        end do	
        ! throw away residual complex values
        out_arr(:,:) = cmplx(real(out_arr(:,:),rp),0.0_rp,rp)
        
        else
          write(*,*) 'WARNING: sub trafo has been called with bad input parameter dir (!=(1 or -1))!'
  end if
  if(benchmarking ==1) call cpu_time(bm_trafo_endtime)
end subroutine

end module
