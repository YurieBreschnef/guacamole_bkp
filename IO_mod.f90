module IO_mod
  use sys_state
  use plans
  use nabla
  use trafo
  implicit none
  contains

  subroutine write_all()
    !if(debuglevel <= 2) write(*,*) '-calling write_all()'
    
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
    call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
    call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)

    call write_u()
    call write_abs_u()

    call write_chem()
    call write_temp()
    call write_buo()

    call write_div()
    call write_vort()

    call write_u_f()
    call write_chem_f()
    call write_chem_f_remap()
    call write_temp_f()
    call write_temp_f_remap()

    call write_u_stat()     ! u-relatetd measures
    call write_E_stat()     ! Energy related measures
    call write_sys_stat()   ! System wide measures

    !if(debuglevel <= 2) write(*,*) '-done with write_all.'
  end subroutine
  subroutine write_vort()
    !write chemical field to file. Note the path variable specifying the destination
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=19),parameter					:: path ='./output/data/vort/'
		!write(suffix,"(I5,A9)") state%step/(steps/maxfiles), ".vort.dat"
		write(suffix,"(I5,A9)") int(state%t/write_intervall), ".vort.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    !calc vorticity in fourier space
    state%s_dummy_f%val(:,:) =( state%ikx%val(:,:)*state%u_f%val(:,:,2) &
                             -state%iky%val(:,:)*state%u_f%val(:,:,1))

    call transform(state%s_dummy_f%val,state%s_dummy%val,-1,shearing,state%t)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_vort!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
                      ,real(state%s_dummy%val(i,j),real_outp_precision)
	  			           
			end do
		end do
    close(20)
  
  end subroutine

  
  subroutine write_u()
    !write velocity field to file in specified realtive path
    integer                             :: io_error = 0
    real(kind = real_outp_precision)    :: out_x,out_y,out_z
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=16),parameter					:: path ='./output/data/u/'
		write(suffix,"(I5,A6)") int(state%t/write_intervall) , ".u.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_u!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
          out_x = real(state%u%val(i,j,1),real_outp_precision)
          out_y = real(state%u%val(i,j,2),real_outp_precision)
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
                      ,out_x,out_y

			end do
		end do
    close(20)
  end subroutine

  subroutine write_chem()
    !write chemical field to file. Note the path variable specifying the destination
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=19),parameter					:: path ='./output/data/chem/'
		write(suffix,"(I5,A9)")  int(state%t/write_intervall), ".chem.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			            ,real(state%chem%val(i,j),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine

  subroutine write_abs_u()
    !write chemical field to file. Note the path variable specifying the destination
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=20),parameter					:: path ='./output/data/abs_u/'
		write(suffix,"(I5,A10)")  int(state%t/write_intervall), ".abs_u.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			            ,real(sqrt(real(state%u%val(i,j,1))**2+real(state%u%val(i,j,2))**2),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine
  subroutine write_temp()
    !write temp-field to file. analogous to write_chem
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=19),parameter					:: path ='./output/data/temp/'
		write(suffix,"(I5,A9)")  int(state%t/write_intervall), ".temp.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_temp!'
		  do i=0,xdim-1
        do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
                  ,real(state%temp%val(i,j),real_outp_precision)
			  !end do
			end do
		end do
    close(20)
  end subroutine

  subroutine write_div()
    ! writes divergence field to file
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=18),parameter					:: path ='./output/data/div/'
		write(suffix,"(I5,A8)")  int(state%t/write_intervall), ".div.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_div!'
    !write div u  (k_i * u_i ==0)condition to dummy
    state%s_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    !write brucker  (k_bar_i * u_i ==0)condition to other dummy
    state%c_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 
    call transform(state%c_dummy_f%val,state%cz_dummy_f%val,-1,1,state%t)
    call dfftw_execute_dft(ifull2D,state%s_dummy_f%val,state%s_dummy%val)
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			           ,real(state%s_dummy%val(i,j),real_outp_precision)&
                     ,real(state%cz_dummy_f%val(i,j),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine

  subroutine write_buo()
    !write buoyancy field to file
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=18),parameter					:: path ='./output/data/buo/'
		write(suffix,"(I5,A8)")  int(state%t/write_intervall), ".buo.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_buo!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			,real(state%temp%val(i,j)*B_therm &
              - B_comp*state%chem%val(i,j),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine
  
  subroutine write_u_f()
    ! writes the absolute amplitude of fourier coefficients at specified loc  to file
    integer                             :: io_error = 0
    real(kind = real_outp_precision)    :: out_x,out_y,out_z
		character(len=1024) 		  					:: filename
    type(sfield)                        :: x_r_dummy
    type(sfield)                        :: x_i_dummy
    type(sfield)                        :: y_r_dummy
    type(sfield)                        :: y_i_dummy
    type(vfield)                        :: dummy
		character(len=50) 		  						:: suffix
		character(len=18),parameter					:: path ='./output/data/u_f/'
		write(suffix,"(I5,A8)")  int(state%t/write_intervall), ".u_f.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    !dummy = state%u_f
		!do i=0,xdim-1
	  !	do j=0,ydim-1
	  ! 		dummy%val(i,j,1) = dummy%val(i,j,1)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
	  ! 		dummy%val(i,j,2) = dummy%val(i,j,2)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
    !  end do
    !end do
    !dummy = rearrange_2Dspectrum(deal_mask(dummy))
    x_r_dummy%val(:,:) = log(abs(real(state%u_f%val(:,:,1),real_outp_precision)))
    x_i_dummy%val(:,:) = log(abs(real(aimag(state%u_f%val(:,:,1)),real_outp_precision)))
    y_r_dummy%val(:,:) = log(abs(real(state%u_f%val(:,:,2),real_outp_precision)))
    y_i_dummy%val(:,:) = log(abs(real(aimag(state%u_f%val(:,:,2)),real_outp_precision)))

    x_r_dummy = rearrange_2Dspectrum(deal_mask(x_r_dummy))
    x_i_dummy = rearrange_2Dspectrum(deal_mask(x_i_dummy))
    y_r_dummy = rearrange_2Dspectrum(deal_mask(y_r_dummy))
    y_i_dummy = rearrange_2Dspectrum(deal_mask(y_i_dummy))

		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_u_f!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,real(x_r_dummy%val(i,j),real_outp_precision),&    !3: real(u_x)
                          real(x_i_dummy%val(i,j),real_outp_precision),&    !4: imag(u_x)
                          real(y_r_dummy%val(i,j),real_outp_precision),&    !5: real(u_y)
                          real(y_i_dummy%val(i,j),real_outp_precision)      !6: imag(u_y)
			end do
		end do
    close(20)
  end subroutine

  subroutine write_chem_f()
    !write fourier chem-field to file. Note path variable
    integer                             :: io_error = 0,remap
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/chem_f/'
		write(suffix,"(I5,A11)")  int(state%t/write_intervall), ".chem_f.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    dummy = state%chem_f
		!do i=0,xdim-1
	  !	do j=0,ydim-1
	  ! 		dummy%val(i,j) = dummy%val(i,j)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
    !  end do
    !end do
    dummy = rearrange_2Dspectrum(deal_mask(dummy))
    !dummy = deal_mask(state%chem_f)

		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem_f!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,log(abs(real(dummy%val(i,j),real_outp_precision))),&
                          log(abs(real(aimag(dummy%val(i,j)),real_outp_precision)))
			end do
		end do
    close(20)
  end subroutine

  subroutine write_chem_f_remap()
    !write fourier chem-field to file. Note path variable
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=27),parameter					:: path ='./output/data/chem_f_remap/'
		write(suffix,"(I5,A17)")  int(state%t/write_intervall), ".chem_f_remap.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    dummy%val = cmplx(1.0_rp,1.0_rp,rp)
    dummy%val = dealiase_field(dummy%val)

    do i =0,xdim-1
      do j =0,ydim-1
        if(dummy%val(i,j) .NE.cmplx(0.0_rp,0.0_rp,rp)) then
            dummy%val(i,j) = cmplx(0.0,0.0)
          else
          dummy%val(i,j) = state%chem_f%val(i,j)
        end if
      end do
    end do

    dummy = rearrange_2Dspectrum((dummy))
    !dummy = deal_mask(state%chem_f)
    !dummy = rearrange_2Dspectrum(deal_mask(dummy))

    !do i =0,xdim-1
    !  do j =0,ydim-1
    !    if((i>xdim/6 .AND.i<(5*xdim/6).AND.(j>ydim/6 .AND.j<(5*ydim/6)))) then
    !    dummy%val(i,j) =cmplx(0.0,0.0)
    !    end if
    !  end do
    !end do


		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem_f_remap!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,abs(real(dummy%val(i,j),real_outp_precision)),&
                          abs(real(aimag(dummy%val(i,j)),real_outp_precision))
			end do
		end do
    close(20)
  end subroutine

  subroutine write_temp_f_remap()
    !write fourier chem-field to file. Note path variable
    integer                             :: io_error = 0,remap
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=27),parameter					:: path ='./output/data/temp_f_remap/'
		write(suffix,"(I5,A17)")  int(state%t/write_intervall), ".temp_f_remap.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    dummy%val(:,:) = cmplx(1.0_rp,1.0_rp,rp)
    dummy%val = dealiase_field(dummy%val)

    do i =0,xdim-1
      do j =0,ydim-1
        if(dummy%val(i,j).NE.cmplx(0.0_rp,0.0_rp,rp)) then
            dummy%val(i,j) = cmplx(0.0,0.0)
          else
            dummy%val(i,j) = state%temp_f%val(i,j)
        end if
      end do
    end do

    dummy = rearrange_2Dspectrum((dummy))
    !dummy = deal_mask(dummy)
    !dummy = rearrange_2Dspectrum(deal_mask(dummy))

    !do i =0,xdim-1
    !  do j =0,ydim-1
    !    if((i>xdim/6 .AND.i<(5*xdim/6).AND.(j>ydim/6 .AND.j<(5*ydim/6)))) then
    !    dummy%val(i,j) =cmplx(0.0,0.0)
    !    end if
    !  end do
    !end do


		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_temp_f_remap!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,(abs(real(dummy%val(i,j),real_outp_precision))),&
                          (abs(real(aimag(dummy%val(i,j)),real_outp_precision)))
			end do
		end do
    close(20)
  end subroutine

  subroutine write_temp_f()
    !write fourier temp-field to file. analogous to write_chem
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/temp_f/'
		write(suffix,"(I5,A11)")  int(state%t/write_intervall), ".temp_f.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		!do i=0,xdim-1
	  !	do j=0,ydim-1
	  ! 		dummy%val(i,j) = dummy%val(i,j)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
    !  end do
    !end do

    dummy = state%temp_f
    dummy = rearrange_2Dspectrum(deal_mask(dummy))

		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_temp_f!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,log(abs(real(dummy%val(i,j),real_outp_precision))),&
                          log(abs(real(aimag(dummy%val(i,j)),real_outp_precision)))
			end do
		end do
    close(20)
  end subroutine
  subroutine write_u_stat()
    !write measured u- diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/u_stat/'
		write(suffix,"(A11)") "u_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_u_stat!'
	 	write(20,*)               state%step,                       & !1
                              state%t,                          & !2
                              measure_vmax(),                   & !3
                              measure_u_rms(),                  & !4
                              maxval(real(state%u%val(:,:,1))), & !5
                              minval(real(state%u%val(:,:,1))), & !6
                              maxval(real(state%u%val(:,:,2))), & !7
                              minval(real(state%u%val(:,:,2)))    !8
    close(20)
    end if
  end subroutine
  subroutine write_E_stat()
    !write measured E-diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/E_stat/'
		write(suffix,"(A11)") "E_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_E_stat!'
	 	  write(20,*) state%step,state%t,measure_Ekin(),measure_Epot(),measure_Ekin()+measure_Epot()
      close(20)
    end if
  end subroutine

  subroutine write_sys_stat()
    !write measured system wide diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: int_dummy_f
    type(sfield)                        :: int_dummy
    type(sfield)                        :: int1_dummy_f
    type(sfield)                        :: int1_dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=23),parameter					:: path ='./output/data/sys_stat/'
		write(suffix,"(A12)") "sys_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_sys_stat!'
      int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                                +state%iky%val(:,:)*state%u_f%val(:,:,2) 
      int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                               +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 
      call dfftw_execute_dft(ifull2D,int_dummy_f%val,int_dummy%val)
      call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

	 	  write(20,*) state%step,                                               & !1
                  state%t,                                                  & !2

                  0.0_rp,&
                  !maxval(real(int_dummy%val,real_outp_precision)),          & !3 maximum of divergence
                  shear,                                                    & !4 strength of shear
                  dt,                                                       & !5
                  maxval(real(int1_dummy%val,real_outp_precision)),         & !6 max brucker div 
                  measure_av(state%s_dummy_f%val)                             !7 average vorticity 

      close(20)
    end if
  end subroutine

  subroutine write_T_stat()
    !write measured T-diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/T_stat/'
		write(suffix,"(A10)") "T_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    !write(*,*) 'sub write_T_stat() has been called'
    if(state%step>=1) then
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_T_stat!'
	 	  write(20,*) state%step,state%t,maxval(real(state%temp%val,real_outp_precision)),measure_av(state%temp%val)&
                   ,minval(real(state%chem%val,real_outp_precision))
      close(20)
    end if
  end subroutine

  subroutine write_C_stat()
    !write measured T-diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/C_stat/'
		write(suffix,"(A11)") "C_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_C_stat!'
	 	  write(20,*) state%step,state%t,maxval(real(state%chem%val,real_outp_precision)),measure_av(state%chem%val)&
                   ,minval(real(state%chem%val,real_outp_precision))
      close(20)
    end if
  end subroutine

  function deal_mask(u_f)
  	! show the cutoff frequencies highlighted in output array. only for visualisation
  	type(sfield)				:: deal_mask
  	type(sfield)				:: u_f
  	complex(kind = rp)											:: max_abs
  	max_abs = cmplx(maxval(abs(real(u_f%val))),maxval(abs(real((aimag(u_f%val))))),rp)
  	deal_mask%val = u_f%val

  	deal_mask%val(xdim/3  ,:        ) 	= max_abs			
  	deal_mask%val(2*xdim/3,:        ) 	= max_abs			
  	deal_mask%val(:       ,2*ydim/3 ) 	= max_abs			
  	deal_mask%val(:       ,ydim/3   ) 	= max_abs
  end function 

  function rearrange_2Dspectrum(arr_f)
	  type(sfield)				:: arr_f
	  type(sfield)				:: rearrange_spec
	  type(sfield)				:: rearrange_2Dspectrum
  
	  do i =0,xdim/2
	  	do j =0,ydim/2
	  		rearrange_spec%val(i,j) = arr_f%val(i+xdim/2-1,j+ydim/2-1)  !lower left
	  	end do
	  end do
	  do i =xdim/2+1,xdim-1
	  	do j =ydim/2+1,ydim-1
	  		rearrange_spec%val(i,j) = arr_f%val(i-xdim/2-1,j-ydim/2-1) 	!upper right
	  	end do
	  end do
	  do i =0,xdim/2
	  	do j =ydim/2+1,ydim-1
	  		rearrange_spec%val(i,j) = arr_f%val(i+xdim/2-1,j-ydim/2-1) ! upper left
	  	end do
	  end do
	  do i =xdim/2+1,xdim-1
	  	do j =0,ydim/2
	  		rearrange_spec%val(i,j) = arr_f%val(i-xdim/2-1,j+ydim/2-1) ! lower right
	  	end do
	  end do
    rearrange_2Dspectrum = rearrange_spec
  end function

end module
