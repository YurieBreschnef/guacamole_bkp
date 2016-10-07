module init
  !module for all actions to be happening at program start
  use sys_state
	use const
  use timestepping
  !use test
	use plans
	use trafo 
  implicit none

  contains
  subroutine init_all()
    !call all special init subroutines  and initialize some variables
    call srand(seed)
    if(debuglevel .GE. 1) write(*,*) '-calling init_all()'
  	state%t	 	    = 0.00_rp
  	dt_2 	        = dt/2.0_rp
  	steps       	= int(tmax/dt,ip)
  	state%step    = 0
    T_rm          = real(Lx/(shear*Ly),rp) /2.0
    threads = omp_get_max_threads ( )     ! find out on how many threads this prog is running
    !threads = 1
    my_thread_id= omp_get_thread_num ( )
    write ( *, '(a,i8)' ) 'The number of processors available = ', omp_get_num_procs ( )
    write ( *, '(a,i8)' ) 'The number of threads available    = ', threads 

    call init_plans()   
    ! ALWAYS INITIALIZE PLANS BEFORE THE ARRAY, since the fftw planning routines overwrite
    ! the input arrays to estimate the best plan execution
    call init_u()
    call init_temp()
    call init_chem()
    call init_k()
    call init_plausibility()
    
    ! TODO write sub to accomodate old state initiation 
    if(debuglevel .GE. 1) write(*,*) '-done with init_all.'
  end subroutine

  subroutine init_plausibility
    if(debuglevel .GE. 1) write(*,*) '_________Plausibility___________________________________'
    if ((B_therm*S_therm-B_comp*s_comp)>=0) then
      if(debuglevel .GE. 1) write(*,*) '-Stable stratification with total density gradient of:',(B_therm*S_therm-B_comp*s_comp)
      else
      if(debuglevel .GE. 1) write(*,*) '-Unstable stratification with total density gradient of:',(B_therm*S_therm-B_comp*s_comp)
    end if
    !if(debuglevel .GE. 1) write(*,*) '-expected finger scale d approximately:',sqrt(sqrt(D_therm*D_visc/(B_therm*S_therm)))
    if(debuglevel .GE. 1) write(*,*) '-TEMP diffusive timescale tau=(dx**2)/diff :X ',((Lx/real(xdim))**2)/D_therm,'|dt:',dt
    if(((Lx/real(xdim))**2)/D_therm < dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
      !TODO make the plausi check the timestepping used, ETD is not as badly restricted as RK4 with diffusive timescale
    end if
    if(debuglevel .GE. 1) write(*,*) '-TEMP diffusive timescale tau=(dy**2)/diff :Y ',((Ly/real(ydim))**2)/D_therm,'|dt:',dt
    if(((Ly/real(ydim))**2)/D_therm < dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
      !TODO make the plausi check the timestepping used, ETD is not as badly restricted as RK4 with diffusive timescale
    end if
    if(debuglevel .GE. 1) write(*,*) '-CHEM diffusive timescale tau=(dx**2)/diff :X ',((Lx/real(xdim))**2)/D_comp,'|dt:',dt
    if(((Lx/real(xdim))**2)/D_comp< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
      !TODO make the plausi check the timestepping used, ETD is not as badly restricted as RK4 with diffusive timescale
    end if
    if(debuglevel .GE. 1) write(*,*) '-CHEM diffusive timescale tau=(dy**2)/diff :Y ',((Ly/real(ydim))**2)/D_comp,'|dt:',dt
    if(((Ly/real(ydim))**2)/D_comp< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
      !TODO make the plausi check the timestepping used, ETD is not as badly restricted as RK4 with diffusive timescale
    end if
    if(debuglevel .GE. 1) write(*,*) '- V   diffusive timescale tau=(dx**2)/diff :X ',((Lx/real(xdim))**2)/D_visc,'|dt:',dt
    if(((Lx/real(xdim))**2)/D_visc< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
      !TODO make the plausi check the timestepping used, ETD is not as badly restricted as RK4 with diffusive timescale
    end if
    if(debuglevel .GE. 1) write(*,*) '- V   diffusive timescale tau=(dy**2)/diff :Y ',((Ly/real(ydim))**2)/D_visc,'|dt:',dt
    if(((Ly/real(ydim))**2)/D_visc< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      stop
      !TODO make the plausi check the timestepping used, ETD is not as badly restricted as RK4 with diffusive timescale
    end if
    write(*,*) 'smallest possible number:', epsilon(1.0_rp)
    write(*,*) '1/eps:', 1.0/epsilon(1.0_rp)
    write(*,*) '1/eps**2:', 1.0/epsilon(1.0_rp)**2

    write(*,*)'_______________________________RUNTIME ESTIMATION______________________________'
    l = benchmarking
    benchmarking =1
    state_nm1 = state  
    call RK4_step()   ! to not compare the first writing step  with the second nonwriting
    call cpu_time(bm_timestepping_starttime)
    call RK4_step()
    call cpu_time(bm_timestepping_endtime)
    call bm_evaluate(.false.)     ! measure how long one step takes and make an estimation 
    write(*,*) 'RK4 step took:', bm_timestepping_time
    if(int(bm_step_time*real(steps))>=3600) then
      write(*,*)'ESTIMATED RUNTIME:', int(bm_timestepping_time*real(steps))/3600,'hours',mod(int(bm_timestepping_time*real(steps)),&
                            int(3600))/60,'min'
    else
      write(*,*)'ESTIMATED RUNTIME:', int(bm_timestepping_time*real(steps))/60,'min',&
            mod(int(bm_timestepping_time*real(steps)),int(60)),'sec'
    end if
    call cpu_time(bm_timestepping_starttime)
    call euler_step()
    call cpu_time(bm_timestepping_endtime)
    call bm_evaluate(.false.)     ! measure how long one step takes and make an estimation 
    write(*,*) '______________________________________________________'
    write(*,*) 'EULER step took:', bm_timestepping_time
    if(int(bm_step_time*real(steps))>=3600) then
      write(*,*)'ESTIMATED RUNTIME:', int(bm_timestepping_time*real(steps))/3600,'hours',mod(int(bm_timestepping_time*real(steps)),&
                            int(3600))/60,'min'
    else
      write(*,*)'ESTIMATED RUNTIME:', int(bm_timestepping_time*real(steps))/60,'min',&
            mod(int(bm_timestepping_time*real(steps)),int(60)),'sec'
    end if
    
    state = state_nm1       ! reset state to init

    benchmarking = l
    if(debuglevel .GE. 1) write(*,*) '_________END of Plausibility___ ________________________'
  end subroutine

  subroutine init_plans()

    !call fftw_init_threads()
    !call fftw_plan_with_nthreads(no_of_threads);
    if(debuglevel .GE. 1) write(*,*) '  -calling init_plans()'

  	call dfftw_plan_dft_1d( x_xf,xdim,x_pen	,x_pen_f,FFTW_FORWARD ,fftw_plan_thoroughness)
  	call dfftw_plan_dft_1d( y_yf,ydim,y_pen	,y_pen_f,FFTW_FORWARD ,fftw_plan_thoroughness)

  	call dfftw_plan_dft_1d( xf_x,xdim,x_pen_f	,x_pen,FFTW_BACKWARD,fftw_plan_thoroughness)
  	call dfftw_plan_dft_1d( yf_y,ydim,y_pen_f	,y_pen,FFTW_BACKWARD,fftw_plan_thoroughness)

  	call dfftw_plan_dft_2d(ifull2D,xdim,ydim,state%temp_f%val,state%temp%val,  &
      FFTW_BACKWARD,fftw_plan_thoroughness)
  	call dfftw_plan_dft_2d( full2D,xdim,ydim,state%temp%val,state%temp_f%val, &
      FFTW_FORWARD ,fftw_plan_thoroughness)

    if(debuglevel <= 1) write(*,*) '  -done with init_plans.'
  end subroutine
  
  subroutine init_u()
    real(kind=rp)                                   ::amp
    integer                                         ::xmodes,ymodes
    !initialize velocity field 
    if(debuglevel .GE. 1) write(*,*) '  -calling init_u()'
    state%u%val = cmplx(0.0_rp,0.0_rp,rp)
    xmodes = 2 !xdim/2
    ymodes = 2 !ydim/16 

    ! source in the middle of the field
    !do i=0,xdim-1
    !  do j=0,ydim-1
    !    state%u%val(i,j,1) = real(i-xdim/2)/abs(real(i-xdim/2))* exp(-(abs(real((i-xdim/2)**2+(j-ydim/2)**2)))) 
    !    state%u%val(i,j,2) = real(j-ydim/2)/abs(real(j-ydim/2))* exp(-(abs(real((i-xdim/2)**2+(j-ydim/2)**2)))) 
    !  end do
    !end do

    !do i=0,xdim-1
    !  do j=0,ydim-1
    !    state%u%val(i,j,1) = sin(real(xmodes) * (real(j)/real(ydim))*2.0_rp*pi)
    !    state%u%val(i,j,2) = sin(real(ymodes) * (real(i)/real(xdim))*2.0_rp*pi)
    !  end do
    !end do
    !do i=0,xdim-1
    !  do j=0,ydim-1
    !    amp = rand()
    !    state%u%val(i,j,1) = real((amp-0.5_rp),rp)
    !    amp = rand()
    !    state%u%val(i,j,2) = real((amp-0.5_rp),rp)
    !  end do
    !end do

    state%u%val = state%u%val *0.00001_rp                             
    call dfftw_execute_dft(full2D,state%u%val(:,:,1),state%u_f%val(:,:,1))
    call dfftw_execute_dft(full2D,state%u%val(:,:,2),state%u_f%val(:,:,2))
    state%u_f%val(0,0,:) = cmplx(0.0_rp,0.0_rp)        ! set constant mode to zero 
    state%u_f%val = state%u_f%val/real(xdim*ydim,rp)   !FFTW NORM
    if(debuglevel .GE. 1) write(*,*) '  -done with init_u.'
  end subroutine

  subroutine init_temp()
    integer                                         ::xpos,ypos
    integer                                         ::xpoints,ypoints
    real(kind=rp)                                   ::amp
    if(debuglevel .GE.1) write(*,*) '  -calling init_temp()'
    !initialize temp field 
    state%temp%val = cmplx(0.0_rp,0.0_rp,rp)
    state%temp_f%val = cmplx(0.0_rp,0.0_rp,rp)

    xpoints = 2 
    ypoints = 2
    do xpos=xdim/xpoints,(xpoints-1)*xdim/xpoints,xdim/xpoints
      do ypos=ydim/ypoints,(ypoints-1)*ydim/ypoints,ydim/ypoints
      amp = rand()
        do i=0,xdim-1
          do j=0,ydim-1
              ! FOURIER INIT
              !state%temp%val(i,j) = state%temp%val(i,j) &
              !+cmplx((amp-0.5_rp)*exp(-( (20.0_rp*real(j-ypos,rp)/real(ydim,rp))**2 &
              !             +(20.0_rp*real(i-xpos,rp)/real(xdim,rp))**2) ),0.0_rp,rp)
              !state%temp_f%val(i,j) = cmplx(1.0_rp,1.0_rp,rp) 
              !-----------------------------------------------------------------------

              state%temp%val(i,j) = state%temp%val(i,j) &
              +cmplx((amp-0.5_rp)*exp(-((30.0_rp*real(j-ypos,rp)/real(ydim,rp))**2 &
                           +(30.0_rp*real(i-xpos,rp)/real(xdim,rp))**2) ),0.0_rp,rp)

              !amp = (rand()-0.5_rp)
              !state%temp%val(i,j) = amp
              !state%temp%val(i,j) =  state%temp%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi*(xdim)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi*(xdim/2)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi*(xdim/4)) 

              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi*(xdim)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi*(xdim/2)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi*(xdim/4)) 

              !state%temp%val(i,j) =  state%temp%val(i,j)+  sin(real(j)/real(ydim)*2.0_rp*pi) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  sin(real(j)/real(ydim)*2.0_rp*pi*(ydim/2)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  sin(real(j)/real(ydim)*2.0_rp*pi*(ydim/4)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(j)/real(ydim)*2.0_rp*pi*(ydim)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(j)/real(ydim)*2.0_rp*pi*(ydim/2)) 
              !state%temp%val(i,j) =  state%temp%val(i,j)+  cos(real(j)/real(ydim)*2.0_rp*pi*(ydim/4)) 
          end do
        end do
      end do
    end do
    state%temp%val = state%temp%val*0.0010_rp

    state%temp_f%val(0,0) = cmplx(0.0_rp,0.0_rp,rp)
    call dfftw_execute_dft(full2D,state%temp%val(:,:),state%temp_f%val(:,:))
    state%temp_f%val = state%temp_f%val/real(xdim*ydim,rp)   !FFTW NORM

    !call dfftw_execute_dft(ifull2D,state%temp_f%val(:,:),state%temp%val(:,:))

    if(debuglevel .GE. 1) write(*,*) '  -done with init_temp.'
  end subroutine

  subroutine init_chem()
    integer                                         ::xpos,ypos
    integer                                         ::xpoints,ypoints
    real(kind=rp)                                   ::amp
    if(debuglevel .GE. 1) write(*,*) '  -calling init_chem()'
    !initialize chemical field 
    state%chem%val = cmplx(0.0_rp,0.0_rp,rp)

    xpoints = 3 
    ypoints = 3
    do xpos=xdim/xpoints,(xpoints-1)*xdim/xpoints,xdim/xpoints
      do ypos=ydim/ypoints,(ypoints-1)*ydim/ypoints,ydim/ypoints
      amp = (rand())
        do i=0,xdim-1
          do j=0,ydim-1
              state%chem%val(i,j) = state%chem%val(i,j) &
              +cmplx((amp-0.5_rp)*exp(-( (30.0_rp*real(j-ypos,rp)/real(ydim,rp))**2 &
                           +(30.0_rp*real(i-xpos,rp)/real(xdim,rp))**2) ),0.0_rp,rp)

!              amp = (rand()-0.5_rp)
!              state%chem%val(i,j) = amp
              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi*(xdim)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi*(xdim/2)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(i)/real(xdim)*2.0_rp*pi*(xdim/4)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi*(xdim)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi*(xdim/2)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  cos(real(i)/real(xdim)*2.0_rp*pi*(xdim/4)) 

              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(j)/real(ydim)*2.0_rp*pi) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(j)/real(ydim)*2.0_rp*pi*(ydim/2)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  sin(real(j)/real(ydim)*2.0_rp*pi*(ydim/4)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  cos(real(j)/real(ydim)*2.0_rp*pi*(ydim)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  cos(real(j)/real(ydim)*2.0_rp*pi*(ydim/2)) 
              !state%chem%val(i,j) =  state%chem%val(i,j)+  cos(real(j)/real(ydim)*2.0_rp*pi*(ydim/4)) 

          end do
        end do
      end do
    end do
    state%chem%val = state%chem%val * 0.0010_rp

    call dfftw_execute_dft(full2D,state%chem%val(:,:),state%chem_f%val(:,:))
    state%chem_f%val = state%chem_f%val/real(xdim*ydim,rp)   !FFTW NORM
    state%temp_f%val(0,0) = cmplx(0.0_rp,0.0_rp,rp)
    if(debuglevel .GE. 1) write(*,*) '  -done with init_chem.'
  end subroutine

  subroutine init_k()
    !make k-values for derivatives
    if(debuglevel .GE. 1) write(*,*) '  -calling init_k()'
    do i=0,xdim-1
      do j=0,ydim-1
          if(i<=xdim/2)then   ! kx
            state%ikx%val(i,:) = cmplx(0.0_rp,(pi*2.0_rp*real(i,rp)/Lx),rp)
          else
            state%ikx%val(i,:) = cmplx(0.0_rp,(pi*2.0_rp*real(i-xdim,rp)/Lx),rp)
          end if

          if(j<=ydim/2)then   !ky
            state%iky%val(:,j) = cmplx(0.0_rp,(pi*2.0_rp*real(j,rp)/Ly),rp)
          else
            state%iky%val(:,j) = cmplx(0.0_rp,(pi*2.0_rp*real(j-ydim,rp)/Ly),rp)
          end if

      end do
    end do
     IF(ALL(state%ikx%val ==0.0_rp).OR.ALL(state%iky%val ==0.0_rp))  then
       write(*,*) 'sub init_k(): ALL ikx or iky are ZERO. BAD. VERY BAD.'
       stop
     end if


    state%ikx_sqr%val = state%ikx%val**2
    state%ikx_sqr%val(0,:) = epsilon(1.0_rp)
    state%iky_sqr%val = state%iky%val**2
    state%iky_sqr%val(:,0) = epsilon(1.0_rp)
    state%iki_sqr%val = state%ikx_sqr%val + state%iky_sqr%val 
    state%iki_sqr%val(0,0) = epsilon(1.0_rp)

    ! also set iki_bar variables because if shearing is of they need to be set at least once
    state%ikx_bar%val(:,:) = state%ikx%val(:,:) 
    state%iky_bar%val(:,:) = state%iky%val(:,:)
    state%ikx_bar_sqr%val(:,:) = state%ikx_bar%val(:,:)**2
    state%iky_bar_sqr%val(:,:) = state%iky_bar%val(:,:)**2
    state%iki_bar_sqr%val(:,:) = state%ikx_bar%val(:,:)**2 + state%iky_bar%val(:,:)**2

    ! used for crossprodukt in fu_Nuk in module pdgl
    state%k_vec%val(:,:,1) = real(imag*state%ikx_bar%val(:,:),rp)
    state%k_vec%val(:,:,2) = real(imag*state%iky_bar%val(:,:),rp)

    write(*,*) 'epsilon:',epsilon(1.0_rp)
    write(*,*) '1/epsilon:',1.0_rp/epsilon(1.0_rp)
    write(*,*) '1/epsilon**2:',1.0_rp/(epsilon(1.0_rp)**2)

    IF(ANY(IsNaN(real(state%ikx%val))))  then
      write(*,*) 'init_k(): NAN detected! in real part of ikx'
      stop
    end if
    IF(ANY(IsNaN(real(state%iky%val))))  then
      write(*,*) 'init_k(): NAN detected! in real part of iky'
      stop
    end if
    IF(ANY(IsNaN(AIMAG(state%ikx%val))))  then
      write(*,*) 'init_k(): NAN detected! in imag part of ikx'
      stop
    end if
    IF(ANY(IsNaN(AIMAG(state%iky%val))))  then
      write(*,*) 'init_k(): NAN detected! in imag part of iky'
      stop
    end if
    call set_ik_bar(sheartime)  
    if(debuglevel .GE. 1) write(*,*) '  -done with init_k.'
  end subroutine
end module
