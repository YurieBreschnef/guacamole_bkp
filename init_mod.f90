module init
  !module for all actions to be happening at program start
  use sys_state
	use const
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

    call init_plans()   
    ! ALWAYS INITIALIZE PLANS BEFORE THE ARRAY, since the fftw planning routines overwrite
    ! the input arrays to estimate the best plan execution
    call init_u()
    call init_temp()
    call init_chem()
    call init_k()
    call init_plausibility()
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
    end if
    if(debuglevel .GE. 1) write(*,*) '-TEMP diffusive timescale tau=(dy**2)/diff :Y ',((Ly/real(ydim))**2)/D_therm,'|dt:',dt
    if(((Ly/real(ydim))**2)/D_therm < dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
    end if
    if(debuglevel .GE. 1) write(*,*) '-CHEM diffusive timescale tau=(dx**2)/diff :X ',((Lx/real(xdim))**2)/D_comp,'|dt:',dt
    if(((Lx/real(xdim))**2)/D_comp< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
    end if
    if(debuglevel .GE. 1) write(*,*) '-CHEM diffusive timescale tau=(dy**2)/diff :Y ',((Ly/real(ydim))**2)/D_comp,'|dt:',dt
    if(((Ly/real(ydim))**2)/D_comp< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
    end if
    if(debuglevel .GE. 1) write(*,*) '- V   diffusive timescale tau=(dx**2)/diff :X ',((Lx/real(xdim))**2)/D_visc,'|dt:',dt
    if(((Lx/real(xdim))**2)/D_visc< dt) then
      write(*,*) ''
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      write(*,*) ''
      stop
    end if
    if(debuglevel .GE. 1) write(*,*) '- V   diffusive timescale tau=(dy**2)/diff :Y ',((Ly/real(ydim))**2)/D_visc,'|dt:',dt
    if(((Ly/real(ydim))**2)/D_visc< dt) then
      write(*,*) 'WARNING: Diffusive timescale is shorter than dt!'
      write(*,*) '_______________________________________________________________________________________________________________'
      stop
    end if
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
    !initialize velocity field 
    if(debuglevel .GE. 1) write(*,*) '  -calling init_u()'
    state%u%val = cmplx(0.0_rp,0.0_rp,rp)
    do i=0,xdim-1
      do j=0,ydim-1

      end do
    end do
    state%u%val = state%u%val /10.0_rp                             
    call dfftw_execute_dft(full2D,state%u%val(:,:,1),state%u_f%val(:,:,1))
    call dfftw_execute_dft(full2D,state%u%val(:,:,2),state%u_f%val(:,:,2))
    state%u_f%val = state%u_f%val/real(xdim*ydim,rp)   !FFTW NORM
    if(debuglevel .GE. 1) write(*,*) '  -done with init_u.'
  end subroutine

  subroutine init_temp()
    integer                                         ::xpos,ypos
    real(kind=rp)                                   ::amp
    if(debuglevel .GE.1) write(*,*) '  -calling init_temp()'
    !initialize temp field 
    state%temp%val = cmplx(0.0_rp,0.0_rp,rp)
    state%temp%val = cmplx(0.0_rp,0.0_rp,rp)

    do xpos=xdim/10,9*xdim/10,xdim/10
      do ypos=ydim/16,15*ydim/16,ydim/16
      amp = rand()
        do i=0,xdim-1
          do j=0,ydim-1
              state%temp%val(i,j) = state%temp%val(i,j) &
              +cmplx((amp-0.5_rp)*exp(-( (40.0_rp*real(j-ypos,rp)/real(ydim,rp))**2 &
                           +(40.0_rp*real(i-xpos,rp)/real(xdim,rp))**2) ),0.0_rp,rp)

              

              !amp = rand()
              !state%temp%val(i,j) = amp
          end do
        end do
      end do
    end do
    state%temp%val = state%temp%val*0.10_rp

    !call s_trafo(state%temp,state%temp_f,0)
    call dfftw_execute_dft(full2D,state%temp%val(:,:),state%temp_f%val(:,:))
    state%temp_f%val = state%temp_f%val/real(xdim*ydim,rp)   !FFTW NORM
    if(debuglevel .GE. 1) write(*,*) '  -done with init_temp.'
  end subroutine

  subroutine init_chem()
    integer                                         ::xpos,ypos
    real(kind=rp)                                   ::amp
    if(debuglevel .GE. 1) write(*,*) '  -calling init_chem()'
    !initialize chemical field 
    state%chem%val = cmplx(0.0_rp,0.0_rp,rp)

    !do xpos=xdim/16,xdim,xdim/16
    !  do ypos=ydim/7,ydim,ydim/7
    !  amp = rand()
        do i=0,xdim-1
          do j=0,ydim-1
    !          state%chem%val(i,j) = state%chem%val(i,j) &
    !          +cmplx((amp-0.5_rp)*exp(-( (20.0_rp*real(j-ypos,rp)/real(ydim,rp))**2 &
    !                       +(30.0_rp*real(i-xpos,rp)/real(xdim,rp))**2) ),0.0_rp,rp)
              amp = rand()
              state%temp%val(i,j) = amp
          end do
        end do
    !  end do
    !end do
    !state%temp%val = state%temp%val/1.0_rp

    state%chem%val = state%chem%val / 10.0
    call dfftw_execute_dft(full2D,state%chem%val(:,:),state%chem_f%val(:,:))
    state%chem_f%val = state%chem_f%val/real(xdim*ydim,rp)   !FFTW NORM
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

    state%ikx_sqr%val = state%ikx%val**2
    state%ikx_sqr%val(0,0) = epsilon(1.0_rp)
    state%iky_sqr%val = state%iky%val**2
    state%iky_sqr%val(0,0) = epsilon(1.0_rp)
    state%iki_sqr%val = state%ikx_sqr%val + state%iky_sqr%val 
    state%iki_sqr%val(0,0) = epsilon(1.0_rp)

    call set_ik_bar(state%t)  
    if(debuglevel .GE. 1) write(*,*) '  -done with init_k.'
  end subroutine
end module
