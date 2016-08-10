module timestepping
  !implements time-stepping scemes that are called by e.g. 'call RK4_step()'
  use sys_state
  use pdgl
  use plans
  implicit none

  contains

  subroutine reset_dt(new_dt)
    real(kind=rp)         ::new_dt
    dt      = new_dt
    dt_2    = dt*(1.0_rp/2.0_rp)           
    dt_3    = dt*(1.0_rp/3.0_rp)       
    dt_4    = dt*(1.0_rp/4.0_rp)       
    dt_8    = dt*(1.0_rp/8.0_rp)       
    dt_34   = dt*(3.0_rp/4.0_rp)       
    dt_29   = dt*(2.0_rp/9.0_rp)       
    dt_49   = dt*(4.0_rp/9.0_rp)       
    dt_724  = dt*(7.0_rp/24.0_rp)         
  end subroutine

  subroutine RK4_adjust_dt()
  !primitive timestep adjustment
  !TODO make proper adaptive RK-timestepping
    dt=0.05_rp/measure_vmax()
    if(dt > dt_max) then
        call reset_dt(dt_max)
    end if
    if(dt < dt_min) then
        call reset_dt(dt_min)
    end if
  end subroutine


  subroutine RK4_step()
  	!performs a timestep with RK4 and stores the new result in u_f,temp_f,chem_f
    if(debuglevel .GE.3) write(*,*)'RK4 sub called'
    call dealiase_all()
    !_____________________k1_________________________________
  	state%u_k1%val = fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
  	state%t_k1%val = ft(state%u_f%val ,state%temp_f%val ,state%t)     
  	state%c_k1%val = fc(state%u_f%val ,state%chem_f%val ,state%t)     
    !_____________________k2_________________________________
  	state%u_k2%val = fu(state%u_f%val   +dt_2*state%u_k1%val,&  !f(u_f,temp_f,chem_f,t)
                        state%temp_f%val+dt_2*state%t_k1%val,&
                        state%chem_f%val+dt_2*state%c_k1%val,&
                        state%t+dt_2)     

  	state%t_k2%val = ft(state%u_f%val   +dt_2*state%u_k1%val,&  !ft(u_f,temp_f,t)
                        state%temp_f%val+dt_2*state%t_k1%val,&
                        state%t+dt_2)     
  	state%c_k2%val = fc(state%u_f%val   +dt_2*state%u_k1%val,&  !ft(u_f,chem_f,t)
                        state%chem_f%val+dt_2*state%c_k1%val,&
                        state%t+dt_2)     
    !_____________________k3_________________________________
  	state%u_k3%val = fu(state%u_f%val   +dt_2*state%u_k2%val,&
                        state%temp_f%val+dt_2*state%t_k2%val,&
                        state%chem_f%val+dt_2*state%c_k2%val,&
                        state%t+dt_2)     
  	state%t_k3%val = ft(state%u_f%val   +dt_2*state%u_k2%val,&
                        state%temp_f%val+dt_2*state%t_k2%val,&
                        state%t+dt_2)     
  	state%c_k3%val = fc(state%u_f%val   +dt_2*state%u_k2%val,&  
                        state%chem_f%val+dt_2*state%c_k2%val,&
                        state%t+dt_2)     
    !_____________________k4_________________________________
  	state%u_k4%val = fu(state%u_f%val   +dt*state%u_k3%val,&
                        state%temp_f%val+dt*state%t_k3%val,&
                        state%chem_f%val+dt*state%c_k3%val,&
                        state%t+dt_2)     
  	state%t_k4%val = ft(state%u_f%val   +dt*state%u_k3%val,&
                        state%temp_f%val+dt*state%t_k3%val,&
                        state%t+dt_2)     

  	state%c_k4%val = fc(state%u_f%val   +dt*state%u_k3%val,&  
                        state%chem_f%val+dt*state%c_k3%val,&
                        state%t+dt_2)     
    !____________________step______________________________'_
  	state%u_f%val     =state%u_f%val     +(dt/6.0_rp)*(      state%u_k1%val&
                                                     +2.0_rp*state%u_k2%val&
                                                     +2.0_rp*state%u_k3%val&
                                                            +state%u_k4%val   )

  	state%temp_f%val  =state%temp_f%val  +(dt/6.0_rp)*(       state%t_k1%val&
                                                      +2.0_rp*state%t_k2%val&
                                                      +2.0_rp*state%t_k3%val&
                                                             +state%t_k4%val)

  	state%chem_f%val  =state%chem_f%val  +(dt/6.0_rp)*(       state%c_k1%val&
                                                      +2.0_rp*state%c_k2%val&
                                                      +2.0_rp*state%c_k3%val&
                                                             +state%c_k4%val)
  	state%t=state%t+dt
  	state%step=state%step+1
  end subroutine

  subroutine euler_step()
  	!performs a timestep with simple euler and stores the new result in u_f,temp_f,chem_f
    if(debuglevel .GE.3) write(*,*)'RK4 sub called'
    call dealiase_all()
  	state%u_f%val    =state%u_f%val    + dt*fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
  	state%temp_f%val =state%temp_f%val + dt*ft(state%u_f%val ,state%temp_f%val ,state%t)     
  	state%chem_f%val =state%chem_f%val + dt*fc(state%u_f%val ,state%chem_f%val ,state%t)     
  	state%t=state%t+dt
  	state%step=state%step+1
  end subroutine

recursive subroutine BogSham_step()
  !implements an adaptive Runge-Kutta type timestepping (after Bogacki-Shampine)
  if(debuglevel .GE.3) write(*,*)'BogSham sub called with dt=',dt
  call dealiase_all()
  !_____________________k1_________________________________
  state%u_k1%val = fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
  state%t_k1%val = ft(state%u_f%val ,state%temp_f%val ,state%t)     
  state%c_k1%val = fc(state%u_f%val ,state%chem_f%val ,state%t)     
  !_____________________k2_________________________________
  state%u_k2%val = fu(state%u_f%val   +dt_2*state%u_k1%val,&
                      state%temp_f%val+dt_2*state%t_k1%val,&
                      state%chem_f%val+dt_2*state%c_k1%val,&
                      state%t+dt_2)     
  state%t_k2%val = ft(state%u_f%val   +dt_2*state%u_k1%val,&
                      state%temp_f%val+dt_2*state%t_k1%val,&
                      state%t+dt_2)     
  state%c_k2%val = fc(state%u_f%val   +dt_2*state%u_k1%val,&  
                      state%chem_f%val+dt_2*state%c_k1%val,&
                      state%t+dt_2)     
  !_____________________k3_________________________________
  state%u_k3%val = fu(state%u_f%val   +dt_34*state%u_k2%val,&  !f(u_f,temp_f,chem_f,t)
                      state%temp_f%val+dt_34*state%t_k2%val,&
                      state%chem_f%val+dt_34*state%c_k2%val,&
                      state%t+dt_34)     
 
  state%t_k3%val = ft(state%u_f%val   +dt_34*state%u_k2%val,&  !ft(u_f,temp_f,t)
                      state%temp_f%val+dt_34*state%t_k2%val,&
                      state%t+dt_34)     
  state%c_k3%val = fc(state%u_f%val   +dt_34*state%u_k2%val,&  !ft(u_f,chem_f,t)
                      state%chem_f%val+dt_34*state%c_k2%val,&
                      state%t+dt_34)     
  !calculate second order step, the higher order step is used for
  !error estimation later on 
  state%dummy_f%val = state%u_f%val         +dt_29*state%u_k1%val&
                                            +dt_3 *state%u_k2%val&
                                            +dt_49*state%u_k3%val
  state%t_dummy_f%val = state%temp_f%val   +dt_29*state%t_k1%val&
                                           +dt_3 *state%t_k2%val&
                                           +dt_49*state%t_k3%val
  state%c_dummy_f%val = state%chem_f%val   +dt_29*state%c_k1%val&
                                           +dt_3 *state%c_k2%val&
                                           +dt_49*state%c_k3%val
 !_____________________k4_________________________________
  state%u_k4%val = fu(state%dummy_f%val,&
                      state%t_dummy_f%val,&
                      state%c_dummy_f%val,&
                      state%t+dt)     
  state%t_k4%val = ft(state%dummy_f%val,&
                      state%t_dummy_f%val,&
                      state%t+dt)     
 
  state%c_k4%val = fc(state%dummy_f%val,&  
                      state%c_dummy_f%val,&
                      state%t+dt)     
  !____________________step______________________________'_
  state%z_dummy_f%val     =state%dummy_f%val +dt_724*(     state%u_k1%val&
                                             +dt_4  *      state%u_k2%val&
                                             +dt_3  *      state%u_k3%val&
                                             +dt_8  *      state%u_k4%val   )
 
  state%tz_dummy_f%val  =state%t_dummy_f%val +dt_724*(     state%t_k1%val&
                                             +dt_4  *      state%t_k2%val&
                                             +dt_3  *      state%t_k3%val&
                                             +dt_8  *      state%t_k4%val   )

  state%cz_dummy_f%val  =state%c_dummy_f%val +dt_724*(     state%c_k1%val&
                                             +dt_4  *      state%c_k2%val&
                                             +dt_3  *      state%c_k3%val&
                                             +dt_8  *      state%c_k4%val   )


  if(debuglevel .GE.3) write(*,*)'BogSham sub testing wether error is too big..'
  if(maxval(real(state%dummy_f%val-state%z_dummy_f%val,rp))>max_step_error) then
    if(debuglevel .GE.3) write(*,*)'BogSham sub error was too big! reset time and recurse! ERROR:',&
      maxval(real(state%dummy_f%val-state%z_dummy_f%val,rp)),'dt:',dt

      if(dt*(1.0_rp-stepwidth_adjustment_rate) >dt_min) then
        ! if stepsize can still be reduced then do so and call recursively
         call reset_dt(dt*(1.0_rp - stepwidth_adjustment_rate))
         call BogSham_step()
      else
        ! if stepsize can not be reduced even further then move on anyway
        state%u_f%val = state%dummy_f%val
        state%temp_f%val = state%t_dummy_f%val
        state%chem_f%val = state%c_dummy_f%val
        state%t=state%t+dt
        state%step=state%step+1
      end if

  else
    if(debuglevel .GE.3) write(*,*)'BogSham sub error was ok increase stepwidth and move on!'
    state%u_f%val = state%dummy_f%val
    state%temp_f%val = state%t_dummy_f%val
    state%chem_f%val = state%c_dummy_f%val
    state%t=state%t+dt
    state%step=state%step+1

    if(dt*(1.0_rp+stepwidth_adjustment_rate) <dt_max) then
      ! check wether the stepwidth is already too high
      call reset_dt(dt*(1.0_rp + stepwidth_adjustment_rate))
    end if
  end if
end subroutine

end module
