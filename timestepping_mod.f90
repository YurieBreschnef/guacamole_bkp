module timestepping
  !implements time-stepping scemes that are called by e.g. 'call RK4_step()','call euler_step()'
  use sys_state
  use pdgl
  use plans
  use remap
  implicit none

  contains
subroutine reset_dt(new_dt)
  ! resets all dt-dependencies to new_dt
  real(kind = rp)               :: new_dt
  dt      = new_dt
  dt_2    = dt*(1.0_rp/ 2.0_rp)           
  dt_3    = dt*(1.0_rp/ 3.0_rp)       
  dt_4    = dt*(1.0_rp/ 4.0_rp)       
  dt_8    = dt*(1.0_rp/ 8.0_rp)       
  dt_34   = dt*(3.0_rp/ 4.0_rp)       
  dt_29   = dt*(2.0_rp/ 9.0_rp)       
  dt_49   = dt*(4.0_rp/ 9.0_rp)       
  dt_724  = dt*(7.0_rp/24.0_rp)         
end subroutine

subroutine RK4_adjust_dt()
!primitive timestep adjustment
  dt=0.05_rp/measure_vmax()
  if(dt > dt_max) then
      call reset_dt(dt_max)
  end if
  if(dt < dt_min) then
      call reset_dt(dt_min)
  end if
end subroutine
!------------------------------------------------------------------------------------------
subroutine RK4_step()
	!performs a timestep with RK4 and stores the new result in u_f,temp_f,chem_f


  if(debuglevel .GE.3) write(*,*)'RK4 sub called'
  !_____________________k1_________________________________
  call set_ik_bar(sheartime) 
	state%u_k1%val = fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
	state%t_k1%val = ft(state%u_f%val ,state%temp_f%val ,state%t)     
	state%c_k1%val = fc(state%u_f%val ,state%chem_f%val ,state%t)     
  !_____________________k2_________________________________
  call set_ik_bar(sheartime) 
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
  call set_ik_bar(sheartime) 
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
  call set_ik_bar(sheartime) 
	state%u_k4%val = fu(state%u_f%val   +dt*state%u_k3%val,&
                      state%temp_f%val+dt*state%t_k3%val,&
                      state%chem_f%val+dt*state%c_k3%val,&
                      state%t+dt)     
	state%t_k4%val = ft(state%u_f%val   +dt*state%u_k3%val,&
                      state%temp_f%val+dt*state%t_k3%val,&
                      state%t+dt)     

	state%c_k4%val = fc(state%u_f%val   +dt*state%u_k3%val,&  
                      state%chem_f%val+dt*state%c_k3%val,&
                      state%t+dt)     
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
	sheartime = sheartime+dt
	state%t=state%t+dt

  call dealiase_all()
  ! REMAPPING
  if(remapping==1 .AND.shearing==1.) then
    call remap_stepwise()
  end if
  call dealiase_all()

	state%step=state%step+1
end subroutine
!------------------------------------------------------------------------------------------
subroutine euler_step()
	!performs a timestep with simple euler and stores the new result in u_f,temp_f,chem_f
  if(debuglevel .GE.3) write(*,*)'euler sub called'

  ! REMAPPING
  if(remapping==1 .AND.shearing==1.) then
    call remap_stepwise()
  end if

  call set_ik_bar(sheartime) 
	state_np1%u_f%val    = state%u_f%val    + dt*fu(state%u_f%val ,state%temp_f%val,state%chem_f%val,sheartime)     
	state_np1%temp_f%val = state%temp_f%val + dt*ft(state%u_f%val ,state%temp_f%val ,sheartime)     
	state_np1%chem_f%val = state%chem_f%val + dt*fc(state%u_f%val ,state%chem_f%val ,sheartime)     

  state%u_f%val    = state_np1%u_f%val
  state%temp_f%val = state_np1%temp_f%val
  state%chem_f%val = state_np1%chem_f%val
  call dealiase_all()

	sheartime = sheartime+dt
	state%t   = state%t+dt
	state%step= state%step+1
end subroutine
!------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine ETD2_step()
!perfomrs ETD2 timestep as in brucker2007
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_RHS_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_exp_qh ! read as exp(q*h)

  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_RHS_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_exp_qh ! read as exp(q*h)

  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_RHS_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_exp_qh ! read as exp(q*h)


  call set_ik_bar(state%t)
  ! set q-values for exponent, note the minus sign in iki_sqr
  u_q(:,:,1)  = D_visc  *state%iki_bar_sqr%val(:,:)
  u_q(:,:,2)  = D_visc  *state%iki_bar_sqr%val(:,:)
  t_q         = D_therm *state%iki_bar_sqr%val
  c_q         = D_comp  *state%iki_bar_sqr%val

  u_exp_qh = exp(u_q*dt)
  t_exp_qh = exp(t_q*dt)
  c_exp_qh = exp(c_q*dt)

  ! calc RHS_n
  u_RHS_n   =  fu_N(state%u_f%val*u_exp_qh,state%temp_f%val*t_exp_qh,state%chem_f%val*c_exp_qh,state%t)
  t_RHS_n   =  ft_N(state%u_f%val*u_exp_qh,state%temp_f%val*t_exp_qh                          ,state%t)
  c_RHS_n   =  fc_N(state%u_f%val*u_exp_qh,                          state%chem_f%val*c_exp_qh,state%t)

  call set_ik_bar(state%t+dt)
  ! set q-values for exponent, note the minus sign in iki_sqr
  u_q(:,:,1)  = -D_visc  *state%iki_bar_sqr%val(:,:)
  u_q(:,:,2)  = -D_visc  *state%iki_bar_sqr%val(:,:)
  t_q         = -D_therm *state%iki_bar_sqr%val
  c_q         = -D_comp  *state%iki_bar_sqr%val
  u_exp_qh = exp(u_q*dt)
  t_exp_qh = exp(t_q*dt)
  c_exp_qh = exp(c_q*dt)
  ! calc RHS_n+1
  u_RHS_np1 =  fu_N(state%u_f%val*u_exp_qh+dt*u_RHS_n,&
                    state%temp_f%val*t_exp_qh+dt*t_RHS_n,&
                    state%chem_f%val*c_exp_qh+dt*c_RHS_n,state%t+dt)
  t_RHS_np1 =  ft_N(state%u_f%val*u_exp_qh+dt*u_RHS_n,&
                    state%temp_f%val*t_exp_qh+dt*t_RHS_n,state%t+dt)
  c_RHS_np1 =  fc_N(state%u_f%val*u_exp_qh+dt*u_RHS_n,&
                    state%chem_f%val*c_exp_qh+dt*c_RHS_n,state%t+dt)

  ! timestep
  state%u_f%val     = -state%u_f%val*u_exp_qh + (dt/2.0_rp)&
                      *(-u_RHS_n*u_exp_qh + u_RHS_np1)
  state%temp_f%val  = -state%temp_f%val*t_exp_qh + (dt/2.0_rp)&
                      *(-t_RHS_n*t_exp_qh + t_RHS_np1)
  state%chem_f%val  = -state%chem_f%val*c_exp_qh + (dt/2.0_rp)&
                      *(-c_RHS_n*c_exp_qh + c_RHS_np1)

  call dealiase_all()

	sheartime = sheartime+dt
	state%t   = state%t+dt
	state%step= state%step+1
end subroutine

end module
