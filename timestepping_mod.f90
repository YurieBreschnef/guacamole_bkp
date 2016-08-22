module timestepping
  !implements time-stepping scemes that are called by e.g. 'call RK4_step()','call euler_step()'
  use sys_state
  use pdgl
  use plans
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
  call dealiase_all()
  !_____________________k1_________________________________
  
  call set_ik_bar(state%t) 
	state%u_k1%val = fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
	state%t_k1%val = ft(state%u_f%val ,state%temp_f%val ,state%t)     
	state%c_k1%val = fc(state%u_f%val ,state%chem_f%val ,state%t)     
  !_____________________k2_________________________________
  call set_ik_bar(state%t+dt_2) 
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
  call set_ik_bar(state%t+dt_2) 
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
  call set_ik_bar(state%t+dt) 
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
	state%t=state%t+dt
	state%step=state%step+1
end subroutine
!------------------------------------------------------------------------------------------
subroutine euler_step()
	!performs a timestep with simple euler and stores the new result in u_f,temp_f,chem_f
  if(debuglevel .GE.3) write(*,*)'RK4 sub called'
  call dealiase_all()
  call set_ik_bar(state%t) 
	state%u_f%val    =state%u_f%val    + dt*fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
	state%temp_f%val =state%temp_f%val + dt*ft(state%u_f%val ,state%temp_f%val ,state%t)     
	state%chem_f%val =state%chem_f%val + dt*fc(state%u_f%val ,state%chem_f%val ,state%t)     
	state%t=state%t+dt
	state%step=state%step+1
end subroutine
!------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine ETD2_step()
!perfomrs ETD2 timestep
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_NL_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_NL_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_NL_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_NL_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_NL_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_NL_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_exp_qh
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_exp_qh
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_exp_qh
  real(kind = rp)                                              :: q,qsqr
  real(kind = rp)                                              :: qh

  !TODO ineffective to recalculate
  !write(*,*) 'ETD: calculating NL-parts..'
  call set_ik_bar(state%t)

  fu_exp_qh(:,:,1) = cmplx(exp(-D_visc* real(state%iki_bar_sqr%val(:,:),rp)*dt),0.0_rp,rp)
  fu_exp_qh(:,:,2) = cmplx(exp(-D_visc* real(state%iki_bar_sqr%val(:,:),rp)*dt),0.0_rp,rp)
  ft_exp_qh(:,:)   = cmplx(exp(-D_therm*real(state%iki_bar_sqr%val(:,:),rp)*dt),0.0_rp,rp)
  fc_exp_qh(:,:)   = cmplx(exp(-D_comp* real(state%iki_bar_sqr%val(:,:),rp)*dt),0.0_rp,rp)
  !NOTE the minus sign hidden in iki_bar

  IF(ALL((real(fu_exp_qh,rp).EQ.(0.0_rp))))           write(*,*) 'WARNING: fu_exp_qh is EQ epsilon'
  IF(ALL((real(fu_exp_qh,rp).LE.epsilon(1.0_rp))))    write(*,*) 'WARNING: fu_exp_qh is LE epsilon'
  IF(ANY((real(fu_exp_qh,rp).LE.epsilon(1.0_rp))))    write(*,*) 'WARNING: fu_exp_qh has LE epsilon'


  fu_NL_n=fu_N(state%u_f%val*fu_exp_qh,state%temp_f%val*ft_exp_qh,state%chem_f%val*fc_exp_qh,state%t)
  ft_NL_n=ft_N(state%u_f%val*fu_exp_qh,state%temp_f%val*ft_exp_qh                           ,state%t)
  fc_NL_n=fc_N(state%u_f%val*fu_exp_qh                           ,state%chem_f%val*fc_exp_qh,state%t)
  !write(*,*)'MAXVAL OF fu_NL_n: Re', MAXVAL(real(fu_NL_n,rp)),'IMAG:',MAXVAL(AIMAG(fu_NL_n)) &
  !          ,'MAXLOC:', MAXLOC(real(fu_NL_n,rp)),MAXLOC(AIMAG(fu_NL_n))

  call set_ik_bar(state%t +dt)
  fu_NL_np1 = fu_N(state%u_f%val*fu_exp_qh+dt*fu_NL_n,state%temp_f%val*ft_exp_qh+dt*ft_NL_n &
    ,state%chem_f%val*fc_exp_qh+dt*fc_NL_n ,state%t+dt)
  ft_NL_np1 = ft_N(state%u_f%val*fu_exp_qh+dt*fu_NL_n,state%temp_f%val*ft_exp_qh+dt*ft_NL_n,state%t+dt)
  fc_NL_np1 = fc_N(state%u_f%val*fu_exp_qh+dt*fu_NL_n,state%chem_f%val*fc_exp_qh+dt*fc_NL_n,state%t+dt)
  call set_ik_bar(state%t)

  IF(ANY(IsNaN(real(fu_NL_n  ))))       write(*,*) 'sub ETD2: NAN detected in array fu_NL_n   '
  IF(ANY(IsNaN(real(fu_NL_np1))))       write(*,*) 'sub ETD2: NAN detected in array fu_NL_np1 '
  IF(ANY(IsNaN(real(ft_NL_n  ))))       write(*,*) 'sub ETD2: NAN detected in array ft_NL_n   '
  IF(ANY(IsNaN(real(ft_NL_np1))))       write(*,*) 'sub ETD2: NAN detected in array ft_NL_np1 '
  IF(ANY(IsNaN(real(fc_NL_n  ))))       write(*,*) 'sub ETD2: NAN detected in array fc_NL_n   '
  IF(ANY(IsNaN(real(fc_NL_np1))))       write(*,*) 'sub ETD2: NAN detected in array fc_NL_np1 '

  IF(ALL((real(fu_NL_n,rp).EQ.0.0_rp)))             write(*,*) 'WARNING: fu_NL     does not contribute to timestepping'
  IF(ALL((real(fu_NL_n,rp).LE.epsilon(1.0_rp))))    write(*,*) 'WARNING: fu_NL     does not contribute to timestepping'
  IF(ALL((real(fu_NL_np1,rp).EQ.0.0_rp)))           write(*,*) 'WARNING: fu_NL_np1 does not contribute to timestepping'
  IF(ALL((real(fu_NL_np1,rp).LE.epsilon(1.0_rp))))  write(*,*) 'WARNING: fu_NL_np1 does not contribute to timestepping'
  IF(ALL((real(ft_NL_n,rp).EQ.0.0_rp)))             write(*,*) 'WARNING: ft_NL     does not contribute to timestepping'
  IF(ALL((real(ft_NL_n,rp).LE.epsilon(1.0_rp))))    write(*,*) 'WARNING: ft_NL     does not contribute to timestepping'
  IF(ALL((real(ft_NL_np1,rp).EQ.0.0_rp)))           write(*,*) 'WARNING: ft_NL_np1 does not contribute to timestepping'
  IF(ALL((real(ft_NL_np1,rp).LE.epsilon(1.0_rp))))  write(*,*) 'WARNING: ft_NL_np1 does not contribute to timestepping'
  IF(ALL((real(fc_NL_n,rp).EQ.0.0_rp)))             write(*,*) 'WARNING: fc_NL     does not contribute to timestepping'
  IF(ALL((real(fc_NL_n,rp).LE.epsilon(1.0_rp))))    write(*,*) 'WARNING: fc_NL     does not contribute to timestepping'
  IF(ALL((real(fc_NL_np1,rp).EQ.0.0_rp)))           write(*,*) 'WARNING: fc_NL_np1 does not contribute to timestepping'
  IF(ALL((real(fc_NL_np1,rp).LE.epsilon(1.0_rp))))  write(*,*) 'WARNING: fc_NL_np1 does not contribute to timestepping'

    do l = 1,2
        state%u_f%val(:,:,l) = state%u_f%val(:,:,l)  *exp(D_visc*state%iki_bar_sqr%val(:,:)*dt) &
                          +dt_2*( fu_NL_n(:,:,l)     *exp(D_visc*state%iki_bar_sqr%val(:,:)*dt) &
                                 +fu_NL_np1(:,:,l)                                     )
  end do

       state%temp_f%val(:,:) = state%temp_f%val(:,:) *exp(D_therm*state%iki_bar_sqr%val(:,:)*dt) &
                          +dt_2*( ft_NL_n(:,:)       *exp(D_therm*state%iki_bar_sqr%val(:,:)*dt) &
                                 +ft_NL_np1(:,:)                                        )
 
       state%chem_f%val(:,:) = state%chem_f%val(:,:) *exp(D_comp*state%iki_bar_sqr%val(:,:)*dt) &
                         +dt_2*( fc_NL_n(:,:)        *exp(D_comp*state%iki_bar_sqr%val(:,:)*dt) &
                                +fc_NL_np1(:,:)                                        )
                                                        !NOTE the sign hidden in iki_sqr!
  IF(ANY((real(state%iki_bar_sqr%val) .GT. 0.0_rp)))  then 
    write(*,*) 'sub ETD2: real part of exponential is greater than 0.sim will explode '
    stop
  end if

  IF(ANY(IsNaN(real(state_np1%chem_f%val))))  then 
    write(*,*) 'sub ETD2: NAN detected in array fu_NL_np1   '
    stop
  end if
  IF(ANY(IsNaN(real(state_np1%temp_f%val))))  then 
    write(*,*) 'sub ETD2: NAN detected in array fu_NL_np1 '
    stop
  end if
  IF(ANY(IsNaN(real(state_np1%u_f%val))))     then 
    write(*,*) 'sub ETD2: NAN detected in array fu_NL_np1 '
    stop
  end if

  state%t           = state%t     +dt
  state%step        = state%step  +1
  !write(*,*) 'sub ETD2: step done.'
end subroutine

end module
