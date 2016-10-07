
subroutine ETD2_step()
!perfomrs ETD2 timestep
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_RHS_nm1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_exp_qh ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_RHS_n_factor 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u_RHS_nm1_factor 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_RHS_nm1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_exp_qh ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_RHS_n_factor 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: t_RHS_nm1_factor 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_RHS_nm1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_exp_qh ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_RHS_n_factor 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: c_RHS_nm1_factor 

  if(state%step ==0) then
    !the very first step of sim is done in euler way because ETD2 needs the past time variable 
    !RHS_nm1 (read as n minus one)
    ! set initial state as the old state
    state_nm1 = state
    call set_ik_bar(sheartime) 
    ! set nm1 variables
    u_RHS_nm1 = fu_N(state_nm1%u_f%val,state_nm1%temp_f%val,state_nm1%chem_f%val,state_nm1%t)
    t_RHS_nm1 = ft_N(state_nm1%u_f%val,state_nm1%temp_f%val                     ,state_nm1%t)
    c_RHS_nm1 = fc_N(state_nm1%u_f%val                     ,state_nm1%chem_f%val,state_nm1%t)
    ! now do first step as euler step
	  state%u_f%val    =state%u_f%val    + dt*fu(state%u_f%val ,state%temp_f%val,state%chem_f%val,state%t) 
	  state%temp_f%val =state%temp_f%val + dt*ft(state%u_f%val ,state%temp_f%val                 ,state%t)     
	  state%chem_f%val =state%chem_f%val + dt*fc(state%u_f%val                  ,state%chem_f%val,state%t)     
	  state%t=state%t+dt
	  state%step=state%step+1
    write(*,*) 'initial timestep done with euler'
    return
  end if
  ! calc  RHS_nm1 and RHS_n (read as n minus one)
  u_RHS_nm1 = fu_N(state_nm1%u_f%val,state_nm1%temp_f%val,state_nm1%chem_f%val,state_nm1%t)
  t_RHS_nm1 = ft_N(state_nm1%u_f%val,state_nm1%temp_f%val                     ,state_nm1%t)
  c_RHS_nm1 = fc_N(state_nm1%u_f%val                     ,state_nm1%chem_f%val,state_nm1%t)
  u_RHS_n   = fu_N(state%u_f%val,state%temp_f%val,state%chem_f%val,state%t)
  t_RHS_n   = ft_N(state%u_f%val,state%temp_f%val                 ,state%t)
  c_RHS_n   = fc_N(state%u_f%val                 ,state%chem_f%val,state%t)
  ! set q-values for exponent, note the minus sign in iki_sqr
  u_q(:,:,1) = D_visc  *state%iki_bar_sqr%val(:,:)
  u_q(:,:,2) = D_visc  *state%iki_bar_sqr%val(:,:)
  t_q = D_therm *state%iki_bar_sqr%val
  c_q = D_comp  *state%iki_bar_sqr%val
  ! calc exponentials for multiplication
  u_exp_qh = exp(u_q*dt)
  t_exp_qh = exp(t_q*dt)
  c_exp_qh = exp(c_q*dt)
  ! set factors to zero (especially the (0,0) index to prevent NAN by division)
  u_RHS_n_factor  = cmplx(0.0_rp,0.0_rp,rp)
  u_RHS_n_factor  = cmplx(0.0_rp,0.0_rp,rp)
  t_RHS_n_factor  = cmplx(0.0_rp,0.0_rp,rp)
  c_RHS_n_factor  = cmplx(0.0_rp,0.0_rp,rp)
  u_RHS_nm1_factor = cmplx(0.0_rp,0.0_rp,rp)
  u_RHS_nm1_factor = cmplx(0.0_rp,0.0_rp,rp)
  t_RHS_nm1_factor = cmplx(0.0_rp,0.0_rp,rp)
  c_RHS_nm1_factor = cmplx(0.0_rp,0.0_rp,rp)
  ! calc factors which go before RHS_n  and RHS_nm1 when timestep is made
  do i=0,xdim-1
   do j=0,ydim-1
    if(.NOT.(i==0.AND.j==0))then
      if(.NOT.(IsNAN(real(1.0_rp/(dt*u_q(i,j,1)**2)))))then
        u_RHS_n_factor(i,j,1)  =((1.0_rp + dt * u_q(i,j,1)) -1.0_rp -2.0_rp*dt*u_q(i,j,1))/(dt*u_q(i,j,1)**2)
        u_RHS_n_factor(i,j,2)  =((1.0_rp + dt * u_q(i,j,2)) -1.0_rp -2.0_rp*dt*u_q(i,j,2))/(dt*u_q(i,j,2)**2)

        t_RHS_n_factor(i,j)    =((1.0_rp + dt * t_q(i,j))   -1.0_rp -2.0_rp*dt*t_q(i,j))  /(dt*t_q(i,j)**2)
        c_RHS_n_factor(i,j)    =((1.0_rp + dt * c_q(i,j))   -1.0_rp -2.0_rp*dt*c_q(i,j))  /(dt*c_q(i,j)**2)

        u_RHS_nm1_factor(i,j,1)=(-u_exp_qh(i,j,1) + 1.0_rp + dt * u_q(i,j,1))/(dt*u_q(i,j,1)**2)
        u_RHS_nm1_factor(i,j,2)=(-u_exp_qh(i,j,2) + 1.0_rp + dt * u_q(i,j,2))/(dt*u_q(i,j,2)**2)

        t_RHS_nm1_factor(i,j)  = (-t_exp_qh(i,j)   + 1.0_rp + dt * t_q(i,j))  /(dt*t_q(i,j)**2)
        c_RHS_nm1_factor(i,j)  = (-c_exp_qh(i,j)   + 1.0_rp + dt * c_q(i,j))  /(dt*c_q(i,j)**2)
      else 
        write(*,*) 'ETD2: sub skipped the calc of RHS-factors because 1/(dt*k**4) is nan!'
      end if
    end if
   end do
  end do

  do i=0,xdim-1
   do j=0,ydim-1
    if(.NOT.(i==0.AND.j==0))then
    if(IsNaN(real(1.0_rp/(dt*u_q(i,j,1)**2))))  then
      write(*,*) 'NAN detected in 1.0/u_q**2 at pos:',i,j,'u_q(i,j,1)=',u_q(i,j,1)
    end if

    if(IsNaN(real(u_RHS_n_factor(i,j,1))))  then
      write(*,*) 'NAN detected in u_RHS_n at pos:',i,j,'u_RHS_n(i,j,1)=',u_RHS_n_factor(i,j,1)
    end if

    if(IsNaN(real(u_RHS_nm1_factor(i,j,1))))  then
      write(*,*) 'NAN detected in u_RHS_n at pos:',i,j,'u_RHS_n(i,j,1)=',u_RHS_nm1_factor(i,j,1)
    end if

    if(IsNaN(real(1.0_rp/(dt*u_q(i,j,2)**2))))  then
      write(*,*) 'NAN detected in 1.0/u_q**2 at pos:',i,j,'u_q(i,j,2)=',u_q(i,j,2)
    end if
    end if
   end do
  end do

  ! set the factors to zero where 1/q_sqr would produce NAN otherwise
  u_RHS_n_factor(0,0,:)   = cmplx(0.0_rp,0.0_rp,rp) 
  t_RHS_n_factor(0,0)     = cmplx(0.0_rp,0.0_rp,rp) 
  c_RHS_n_factor(0,0)     = cmplx(0.0_rp,0.0_rp,rp) 
  u_RHS_nm1_factor(0,0,:) = cmplx(0.0_rp,0.0_rp,rp) 
  t_RHS_nm1_factor(0,0)   = cmplx(0.0_rp,0.0_rp,rp) 
  c_RHS_nm1_factor(0,0)   = cmplx(0.0_rp,0.0_rp,rp) 

  !set the current state to be the new old step
  state_nm1 = state
  ! make step
  state%u_f%val   =state%u_f%val    *u_exp_qh + u_RHS_n*u_RHS_n_factor + u_RHS_nm1*u_RHS_nm1_factor
  state%temp_f%val=state%temp_f%val *t_exp_qh + t_RHS_n*t_RHS_n_factor + t_RHS_nm1*t_RHS_nm1_factor
  state%chem_f%val=state%chem_f%val *c_exp_qh + c_RHS_n*c_RHS_n_factor + c_RHS_nm1*c_RHS_nm1_factor

  state%t           = state%t     +dt
  state%step        = state%step  +1
  write(*,*) 'ETD step done.'
end subroutine
