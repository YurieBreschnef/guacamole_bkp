module sys_state
	! state-module: contains all information characterizing a system state at a given time
  ! holds fields, current time and step
	use const
  implicit none

  type sfield
    !scalar field
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1):: val = 0.0_rp
    ! Chem/Temp-field value (i,j), amplitude stored in %val
  end type

  type vfield
    !vector field
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1,2):: val = 0.0_rp  
    !u(i,j,spatial_dir) acces by state%u%val(i,j, spatial direction whished)
  end type

  type kfield
    !field of k-values for derivaties in spectral domain
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1):: val = 0.0_rp  
  end type

  type system_state
    ! defines the physical system state in its entirety. Many functions and subs
    ! will use this type as return type
    type(vfield)                      :: u
    type(vfield)                      :: dummy
    type(vfield)                      :: dummy_f
    type(vfield)                      :: z_dummy_f
    type(vfield)                      :: u_f
    type(vfield)                      :: u_k1,u_k2,u_k3,u_k4	

    type(sfield)                      :: temp
    type(sfield)                      :: temp_f
    type(sfield)                      :: t_k1,t_k2,t_k3,t_k4	
    type(sfield)                      :: chem  
    type(sfield)                      :: chem_f

    type(sfield)                      :: s_dummy
    type(sfield)                      :: s_dummy_f

    type(sfield)                      :: c_dummy_f
    type(sfield)                      :: cz_dummy_f

    type(sfield)                      :: t_dummy_f
    type(sfield)                      :: tz_dummy_f
    type(sfield)                      :: c_k1,c_k2,c_k3,c_k4

  	type(kfield)                      :: ikx,iky! k's for deriv
  	type(kfield)                      :: ikx_sqr,iky_sqr! k's for deriv
  	type(kfield)                      :: iki_sqr

  	type(kfield)                      :: ikx_bar,iky_bar! k's for deriv
  	type(kfield)                      :: ikx_bar_sqr,iky_bar_sqr! k's for deriv
  	type(kfield)                      :: iki_bar_sqr

  	real(kind = rp)								   	:: t	  = 0.0_rp	! time variable
  	integer									    	   	:: step = 0 			! acute step of sim
  end type

  !MAIN STATE VARIABLE:
  type(system_state)                                          ::state


contains

  function measure_Ekin()
    !measures a absolute value for kinetic energy within system.
    !TODO - devide by relative density
    real(kind=real_outp_precision)            ::measure_Ekin
    real(kind=rp)                             ::Ekin
    Ekin = 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      Ekin = Ekin+(real(state%u%val(i,j,1))**2 &
                  +real(state%u%val(i,j,2))**2)
    end do
    end do
    measure_Ekin = real(Ekin,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_vmax()
    ! returns the highest absolute velocitiy within the system
    real(kind=real_outp_precision)            ::measure_vmax
    real(kind=rp)                             ::vmax
    vmax =0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
    do l=1,2
        if(vmax < real(state%u%val(i,j,l),rp)) vmax = real(state%u%val(i,j,l),rp)
    end do
    end do
    end do
    measure_vmax = vmax
  end function

  function measure_Epot()
    !measures a absolute value for kinetic energy within system.
    !TODO - devide by relative density
    real(kind=real_outp_precision)            ::measure_Epot
    real(kind=rp)                             ::Epot
    Epot= 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      Epot= Epot+( B_therm*state%temp%val(i,j) - B_comp*state%chem%val(i,j))*real(j)/real(ydim)*Ly
    end do
    end do
    measure_Epot = real(Epot,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_av_temp()
    !measures a relative value for temp  within system.
    real(kind=real_outp_precision)            ::measure_av_temp
    real(kind=rp)                             ::tmp
    tmp= 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      tmp= tmp+(state%temp%val(i,j)) 
    end do
    end do
    measure_av_temp= real(tmp,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_av_chem()
    !measures a relative value for temp  within system.
    real(kind=real_outp_precision)            ::measure_av_chem
    real(kind=rp)                             ::chem
    chem= 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      chem= chem+(state%chem%val(i,j)) 
    end do
    end do
    measure_av_chem= real(chem,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_u_rms()
    !measures a relative value for temp  within system.
    real(kind=real_outp_precision)            ::measure_u_rms
    real(kind=rp)                             ::u_rms
    u_rms= 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
    do k=1,2
      u_rms= u_rms+(state%u%val(i,j,1)**2 + state%u%val(i,j,2)**2)
    end do
    end do
    end do
    measure_u_rms= sqrt(real((u_rms),real_outp_precision)/real(xdim*ydim,real_outp_precision))
  end function


  subroutine set_ik_bar(ktime)
    ! resets the ik_bar wave vectors for a given time
    ! note how the x component is not really used, but still there for generality and possible
    ! future changes
    real(kind = rp)                              :: ktime
    !TODO very ineffective to reset k's if shearing is off.
    state%ikx_bar%val(:,:) = state%ikx%val(:,:) 
    state%iky_bar%val(:,:) = state%iky%val(:,:) 
    state%ikx_bar_sqr%val(:,:) = state%ikx_bar%val(:,:)**2
    state%ikx_bar_sqr%val(0,0) = epsilon(1.0_rp)
    state%iky_bar_sqr%val(:,:) = state%iky_bar%val(:,:)**2
    state%iky_bar_sqr%val(0,0) = epsilon(1.0_rp)
    state%iki_bar_sqr%val(:,:) = state%ikx_bar%val(:,:)**2 + state%iky_bar%val(:,:)**2
    state%iki_bar_sqr%val(0,0) = epsilon(1.0_rp)
    if(shearing ==1) then
      state%ikx_bar%val(:,:) = state%ikx%val(:,:) 
      state%iky_bar%val(:,:) = state%iky%val(:,:) - shear*ktime*state%ikx%val(:,:)

      state%ikx_bar_sqr%val(:,:) = state%ikx_bar%val(:,:)**2
      state%iky_bar_sqr%val(:,:) = state%iky_bar%val(:,:)**2

      state%iki_bar_sqr%val(:,:) = state%ikx_bar%val(:,:)**2 + state%iky_bar%val(:,:)**2
    end if
  end subroutine

end module
