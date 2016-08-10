module pdgl
  ! defines the evolution equations for the fields in fourier space:
  ! du/dt = fu(...)
  ! dtemp/dt = ft(...)
  ! dchem/dt = fc(...)
  use sys_state
  use plans
  use trafo
  use nabla
  implicit none
  
contains
!---------------------------Fu---------------------------------------------------------------
function fu(u_f,temp_f,chem_f,t)
  !rhs equation for velocity field. is subdivided in several terms
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: temp_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: chem_f 
  real(kind = rp)                                           ,intent(in) :: t
  !if(debuglevel .GE. 3) write(*,*) 'function fu called'
  fu = cmplx(0.0_rp,0.0_rp,rp)
  fu = fu + fu_Nuk(u_f,t)                 !Nonlinear part
  fu = fu + fu_diff(u_f,t)                !DIFFUSION
  fu = fu + fu_buo(u_f,temp_f,chem_f,t)   !BUOYANCY 
  fu(:,:,1) = dealiase_field(fu(:,:,1))
  fu(:,:,2) = dealiase_field(fu(:,:,2))
  !_______________________________________________________________
end function 
!----------------------------------------
function fu_diff(u_f,t)
  ! calculates the diffusion term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)                     :: fu_diff
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in)          :: u_f
  real(kind = rp)                                           ,intent(in) :: t
  integer                                                               :: dir
  do dir = 1,2
     fu_diff(:,:,dir) = D_visc*( state%iki_sqr%val(:,:))*u_f(:,:,dir)
  !Note that the minus sign is intrinsicly included by the squared (ikx**2 + iky**2)
  end do
end function
!----------------------------------------
function fu_buo(u_f,temp_f,chem_f,t)
  ! calculates the bouoancy term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: temp_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: chem_f 
  real(kind = rp)                                  ,intent(in) :: t
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_buo
  fu_buo(:,:,:) =cmplx(0.0_rp,0.0_rp)
  do i=0,xdim-1
    do j=0,ydim-1
        if (.NOT.(i==0.AND.j==0)) then
          fu_buo(i,j,1) =fu_buo(i,j,1)+B_therm*&
              (-state%ikx%val(i,j)*state%iky%val(i,j)*temp_f(i,j))/state%iki_sqr%val(i,j)
          fu_buo(i,j,2) =fu_buo(i,j,2)+B_therm*&
              (-state%iky%val(i,j)*state%iky%val(i,j)*temp_f(i,j))/state%iki_sqr%val(i,j)
          fu_buo(i,j,2) =fu_buo(i,j,2)+B_therm*temp_f(i,j)


          fu_buo(i,j,1) =fu_buo(i,j,1)-B_comp*&
              (-state%ikx%val(i,j)*state%iky%val(i,j)*chem_f(i,j))/state%iki_sqr%val(i,j)
          fu_buo(i,j,2) =fu_buo(i,j,2)-B_comp*&
              (-state%iky%val(i,j)*state%iky%val(i,j)*chem_f(i,j))/state%iki_sqr%val(i,j)
          fu_buo(i,j,2) =fu_buo(i,j,2)-B_comp*chem_f(i,j)

        end if
    end do
  end do
end function
!----------------------------------------
function fu_Nuk(u_f,t)
  ! calculates the advection term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp)                                  ,intent(in) :: t
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_Nuk
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: Nuk_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: Nuk 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: div_Nuk_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: k_vec
  IF(ANY(IsNaN(real(u_f))))  then
    write(*,*) 'func fu_Nuk(): NAN detected in input array'
    stop
  end if

  k_vec(:,:,1) = state%ikx%val(:,:)
  k_vec(:,:,2) = state%iky%val(:,:)
  Nuk_f = crossp(k_vec,u)

  call transform(Nuk_f(:,:,1),Nuk(:,:,1),-1,shearing,state%t)
  call transform(Nuk_f(:,:,2),Nuk(:,:,2),-1,shearing,state%t)
  call transform(u_f(:,:,1),u(:,:,1),-1,shearing,state%t)
  call transform(u_f(:,:,2),u(:,:,2),-1,shearing,state%t)
  !call dfftw_execute_dft(ifull2D,Nuk_f(:,:,1),Nuk(:,:,1))
  !call dfftw_execute_dft(ifull2D,Nuk_f(:,:,2),Nuk(:,:,2))
  !call dfftw_execute_dft(ifull2D,u_f(:,:,1),u(:,:,1))
  !call dfftw_execute_dft(ifull2D,u_f(:,:,2),u(:,:,2))

  Nuk = crossp(-Nuk,u)

  call transform(Nuk(:,:,1),Nuk_f(:,:,1),1,shearing,state%t)
  call transform(Nuk(:,:,2),Nuk_f(:,:,2),1,shearing,state%t)
  !call dfftw_execute_dft(full2D,Nuk(:,:,1),Nuk_f(:,:,1))
  !call dfftw_execute_dft(full2D,Nuk(:,:,2),Nuk_f(:,:,2))

  div_Nuk_f = (state%ikx%val(:,:)*Nuk_f(:,:,1) &
              +state%iky%val(:,:)*Nuk_f(:,:,2) )

  do i=0,xdim-1 
    do j=0,ydim-1 
        if (.NOT.(i==0.AND.j==0)) then
        fu_Nuk(i,j,1) = Nuk_f(i,j,1)-state%ikx%val(i,j)&
                          *div_Nuk_f(i,j)/state%iki_sqr%val(i,j)
        fu_Nuk(i,j,2) = Nuk_f(i,j,2)-state%iky%val(i,j)&
                          *div_Nuk_f(i,j)/state%iki_sqr%val(i,j)
        end if
    end do
  end do


  IF(ANY(IsNaN(real(fu_Nuk))))  then
    write(*,*) 'func fu_Nuk(): NAN detected in output array'
    stop
  end if
end function
!---------------------------F_temp------------------------------------------------------------
function ft(u_f,temp_f,t)
  !rhs of temperature equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)     ,intent(in):: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2) ,intent(in):: u_f
  real(kind = rp),intent(in)                                   :: t
  ft = cmplx(0.0,0.0,rp)
  ft = ft + ft_adv(u_f,temp_f,t)  !ADVECTION
  ft = ft + ft_diff(temp_f)       !DIFFUSION
  ft = ft + ft_strat(u_f,temp_f)  !BACKGROUND stratification
  ft = dealiase_field(ft)
end function 
!--------------------------------------
function ft_adv(u_f,temp_f,t)
  implicit none
  ! temperature advection term  [nabla dot (temp*u)]
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_adv
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                   :: t

  call transform(u_f(:,:,1),state%dummy%val(:,:,1),-1,shearing,state%t)
  call transform(u_f(:,:,2),state%dummy%val(:,:,2),-1,shearing,state%t)
  call transform(temp_f,state%s_dummy%val,-1,shearing,state%t)
  !call dfftw_execute_dft(ifull2D,u_f(:,:,1) ,state%dummy%val(:,:,1))
  !call dfftw_execute_dft(ifull2D,u_f(:,:,2) ,state%dummy%val(:,:,2))
  !call dfftw_execute_dft(ifull2D,temp_f(:,:),state%s_dummy%val(:,:))

  !                                             temp * u          (realspace)
  state%dummy%val(:,:,1) = state%s_dummy%val(:,:)*state%dummy%val(:,:,1)
  state%dummy%val(:,:,2) = state%s_dummy%val(:,:)*state%dummy%val(:,:,2)

  !now trafo temp*u back to fourier space
  call transform(state%dummy%val(:,:,1),state%dummy_f%val(:,:,1), 1,shearing,state%t)
  call transform(state%dummy%val(:,:,2),state%dummy_f%val(:,:,2), 1,shearing,state%t)

  !call dfftw_execute_dft(full2D,state%dummy%val(:,:,1),state%dummy_f%val(:,:,1))
  !call dfftw_execute_dft(full2D,state%dummy%val(:,:,2),state%dummy_f%val(:,:,2))
  !state%dummy_f%val = state%dummy_f%val/real(xdim*ydim,rp) !FFTW NORM
! advection 
 ft_adv(:,:) =-(   state%ikx%val(:,:) * state%dummy_f%val(:,:,1)  &  
                  +state%iky%val(:,:) * state%dummy_f%val(:,:,2) )
end function
!--------------------------------------
function ft_diff(temp_f)
  ! temperature diffusion term in spectral domain. iki_sqr = (ikx**2 +iky**2)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_diff
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: temp_f 
  ft_diff(:,:)  = D_therm*( state%iki_sqr%val(:,:)*temp_f(:,:)) 
  !Note that the minus sign is intrinsicly included by the squared (ikx**2 + iky**2)
end function
!--------------------------------------
function ft_strat(u_f,temp_f)
  ! influence of temperature background stratification
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,2),intent(in)   :: u_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_strat
  ft_strat(:,:)  = -S_therm*u_f(:,:,2)
end function

!---------------------------F_chem------------------------------------------------------------
function fc(u_f,chem_f,t)
  !rhs of compositional equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                            :: t
  fc = cmplx(0.0,0.0,rp)
  fc = fc + fc_adv(u_f,chem_f,t)      !ADVECTION
  fc = fc + fc_diff(chem_f)           !DIFFUSION
  fc = fc + fc_strat(u_f,chem_f)      !BACKGROUND STRATIFICATION
  fc = dealiase_field(fc)
  !IF(ALL(fc ==cmplx(0.0_rp,0.0,rp)))  then
  !  write(*,*) 'func fc(): all output values are zero! '
  !  !stop
  !end if
end function 
!--------------------------------------
function fc_adv(u_f,chem_f,t)
  implicit none
  ! chem advection term  [nabla dot (chem*u)]
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_adv
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                   :: t

  call transform(u_f(:,:,1),state%dummy%val(:,:,1),-1,shearing,state%t)
  call transform(u_f(:,:,2),state%dummy%val(:,:,2),-1,shearing,state%t)
  call transform(chem_f,state%s_dummy%val,-1,shearing,state%t)
  !call dfftw_execute_dft(ifull2D,u_f(:,:,1) ,state%dummy%val(:,:,1))
  !call dfftw_execute_dft(ifull2D,u_f(:,:,2) ,state%dummy%val(:,:,2))
  !call dfftw_execute_dft(ifull2D,chem_f(:,:),state%s_dummy%val(:,:))

  !                                             chem* u          (realspace)
  state%dummy%val(:,:,1) = state%s_dummy%val(:,:)*state%dummy%val(:,:,1)
  state%dummy%val(:,:,2) = state%s_dummy%val(:,:)*state%dummy%val(:,:,2)

  !now trafo chem*u back to fourier space
  call transform(state%dummy%val(:,:,1),state%dummy_f%val(:,:,1), 1,shearing,state%t)
  call transform(state%dummy%val(:,:,2),state%dummy_f%val(:,:,2), 1,shearing,state%t)
  !call dfftw_execute_dft(full2D,state%dummy%val(:,:,1),state%dummy_f%val(:,:,1))
  !call dfftw_execute_dft(full2D,state%dummy%val(:,:,2),state%dummy_f%val(:,:,2))
  !state%dummy_f%val = state%dummy_f%val/real(xdim*ydim,rp) !FFTW NORM

! advection 
 fc_adv(:,:) =-( state%ikx%val(:,:) * state%dummy_f%val(:,:,1)  &  
                +state%iky%val(:,:) * state%dummy_f%val(:,:,2) )
end function
!--------------------------------------
function fc_diff(chem_f)
  ! compositional diffusion term in spectral domain. iki_sqr = (ikx**2 +iky**2)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_diff
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  fc_diff(:,:)  = D_therm*( state%iki_sqr%val(:,:)*chem_f(:,:))
  !Note that the minus sign is intrinsicly included by the squared (ikx**2 + iky**2)
end function
!--------------------------------------
function fc_strat(u_f,chem_f)
  ! influence of chemical background stratification
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,2),intent(in)   :: u_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_strat
  fc_strat(:,:)  = -S_comp*u_f(:,:,2)
end function
end module
