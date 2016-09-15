module pdgl
  ! defines the evolution equations for the fields in fourier space:
  ! du/dt = fu(...)
  ! dtemp/dt = ft(...)
  ! dchem/dt = fc(...)
  use sys_state
  use plans
  use trafo
  use nabla
  use benchmark
  implicit none
  
contains
!---------------------------Fu---------------------------------------------------------------
function fu(u_f,temp_f,chem_f,t)
  !rhs equation for velocity field. is subdivided in several terms
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: temp_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: chem_f 
  real(kind = rp)                                  ,intent(in) :: t
  if(benchmarking ==1) call cpu_time(bm_fu_starttime)
  !if(debuglevel .GE. 3) write(*,*) 'function fu called'
  call set_ik_bar(t)
  fu = cmplx(0.0_rp,0.0_rp,rp)
  fu = fu + fu_Nuk(u_f,t)                 !Nonlinear part
  fu = fu + fu_diff(u_f,t)                !DIFFUSION
  fu = fu + fu_buo(u_f,temp_f,chem_f,t)   !BUOYANCY 
  fu = fu + fu_shear(u_f,t)               !SHEAR

  !fu(:,:,1) = dealiase_field(fu(:,:,1))
  !fu(:,:,2) = dealiase_field(fu(:,:,2))
  if(benchmarking ==1) call cpu_time(bm_fu_endtime)
end function 
!----------------------------------------
function fu_L(u_f,t)
  !linear operator of function fu = L+ N
  ! maybe used for ETD
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_L
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp)                                  ,intent(in) :: t
  !if(debuglevel .GE. 3) write(*,*) 'function fu called'
  call set_ik_bar(t)
  fu_L = cmplx(0.0_rp,0.0_rp,rp)
  fu_L = fu_L + fu_diff(u_f,t)                !DIFFUSION
  !fu_L(:,:,1) = dealiase_field(fu_L(:,:,1))
  !fu_L(:,:,2) = dealiase_field(fu_L(:,:,2))
end function 
!----------------------------------------
function fu_N(u_f,temp_f,chem_f,t)
  !nonlinear operator of function fu = L+ N
  !maybe used for ETD
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_N
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: temp_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: chem_f 
  real(kind = rp)                                  ,intent(in) :: t
  IF(ALL((real(u_f,rp).EQ.0.0_rp)))  then
    write(*,*) 'func fu_N():WARNING:  all zeroes detected in input array u_f'
  end if
  IF(ALL((real(chem_f,rp).EQ.0.0_rp)))  then
    write(*,*) 'func fu_N():WARNING:  all zeroes detected in input array chem_f'
  end if
  IF(ALL((real(temp_f,rp).EQ.0.0_rp)))  then
    write(*,*) 'func fu_N():WARNING: all zeroes detected in input array temp_f'
  end if

  call set_ik_bar(t)
  fu_N = cmplx(0.0_rp,0.0_rp,rp)
  fu_N = fu_N + fu_Nuk(u_f,t)                 !Nonlinear part
  fu_N = fu_N + fu_buo(u_f,temp_f,chem_f,t)   !BUOYANCY 
  if(shearing==1) then
    fu_N = fu_N + fu_shear(u_f,t)               !SHEAR
  end if

  !fu_N(:,:,1) = dealiase_field(fu_N(:,:,1))
  !fu_N(:,:,2) = dealiase_field(fu_N(:,:,2))

  IF(ALL((real(fu_N,rp).EQ.0.0_rp)))write(*,*) 'WARNING: fu_N does not contribute to pdgl!'
end function 
!----------------------------------------
function fu_shear(u_f,t)
  ! calculates the shear-related term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp)                                  ,intent(in) :: t
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_shear
  if(benchmarking ==1) call cpu_time(bm_fu_shear_starttime)
  fu_shear(:,:,:) =cmplx(0.0_rp,0.0_rp)
  if(shearing ==1) then
    do i=0,xdim-1
      do j=0,ydim-1
          if (.NOT.(i==0.AND.j==0)) then
             fu_shear(i,j,1) = -shear*state%u_f%val(i,j,2)
             fu_shear(i,j,2) = cmplx(0.0_rp,0.0_rp)
             fu_shear(i,j,1) = fu_shear(i,j,1)+2.0_rp*((state%ikx_bar%val(i,j)*state%ikx_bar%val(i,j))&
                                          /state%iki_bar_sqr%val(i,j))*shear*state%u_f%val(i,j,2)
             fu_shear(i,j,2) = fu_shear(i,j,2)+2.0_rp*((state%iky_bar%val(i,j)*state%ikx_bar%val(i,j))&
                                          /state%iki_bar_sqr%val(i,j))*shear*state%u_f%val(i,j,2)
             !NOTE: minus sign is due to imag included in ikx,iky and their multiplikation
          end if
      end do
    end do
  end if
  if(benchmarking ==1) call cpu_time(bm_fu_shear_endtime)
end function
!----------------------------------------
function fu_diff(u_f,t)
  ! calculates the diffusion term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)                     :: fu_diff
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in)          :: u_f
  real(kind = rp)                                           ,intent(in) :: t
  integer                                                               :: dir
  if(benchmarking ==1) call cpu_time(bm_fu_diff_starttime)
  call set_ik_bar(t)
  fu_diff(:,:,1) = D_visc*(state%iki_bar_sqr%val(:,:))*u_f(:,:,1)
  fu_diff(:,:,2) = D_visc*(state%iki_bar_sqr%val(:,:))*u_f(:,:,2)
  !note minus sign intrinsic in (i*k)**2
  if(benchmarking ==1) call cpu_time(bm_fu_diff_endtime)
end function
!----------------------------------------
function fu_buo(u_f,temp_f,chem_f,t)
  ! calculates the bouoancy term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: temp_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)    ,intent(in) :: chem_f 
  real(kind = rp)                                  ,intent(in) :: t
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_buo
  if(benchmarking ==1) call cpu_time(bm_fu_buo_starttime)
  fu_buo(:,:,:) =cmplx(0.0_rp,0.0_rp)
  do i=0,xdim-1
    do j=0,ydim-1
        if (.NOT.(i==0.AND.j==0)) then
          !TEMPERATURE PART
          fu_buo(i,j,1) =fu_buo(i,j,1)+B_therm*&
              (-state%ikx_bar%val(i,j)*(state%iky_bar%val(i,j)*temp_f(i,j)))/state%iki_bar_sqr%val(i,j)
          fu_buo(i,j,2) =fu_buo(i,j,2)+B_therm*&
              (-state%iky_bar%val(i,j)*(state%iky_bar%val(i,j)*temp_f(i,j)))/state%iki_bar_sqr%val(i,j)

          fu_buo(i,j,2) =fu_buo(i,j,2)+B_therm*temp_f(i,j)
          !CHEMICAL PART
          fu_buo(i,j,1) =fu_buo(i,j,1)-B_comp*&
              (-state%ikx_bar%val(i,j)*(state%iky_bar%val(i,j)*chem_f(i,j)))/state%iki_bar_sqr%val(i,j)
          fu_buo(i,j,2) =fu_buo(i,j,2)-B_comp*&
              (-state%iky_bar%val(i,j)*(state%iky_bar%val(i,j)*chem_f(i,j)))/state%iki_bar_sqr%val(i,j)

          fu_buo(i,j,2) =fu_buo(i,j,2)-B_comp*chem_f(i,j)
          !note the minus signs within ik- variables
        end if
    end do
  end do
  if(benchmarking ==1) call cpu_time(bm_fu_buo_endtime)
end function
!----------------------------------------
function fu_Nuk(u_f,t)
  ! calculates the advection term of velocity evolution equation in fourier space
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp)                                  ,intent(in) :: t
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: fu_Nuk
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: Nuk_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: omega
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: omega_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: Nuk 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: div_Nuk_f
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            :: u
  if(benchmarking ==1) call cpu_time(bm_fu_Nuk_starttime)
  !IF(ANY(IsNaN(real(u_f))))  then
  !  write(*,*) 'func fu_Nuk(): NAN detected in input array u_f'
  !  stop
  !end if
  call set_ik_bar(t)
  !do crossproduct
  omega_f = crossp(state%k_vec%val,u_f)
  !IF(ANY(IsNaN(real(Nuk_f))))  then
  !  write(*,*) 'func fu_Nuk(): NAN detected after crossp of nabla x u'
  !  stop
  !end if

  call transform(omega_f(:,:,1),omega(:,:,1),-1,shearing,t)
  call transform(omega_f(:,:,2),omega(:,:,2),-1,shearing,t)
  call transform(u_f(:,:,1),u(:,:,1),-1,shearing,t)
  call transform(u_f(:,:,2),u(:,:,2),-1,shearing,t)

  !do crossproduct
  Nuk = -crossp(omega,u)
  !IF(ANY(IsNaN(real(Nuk))))  then
  !  write(*,*) 'func fu_Nuk(): NAN detected in Nuk after crossp in realspace'
  !  stop
  !end if

  call transform(Nuk(:,:,1),Nuk_f(:,:,1),1,shearing,t)
  call transform(Nuk(:,:,2),Nuk_f(:,:,2),1,shearing,t)
  !IF(ANY(IsNaN(real(Nuk_f))))  then
  !  write(*,*) 'func fu_Nuk(): NAN detected in Nuk_f after trafo'
  !  stop
  !end if

 div_Nuk_f = (state%ikx_bar%val(:,:)*Nuk_f(:,:,1) &
             +state%iky_bar%val(:,:)*Nuk_f(:,:,2) )

  !IF(ANY(IsNaN(real(Nuk_f))))  then
  !  write(*,*) 'func fu_Nuk(): NAN detected before k mult '
  !  !stop
  !end if

  do i=0,xdim-1 
    do j=0,ydim-1 
        if (.NOT.(i==0.AND.j==0)) then
        !IF(ANY(IsNaN(real(fu_Nuk(i,j,1)))))  then
        !  write(*,*) 'func fu_Nuk(): NAN detected at pos (bef):',i,j,1
        !  stop
        !end if
        !if(state%iki_bar_sqr%val(i,j).EQ.cmplx(0.0_rp,0.0_rp)) then
        !  write(*,*) 'func fu_Nuk(): zero in k-vec NAN imminent at pos (bef):',i,j
        !  stop
        !end if 
        IF(.NOT.(IsNaN(real(1.0_rp/state%iki_bar_sqr%val(i,j)))))  then
          fu_Nuk(i,j,1) =Nuk_f(i,j,1)-state%ikx_bar%val(i,j)&
                            *div_Nuk_f(i,j)/state%iki_bar_sqr%val(i,j)
          fu_Nuk(i,j,2) =Nuk_f(i,j,2)-state%iky_bar%val(i,j)&
                            *div_Nuk_f(i,j)/state%iki_bar_sqr%val(i,j)
          !note minus sign resulting from i factor in ikx,iky..
        else
          write(*,*) 'func fu_Nuk(): NAN detected (before 1.0_rp/iki_bar_sqr)at  at pos:',i,j
          stop
        end if
        end if
    end do
  end do


  !IF(ANY(IsNaN(real(fu_Nuk))))  then
  !  write(*,*) 'func fu_Nuk(): NAN detected in output array'
  !  stop
  !end if
  if(benchmarking ==1) call cpu_time(bm_fu_Nuk_endtime)
end function
!---------------------------F_temp------------------------------------------------------------

function ft(u_f,temp_f,t)
  !rhs of temp equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                            :: t
  if(benchmarking ==1) call cpu_time(bm_ft_starttime)
  call set_ik_bar(t)
  ft = cmplx(0.0,0.0,rp)
  ft = ft + ft_adv(u_f,temp_f,t)      !ADVECTION
  ft = ft + ft_diff(temp_f)           !DIFFUSION
  ft = ft + ft_strat(u_f,temp_f)      !BACKGROUND STRATIFICATION
  !ft = dealiase_field(ft)
  !IF(ALL(ft ==cmplx(0.0_rp,0.0,rp)))  then
  !  write(*,*) 'func ft(): all output values are zero! '
  !  !stop
  !end if
  if(benchmarking ==1) call cpu_time(bm_ft_endtime)
end function 
!--------------------------------------
function ft_L(temp_f,t)
  !linear part of temp evolution equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_L
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)     ,intent(in):: temp_f 
  real(kind = rp),intent(in)                                   :: t
  call set_ik_bar(t)
  !IF(ANY(IsNaN(real(temp_f))))  then
  !  write(*,*) 'func ft_L(): NAN detected in input array'
  !  stop
  !end if
  ft_L = cmplx(0.0,0.0,rp)
  ft_L = ft_L + ft_diff(temp_f)       !DIFFUSION
  !ft_L = dealiase_field(ft_L)
end function 
!--------------------------------------
function ft_N(u_f,temp_f,t)
  ! nonlinear part of temp evolution equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_N
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)     ,intent(in):: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2) ,intent(in):: u_f
  real(kind = rp),intent(in)                                   :: t
  call set_ik_bar(t)
  !IF(ALL((real(u_f,rp).EQ.0.0_rp)))  then
  !  write(*,*) 'func ft_N(): WARNING: ALL ZEROES detected in input array u_f'
  !end if
  !IF(ALL((real(temp_f,rp).EQ.0.0_rp)))  then
  !  write(*,*) 'func ft_N(): WARNING ALL ZEROES detected in input array temp_f'
  !end if
  !IF(ANY(IsNaN(real(u_f))).OR.ANY(IsNaN(real(temp_f))))  then
  !  write(*,*) 'func ft_N(): NAN detected in input array'
  !  stop
  !end if
  ft_N = cmplx(0.0,0.0,rp)
  ft_N = ft_N + ft_adv(u_f,temp_f,t)  !ADVECTION
  ft_N = ft_N + ft_strat(u_f,temp_f)  !BACKGROUND stratification
  !ft_N = dealiase_field(ft_N)
  IF(ALL((real(ft_N,rp).EQ.0.0_rp)))write(*,*) 'WARNING: ft_N does not contribute to pdgl!'
end function 
!--------------------------------------
function ft_adv(u_f,temp_f,t)
  implicit none
  ! temp advection term  [nabla dot (temp*u)]
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_adv
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                   :: t
  if(benchmarking ==1) call cpu_time(bm_ft_adv_starttime)
  !IF(ANY(IsNaN(real(u_f))).OR.ANY(IsNaN(real(temp_f))))  then
  !  write(*,*) 'func ft_adv(): NAN detected in input array'
  !  stop
  !end if

  call transform(u_f(:,:,1),state%dummy%val(:,:,1),-1,shearing,t)
  call transform(u_f(:,:,2),state%dummy%val(:,:,2),-1,shearing,t)
  call transform(temp_f,state%s_dummy%val,-1,shearing,t)

  !                                            temp_f* u          (realspace)
  state%dummy%val(:,:,1) = state%s_dummy%val(:,:)*state%dummy%val(:,:,1)
  state%dummy%val(:,:,2) = state%s_dummy%val(:,:)*state%dummy%val(:,:,2)

  !now trafo temp*u back to fourier space
  call transform(state%dummy%val(:,:,1),state%dummy_f%val(:,:,1), 1,shearing,t)
  call transform(state%dummy%val(:,:,2),state%dummy_f%val(:,:,2), 1,shearing,t)

! advection 
 ft_adv(:,:) =-( state%ikx_bar%val(:,:) * state%dummy_f%val(:,:,1)  &  
                +state%iky_bar%val(:,:) * state%dummy_f%val(:,:,2) )
  if(benchmarking ==1) call cpu_time(bm_ft_adv_endtime)
end function
!--------------------------------------
function ft_diff(temp_f)
  ! compositional diffusion term in spectral domain. iki_sqr = (ikx**2 +iky**2)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_diff
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: temp_f 
  if(benchmarking ==1) call cpu_time(bm_ft_diff_starttime)
  !IF(ANY(IsNaN(real(temp_f))))  then
  !  write(*,*) 'func ft_diff(): NAN detected in input array'
  !  stop
  !end if
  !ft_diff(:,:)  = D_therm*( state%iki_sqr%val(:,:)*temp_f(:,:))
  ft_diff(:,:)  = D_therm*( state%iki_bar_sqr%val(:,:)*temp_f(:,:))
  !Note that the minus sign is intrinsicly included by the squared (ikx**2 + iky**2)
  if(benchmarking ==1) call cpu_time(bm_ft_diff_endtime)
end function
!--------------------------------------
function ft_strat(u_f,temp_f)
  ! influence of chemical background stratification
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,2),intent(in)   :: u_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: temp_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: ft_strat
  if(benchmarking ==1) call cpu_time(bm_ft_strat_starttime)
  ft_strat(:,:)  = -S_therm*u_f(:,:,2)
  if(benchmarking ==1) call cpu_time(bm_ft_strat_endtime)
end function
!---------------------------F_chem------------------------------------------------------------
function fc(u_f,chem_f,t)
  !rhs of compositional equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                            :: t
  if(benchmarking ==1) call cpu_time(bm_fc_starttime)
  call set_ik_bar(t)
  fc = cmplx(0.0,0.0,rp)
  fc = fc + fc_adv(u_f,chem_f,t)      !ADVECTION
  fc = fc + fc_diff(chem_f)           !DIFFUSION
  fc = fc + fc_strat(u_f,chem_f)      !BACKGROUND STRATIFICATION
  !fc = dealiase_field(fc)
  !IF(ALL(fc ==cmplx(0.0_rp,0.0,rp)))  then
  !  write(*,*) 'func fc(): all output values are zero! '
  !  !stop
  !end if
  if(benchmarking ==1) call cpu_time(bm_fc_endtime)
end function 
!--------------------------------------
function fc_L(chem_f,t)
  !linear part of chem evolution equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_L
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)     ,intent(in):: chem_f 
  real(kind = rp),intent(in)                                   :: t
  call set_ik_bar(t)
  !IF(ANY(IsNaN(real(chem_f))))  then
  !  write(*,*) 'func fc_L(): NAN detected in input array'
  !  stop
  !end if
  fc_L = cmplx(0.0,0.0,rp)
  fc_L = fc_L + fc_diff(chem_f)       !DIFFUSION
  !fc_L = dealiase_field(fc_L)
end function 
!--------------------------------------
function fc_N(u_f,chem_f,t)
  ! nonlinear part of chem evolution equation
  !TODO make transforms as effective as possible
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_N
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)     ,intent(in):: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2) ,intent(in):: u_f
  real(kind = rp),intent(in)                                   :: t
  call set_ik_bar(t)
  !IF(ALL((real(u_f,rp).EQ.0.0_rp)))  then
  !  write(*,*) 'func fc_N(): WARNING: ALL ZEROES detected in input array u_f'
  !end if
  !IF(ALL((real(chem_f,rp).EQ.0.0_rp)))  then
  !  write(*,*) 'func fc_N(): WARNING ALL ZEROES detected in input array chem_f'
  !end if
  !IF(ANY(IsNaN(real(u_f))).OR.ANY(IsNaN(real(chem_f))))  then
  !  write(*,*) 'func fc_N(): NAN detected in input array'
  !  stop
  !end if
  fc_N = cmplx(0.0,0.0,rp)
  fc_N = fc_N + fc_adv(u_f,chem_f,t)  !ADVECTION
  fc_N = fc_N + fc_strat(u_f,chem_f)  !BACKGROUND stratification
  !fc_N = dealiase_field(fc_N)
  IF(ALL((real(fc_N,rp).EQ.0.0_rp)))write(*,*) 'WARNING: fc_N does not contribute to pdgl!'
end function 
!--------------------------------------
function fc_adv(u_f,chem_f,t)
  implicit none
  ! chem advection term  [nabla dot (chem*u)]
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_adv
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2),intent(in) :: u_f
  real(kind = rp),intent(in)                                   :: t
  if(benchmarking ==1) call cpu_time(bm_fc_adv_starttime)
  !IF(ANY(IsNaN(real(u_f))).OR.ANY(IsNaN(real(chem_f))))  then
  !  write(*,*) 'func fc_adv(): NAN detected in input array'
  !  stop
  !end if

  call transform(u_f(:,:,1),state%dummy%val(:,:,1),-1,shearing,t)
  call transform(u_f(:,:,2),state%dummy%val(:,:,2),-1,shearing,t)
  call transform(chem_f,state%s_dummy%val,-1,shearing,t)

  !                                             chem* u          (realspace)
  state%dummy%val(:,:,1) = state%s_dummy%val(:,:)*state%dummy%val(:,:,1)
  state%dummy%val(:,:,2) = state%s_dummy%val(:,:)*state%dummy%val(:,:,2)

  !now trafo chem*u back to fourier space
  call transform(state%dummy%val(:,:,1),state%dummy_f%val(:,:,1), 1,shearing,t)
  call transform(state%dummy%val(:,:,2),state%dummy_f%val(:,:,2), 1,shearing,t)

! advection 
 fc_adv(:,:) =-( state%ikx_bar%val(:,:) * state%dummy_f%val(:,:,1)  &  
                +state%iky_bar%val(:,:) * state%dummy_f%val(:,:,2) )
  if(benchmarking ==1) call cpu_time(bm_fc_adv_endtime)
end function
!--------------------------------------
function fc_diff(chem_f)
  ! compositional diffusion term in spectral domain. iki_sqr = (ikx**2 +iky**2)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_diff
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  if(benchmarking ==1) call cpu_time(bm_fc_diff_starttime)
  !IF(ANY(IsNaN(real(chem_f))))  then
  !  write(*,*) 'func fc_diff(): NAN detected in input array'
  !  stop
  !end if
  !fc_diff(:,:)  = D_comp*( state%iki_sqr%val(:,:)*chem_f(:,:))
  fc_diff(:,:)  = D_comp*( state%iki_bar_sqr%val(:,:)*chem_f(:,:))
  !Note that the minus sign is intrinsicly included by the squared (ikx**2 + iky**2)
  if(benchmarking ==1) call cpu_time(bm_fc_diff_endtime)
end function
!--------------------------------------
function fc_strat(u_f,chem_f)
  ! influence of chemical background stratification
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,2),intent(in)   :: u_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1),intent(in)     :: chem_f 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)                :: fc_strat
  if(benchmarking ==1) call cpu_time(bm_fc_strat_starttime)
  fc_strat(:,:)  = -S_comp*u_f(:,:,2)
  if(benchmarking ==1) call cpu_time(bm_fc_strat_endtime)
end function
end module
