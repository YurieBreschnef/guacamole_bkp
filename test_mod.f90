module test
  !Module is used for system tests, subroutine tests and physicality checks. everything that
  !makes sure the program is running correctly goes in here 
  ! TODO: add test for flat case (zdim<4)
  ! TODO: add sub to compare pencil-decomp and fftw3D trafo
  ! TODO: add tests for pdgl-terms  (diff etc.)
  use sys_state
  use const
  use trafo
  use timestepping
  use init
  implicit none

  contains

  subroutine test_all()
  if(debuglevel.LE.1) write(*,*) '__________________SELFTEST:__________________________'
    !call test_trafo()
    call test_derivatives()
    call div_tester()
    !call test_periodicity()
  if(debuglevel.LE.1) write(*,*) '__________________SELFTEST DONE!_____________________'
  end subroutine

  !------------------------------------------------------------------------------------------

  subroutine NAN_check(caller)
    !checks a given system state for NAN's ; caller should look like "subroutine blabla"
    !or "function foobar"
		character(*),intent(inout) 		  					:: caller
    IF(ANY(IsNaN(real(state%u%val))))       write(*,*) caller,': NAN detected in array u '
    IF(ANY(IsNaN(real(state%u_f%val))))     write(*,*) caller,': NAN detected in array u_f '
    IF(ANY(IsNaN(real(state%temp%val))))    write(*,*) caller,': NAN detected in array temp '
    IF(ANY(IsNaN(real(state%temp_f%val))))  write(*,*) caller,': NAN detected in array temp_f '
    IF(ANY(IsNaN(real(state%chem%val))))    write(*,*) caller,': NAN detected in array chem'
    IF(ANY(IsNaN(real(state%chem_f%val))))  write(*,*) caller,': NAN detected in array chem_f'
  end subroutine 


 subroutine div_tester()
    type(system_state)    :: init_dummy
    type(sfield)          :: int_dummy_f
    type(sfield)          :: int_dummy
    type(sfield)          :: int1_dummy_f
    type(sfield)          :: int1_dummy
    real(kind=rp)         :: dummy_dt = 1.0e-6
    real(kind=rp)         :: maxdiv_after,maxdiv_before
    real(kind=rp)         :: b_maxdiv_after,b_maxdiv_before
    init_dummy = state
    !-----------------
    if(debuglevel.LE.1) write(*,*) '  ______________div-test__________________'
    !-------------------------------------------------------------------
    call set_ik_bar(state%t)
    !calc current maxdiv
    int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

    call transform(int_dummy_f%val,int_dummy%val,-1,1,state%t) 
    call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

    maxdiv_before   = maxval(real(int_dummy%val,real_outp_precision))
    b_maxdiv_before   = maxval(real(int1_dummy%val,real_outp_precision))

    !calc maxdiv after one euler step
    state%u_f%val = state%u_f%val + fu(state%u_f%val,state%temp_f%val,state%chem_f%val,state%t) *dummy_dt

    int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

    call transform(int_dummy_f%val,int_dummy%val,-1,0,state%t) 
    call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

    maxdiv_after = maxval(real(int_dummy%val,real_outp_precision))
    b_maxdiv_after = maxval(real(int1_dummy%val,real_outp_precision))

    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -COMPLETE'
    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() increased REGULAR div(u) by       :',&
                                         maxdiv_after-maxdiv_before,'| maxdiv_before:',&
                                          maxdiv_before,'| maxdiv_after:', maxdiv_after
    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() increased BRUCKER-div(u) by       :',&
                                         b_maxdiv_after-b_maxdiv_before,'| maxdiv_before:',&
                                          b_maxdiv_before,'| maxdiv_after:', b_maxdiv_after
    if(debuglevel.LE.1) write(*,*) ''
    ! reset to init state
    state = init_dummy 
    state%u_f%val = state%u_f%val + fu_Nuk(state%u_f%val,state%t) *dummy_dt

    int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

    call transform(int_dummy_f%val,int_dummy%val,-1,0,state%t) 
    call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

    maxdiv_after = maxval(real(int_dummy%val,real_outp_precision))
    b_maxdiv_after = maxval(real(int1_dummy%val,real_outp_precision))

    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -N_uk increased REGULAR div(u) by :',&
                                         maxdiv_after-maxdiv_before,'| maxdiv_before:',&
                                          maxdiv_before,'| maxdiv_after:', maxdiv_after
    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -N_uk increased BRUCKER-div(u) by :',&
                                         b_maxdiv_after-b_maxdiv_before,'| maxdiv_before:',&
                                          b_maxdiv_before,'| maxdiv_after:', b_maxdiv_after
    if(debuglevel.LE.1) write(*,*) ''
    
    ! reset to init state
    state = init_dummy 
    state%u_f%val = state%u_f%val + fu_diff(state%u_f%val,state%t) *dummy_dt

    int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

    call transform(int_dummy_f%val,int_dummy%val,-1,0,state%t) 
    call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

    maxdiv_after = maxval(real(int_dummy%val,real_outp_precision))
    b_maxdiv_after = maxval(real(int1_dummy%val,real_outp_precision))

    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -diff increased REGULAR div(u) by :',&
                                         maxdiv_after-maxdiv_before,'| maxdiv_before:',&
                                          maxdiv_before,'| maxdiv_after:', maxdiv_after
    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -diff increased BRUCKER-div(u) by :',&
                                         b_maxdiv_after-b_maxdiv_before,'| maxdiv_before:',&
                                          b_maxdiv_before,'| maxdiv_after:', b_maxdiv_after
    if(debuglevel.LE.1) write(*,*) ''
    ! reset to init state
    state = init_dummy 
    state%u_f%val = state%u_f%val + fu_buo(state%u_f%val,state%temp_f%val,state%chem_f%val,state%t) *dummy_dt

    int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

    call transform(int_dummy_f%val,int_dummy%val,-1,0,state%t) 
    call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

    maxdiv_after = maxval(real(int_dummy%val,real_outp_precision))
    b_maxdiv_after = maxval(real(int1_dummy%val,real_outp_precision))

    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -buo  increased REGULAR div(u) by :',&
                                         maxdiv_after-maxdiv_before,'| maxdiv_before:',&
                                          maxdiv_before,'| maxdiv_after:', maxdiv_after
    if(debuglevel.LE.1) write(*,*) '  -div-test: fu() -buo  increased BRUCKER-div(u) by :',&
                                         b_maxdiv_after-b_maxdiv_before,'| maxdiv_before:',&
                                          b_maxdiv_before,'| maxdiv_after:', b_maxdiv_after
    if(debuglevel.LE.1) write(*,*) ''

    ! reset to init state
    state = init_dummy 
    state%u_f%val = state%u_f%val + fu_shear(state%u_f%val,state%t) *dummy_dt

    int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%u_f%val(:,:,2) 
    int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

    call transform(int_dummy_f%val,int_dummy%val,-1,0,state%t) 
    call transform(int1_dummy_f%val,int1_dummy%val,-1,1,state%t) 

    maxdiv_after = maxval(real(int_dummy%val,real_outp_precision))
    b_maxdiv_after = maxval(real(int1_dummy%val,real_outp_precision))

    if(debuglevel.LE.1) write(*,*) '  -div-test: fu()-shear increased REGULAR div(u) by :',&
                                         maxdiv_after-maxdiv_before,'| maxdiv_before:',&
                                          maxdiv_before,'| maxdiv_after:', maxdiv_after
    if(debuglevel.LE.1) write(*,*) '  -div-test: fu()-shear increased BRUCKER-div(u) by :',&
                                         b_maxdiv_after-b_maxdiv_before,'| maxdiv_before:',&
                                          b_maxdiv_before,'| maxdiv_after:', b_maxdiv_after
    if(debuglevel.LE.1) write(*,*) ''
    

    !-------------------------------------------------------------------
    if(debuglevel.LE.1) write(*,*) '  -div-test: done.'
    if(debuglevel.LE.1) write(*,*) '  ________________________________________'
    ! reset to previous state
    state = init_dummy
 end subroutine

  subroutine test_derivatives()
    !tests the derivatives of the code in all three directions
    integer               :: direction=0
    type(system_state)    :: init_dummy
    init_dummy = state

  if(debuglevel.LE.1) write(*,*) '  -derivative-test: filling arrays for derivative testing..'
    do i =0,xdim-1
      do j =0,ydim-1
            state%temp%val(i,j) = cmplx(psi1(i,j),0.0_rp,rp)
          end do
      end do
    if(debuglevel.LE.1) write(*,*) '  -derivative-test: filling done.'



    if(debuglevel.LE.1) write(*,*) '  -derivative-test: transforming..'
    !call s_trafo( state%temp  ,state%temp_f  ,0)
    call dfftw_execute_dft(full2D,state%temp%val,state%temp_f%val)
    state%temp_f%val = state%temp_f%val/real(xdim*ydim,rp)
    if(debuglevel.LE.1) write(*,*) '  -derivative-test: transforming done.'



    if(debuglevel.LE.1) write(*,*) '  -derivative-test: multiplying k-values..'
    state%temp_f%val = state%temp_f%val*(state%ikx%val + state%iky%val)
    if(debuglevel.LE.1) write(*,*) '  -derivative-test: multiplying multiplying done..'


    if(debuglevel.LE.1) write(*,*) '  -derivative-test: inverse transforming..'
    !call is_trafo( state%temp_f  ,state%temp  ,0)
    call dfftw_execute_dft(ifull2D,state%temp_f%val,state%temp%val)
    if(debuglevel.LE.1) write(*,*) '  -derivative-test: inverse transforming done.'

    do i =0,xdim-1
      do j =0,ydim-1
            state%temp%val(i,j) = state%temp%val(i,j) - cmplx(div_psi1(i,j),0.0_rp,rp)
      end do
    end do

    if(debuglevel.LE.1) write(*,*) '  -derivative-test: MAX ERROR:',maxval(real(state%temp%val))
    if(maxval(real(state%temp%val)).GE.1.0e-10) then
        write(*,*) 'ERROR: derivative error is too high!'
        stop
    end if 
    if(debuglevel.LE.1) write(*,*) '  -derivative-test: MAX ERROR LOC:',maxloc(real(state%temp%val))
    if(debuglevel.LE.1) write(*,*) '  -derivative-test: done.'
  state = init_dummy
  end subroutine

  function psi1(q,r)
    !trigonometric test function for derivative-testing etc.
    use const
    implicit none
    integer,intent(in)        :: q,r
    real(kind = rp)			      :: psi1
    psi1 = 0.0_rp
    psi1 =   psi1 &
            +4.0_rp*sin(1.0_rp*(real(q,rp)/real(xdim,rp))*2.0_rp*pi)&
            +8.0_rp*cos(2.0_rp*(real(r,rp)/real(ydim,rp))*2.0_rp*pi)
  end function

  function div_psi1(q,r)
    !divergence of test function psi
    use const
    implicit none
    integer,intent(in)        :: q,r
    real(kind = rp)			      :: div_psi1
    div_psi1 = 0.0_rp
    div_psi1 =   div_psi1 &
            +4.0_rp*(2.0_rp*pi/Lx)*cos(1.0_rp*(real(q,rp)/real(xdim,rp))*2.0_rp*pi)&
           -16.0_rp*(2.0_rp*pi/Ly)*sin(2.0_rp*(real(r,rp)/real(ydim,rp))*2.0_rp*pi)
  end function
end module
