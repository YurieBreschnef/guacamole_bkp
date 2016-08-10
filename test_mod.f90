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
    !call test_periodicity()
  if(debuglevel.LE.1) write(*,*) '__________________SELFTEST DONE!_____________________'
  end subroutine

  !------------------------------------------------------------------------------------------

  subroutine NAN_check(caller)
    !checks a given system state for NAN's ; caller should look like "subroutine blabla"
    !or "function foobar"
		character(len=16),intent(inout) 		  					:: caller
    IF(ANY(IsNaN(real(state%u%val))))       write(*,*) caller,': NAN detected in array u '
    IF(ANY(IsNaN(real(state%u_f%val))))     write(*,*) caller,': NAN detected in array u_f '
    IF(ANY(IsNaN(real(state%temp%val))))    write(*,*) caller,': NAN detected in array temp '
    IF(ANY(IsNaN(real(state%temp_f%val))))  write(*,*) caller,': NAN detected in array temp_f '
    IF(ANY(IsNaN(real(state%chem%val))))    write(*,*) caller,': NAN detected in array chem'
    IF(ANY(IsNaN(real(state%chem_f%val))))  write(*,*) caller,': NAN detected in array chem_f'
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
