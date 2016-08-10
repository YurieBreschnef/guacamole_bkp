module exit_mod
  !module for all actions to be happening at program end
  use sys_state
	use const
	use plans
  implicit none

  contains
  subroutine exit_all()
    !call all special init subs
    if(debuglevel <= 1) write(*,*) '-calling exit_all()'
    call exit_plans()
    !call exit_data()     !TODO deallocate cleanly!
    if(debuglevel <= 1) write(*,*) '-done with exit_all.'
  end subroutine

  subroutine exit_plans()
    if(debuglevel <= 1) write(*,*) '  -calling exit_plans()'
	  call dfftw_destroy_plan(x_xf)
	  call dfftw_destroy_plan(xf_x)
	  call dfftw_destroy_plan(y_yf)
	  call dfftw_destroy_plan(yf_y)
    call dfftw_destroy_plan(full2D)
    call dfftw_destroy_plan(ifull2D)
    if(debuglevel <= 1) write(*,*) '  -done with exit_plans.'
  end subroutine
  
end module
