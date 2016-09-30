program guacamole 
  !main program part, the main timestepping loop is implemented here, calling the actual
  !stepping routines in the timestepping_mod, and writing to files by calling write_all()
  !from IO_module.

  !use omp_lib
  !use iso_c_binding
  use plans
  use sys_state
  use init
  use exit_mod
  use timestepping
  use IO_mod
  use trafo
  use test
  use benchmark
  implicit none
  if(debuglevel .GE. 1) write(*,*) '__________________START____________________________________'

  call init_all()
  call test_all()
  ! write initial state to file (name: 0.type.dat)
  call write_all()
   last_written = last_written+write_intervall

  if(debuglevel .GE.1) write(*,*) '__________________TIMESTEPPING_____________________________'
  do main_stp= 0,steps
    !if(benchmarking ==1) call cpu_time(bm_step_starttime)
    if(benchmarking ==1) bm_step_starttime=  omp_get_wtime()

    ! STATISTICS WRITING----------------------------------------------------------------------
    !if(benchmarking ==1) call cpu_time(bm_statwrite_starttime)
    if(benchmarking ==1) bm_statwrite_starttime=  omp_get_wtime()
      if(mod(state%step,(measure_every)).EQ.0) then
        call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
        call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
        call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
        call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)
        call write_u_stat()
        call write_E_stat()
        call write_T_stat()
        call write_C_stat()
        call write_sys_stat()
      end if 
    !if(benchmarking ==1) call cpu_time(bm_statwrite_endtime)
    if(benchmarking ==1) bm_statwrite_endtime=  omp_get_wtime()
    ! ----------------------------------------------------------------------------------------


    ! FILE WRITING----------------------------------------------------------------------------
    !if(benchmarking ==1) call cpu_time(bm_filewrite_starttime)
    if(benchmarking ==1) bm_filewrite_starttime=  omp_get_wtime()
    !if(benchmarking ==1) call cpu_time(bm_filewrite_endtime)
    if(benchmarking ==1) bm_filewrite_endtime=  omp_get_wtime()
      if(state%t > last_written) then
        !if(benchmarking ==1) call cpu_time(bm_filewrite_starttime)
        if(benchmarking ==1) bm_filewrite_starttime=  omp_get_wtime()
        call write_all()
        !write(*,*) 'MAXVAL:', maxval(real(state%u%val(:,:,:,1)))
        last_written = last_written+write_intervall
        !if(benchmarking ==1) call cpu_time(bm_filewrite_endtime)
        if(benchmarking ==1) bm_filewrite_endtime=  omp_get_wtime()
      end if 
    ! ----------------------------------------------------------------------------------------

    ! CONSOLE OUTPUT----------------------------------------------------------------------------
  	if(mod(state%step,(steps/1000)).EQ.0) then
        write(*,*) (state%step/(steps/1000)) ,'permille|step:',main_stp, &
    '|t:',state%t,'| dt:',dt,'|shearing:',shearing,'|sheartime:',sheartime,'T_rm',T_rm
                 
    !    call div_tester()
    end if
    ! ----------------------------------------------------------------------------------------



    ! TIMESTEPPING ---------------------------------------------------------------------------
    !if(benchmarking ==1) call cpu_time(bm_timestepping_starttime)
    if(benchmarking ==1) bm_timestepping_starttime=  omp_get_wtime()
      !call RK4_adjust_dt()
      call RK4_step()
      !call euler_step()
      !call div_tester()
      !call ETD2_step()
    !if(benchmarking ==1) call cpu_time(bm_timestepping_endtime)
    if(benchmarking ==1) bm_timestepping_endtime=  omp_get_wtime()
    ! ----------------------------------------------------------------------------------------
    !if(benchmarking ==1) call cpu_time(bm_step_endtime)
    if(benchmarking ==1) bm_step_endtime=  omp_get_wtime()

    !BENCHMARKING------------------------------------
  	!if(mod(state%step-1,(steps/1000)).EQ.0) then
    if(benchmarking ==1) call bm_evaluate(.true.)
    !end if
  end do

  if(debuglevel <= 1) write(*,*) '__________________END OF TIMESTEPPING______________________'
  if(debuglevel <= 1) write(*,*) '__________________END______________________________________'
  call exit_all()
end program guacamole 



    !if(state%t >1.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear =  0.1
    !end if

    !if(state%t >2.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear = 0.05_rp
    !end if

    !if(state%t >3.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear =0.05+ (state%t-3.0*tmax/5.0_rp)/(tmax/(5.0)) *0.05_rp
    !end if

    !if(state%t >4.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear =0.10
    !end if
