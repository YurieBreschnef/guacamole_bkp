module benchmark
use const
use sys_state

real(kind = rp)               :: bm_step_starttime
real(kind = rp)               :: bm_step_endtime
real(kind = rp)               :: bm_step_time

real(kind = rp)               :: bm_statwrite_starttime
real(kind = rp)               :: bm_statwrite_endtime
real(kind = rp)               :: bm_statwrite_time

real(kind = 8)               :: bm_trafo_starttime
real(kind = 8)               :: bm_trafo_endtime
real(kind = 8)               :: bm_trafo_time

real(kind = rp)               :: bm_filewrite_starttime
real(kind = rp)               :: bm_filewrite_endtime
real(kind = rp)               :: bm_filewrite_time

real(kind = rp)               :: bm_timestepping_starttime
real(kind = rp)               :: bm_timestepping_endtime
real(kind = rp)               :: bm_timestepping_time

real(kind = rp)               :: bm_fu_starttime
real(kind = rp)               :: bm_fu_endtime
real(kind = rp)               :: bm_fu_time
real(kind = rp)               :: bm_fu_Nuk_starttime
real(kind = rp)               :: bm_fu_Nuk_endtime
real(kind = rp)               :: bm_fu_Nuk_time
real(kind = rp)               :: bm_fu_buo_starttime
real(kind = rp)               :: bm_fu_buo_endtime
real(kind = rp)               :: bm_fu_buo_time
real(kind = rp)               :: bm_fu_diff_starttime
real(kind = rp)               :: bm_fu_diff_endtime
real(kind = rp)               :: bm_fu_diff_time
real(kind = rp)               :: bm_fu_shear_starttime
real(kind = rp)               :: bm_fu_shear_endtime
real(kind = rp)               :: bm_fu_shear_time

real(kind = rp)               :: bm_ft_starttime
real(kind = rp)               :: bm_ft_endtime
real(kind = rp)               :: bm_ft_time
real(kind = rp)               :: bm_ft_adv_starttime
real(kind = rp)               :: bm_ft_adv_endtime
real(kind = rp)               :: bm_ft_adv_time
real(kind = rp)               :: bm_ft_diff_starttime
real(kind = rp)               :: bm_ft_diff_endtime
real(kind = rp)               :: bm_ft_diff_time
real(kind = rp)               :: bm_ft_strat_starttime
real(kind = rp)               :: bm_ft_strat_endtime
real(kind = rp)               :: bm_ft_strat_time

real(kind = rp)               :: bm_fc_starttime
real(kind = rp)               :: bm_fc_endtime
real(kind = rp)               :: bm_fc_time
real(kind = rp)               :: bm_fc_adv_starttime
real(kind = rp)               :: bm_fc_adv_endtime
real(kind = rp)               :: bm_fc_adv_time
real(kind = rp)               :: bm_fc_diff_starttime
real(kind = rp)               :: bm_fc_diff_endtime
real(kind = rp)               :: bm_fc_diff_time
real(kind = rp)               :: bm_fc_strat_starttime
real(kind = rp)               :: bm_fc_strat_endtime
real(kind = rp)               :: bm_fc_strat_time



real(kind = rp)               :: bm_percent_unaccounted

contains

subroutine bm_evaluate(write_to_console)
  logical,intent(in)               :: write_to_console
  bm_step_time  =         bm_step_endtime-bm_step_starttime
  bm_trafo_time=         bm_trafo_endtime-bm_trafo_starttime
  bm_fu_time    =         bm_fu_endtime-bm_fu_starttime
  bm_fu_Nuk_time    =     bm_fu_Nuk_endtime-bm_fu_Nuk_starttime
  bm_fu_buo_time    =     bm_fu_buo_endtime-bm_fu_buo_starttime
  bm_fu_diff_time    =     bm_fu_diff_endtime-bm_fu_diff_starttime
  bm_fu_shear_time    =     bm_fu_shear_endtime-bm_fu_shear_starttime

  bm_ft_time    =         bm_ft_endtime-bm_ft_starttime
  bm_ft_adv_time    =     bm_ft_adv_endtime-bm_ft_adv_starttime
  bm_ft_diff_time   =     bm_ft_diff_endtime-bm_ft_diff_starttime
  bm_ft_strat_time   =     bm_ft_strat_endtime-bm_ft_strat_starttime

  bm_fc_time    =         bm_fc_endtime-bm_fc_starttime
  bm_fc_adv_time    =     bm_fc_adv_endtime-bm_fc_adv_starttime
  bm_fc_diff_time   =     bm_fc_diff_endtime-bm_fc_diff_starttime
  bm_fc_strat_time   =    bm_fc_strat_endtime-bm_fc_strat_starttime

  bm_timestepping_time =  bm_timestepping_endtime-bm_timestepping_starttime
  bm_statwrite_time =     bm_statwrite_endtime-bm_statwrite_starttime
  bm_filewrite_time =     bm_filewrite_endtime-bm_filewrite_starttime
  bm_timestepping_time=   bm_timestepping_endtime-bm_timestepping_starttime

  if(write_to_console) then
  write(*,*) '______________________________BENCHMARK:_______________________________________'
  write(*,*) 'total step:                :',bm_step_time               ,'sec,'
  write(*,*) '  -statwrite               :  ',bm_statwrite_time          ,'sec,',int(100.0_rp*bm_statwrite_time/bm_step_time),'%'
  write(*,*) '  -filewrite               :  ',bm_filewrite_time,'sec,',int(100.0_rp*bm_filewrite_time/bm_step_time),'%'
  write(*,*) '  -timestepping:           :  ',bm_timestepping_time,'sec,',int(100.0_rp*bm_timestepping_time/bm_step_time),'%'
  write(*,*) '    -function fu:          :    ',bm_fu_time,'sec,',int(100.0_rp*bm_fu_time/bm_step_time),'%'
  write(*,*) '       -function fu_Nuk    :      ',bm_fu_Nuk_time,'sec,',int(100.0_rp*bm_fu_Nuk_time/bm_fu_time),'%'
  write(*,*) '       -function fu_buo    :      ',bm_fu_buo_time,'sec,',int(100.0_rp*bm_fu_buo_time/bm_fu_time),'%'
  write(*,*) '       -function fu_diff   :      ',bm_fu_diff_time,'sec,',int(100.0_rp*bm_fu_diff_time/bm_fu_time),'%'
  write(*,*) '       -function fu_shear  :      ',bm_fu_shear_time,'sec,',int(100.0_rp*bm_fu_shear_time/bm_fu_time),'%'

  write(*,*) '    -function ft:          :    ',bm_ft_time,'sec,',int(100.0_rp*bm_ft_time/bm_step_time),'%'
  write(*,*) '       -function ft_adv    :      ',bm_ft_adv_time,'sec,',int(100.0_rp*bm_ft_adv_time/bm_ft_time),'%'
  write(*,*) '       -function ft_diff   :      ',bm_ft_diff_time,'sec,',int(100.0_rp*bm_ft_diff_time/bm_ft_time),'%'
  write(*,*) '       -function ft_strat  :      ',bm_ft_strat_time,'sec,',int(100.0_rp*bm_ft_strat_time/bm_ft_time),'%'

  write(*,*) '    -function fc:          :    ',bm_fc_time,'sec,',int(100.0_rp*bm_fc_time/bm_step_time),'%'
  write(*,*) '       -function fc_adv    :      ',bm_fc_adv_time,'sec,',int(100.0_rp*bm_fc_adv_time/bm_fc_time),'%'
  write(*,*) '       -function fc_diff   :      ',bm_fc_diff_time,'sec,',int(100.0_rp*bm_fc_diff_time/bm_fc_time),'%'
  write(*,*) '       -function fc_strat  :      ',bm_fc_strat_time,'sec,',int(100.0_rp*bm_fc_strat_time/bm_fc_time),'%'
  write(*,*) ''
  write(*,*) '    -single trafo          :',bm_trafo_time,'sec,',int(100.0_rp*bm_trafo_time/bm_step_time),'%'
  write(*,*) 'percent unnacounted        :',1.0_rp - (bm_statwrite_time+bm_filewrite_time+bm_timestepping_time)/bm_step_time,'%'
  ! note that ft and fc will take the same ammount of computing time
  end if
end subroutine

end module
