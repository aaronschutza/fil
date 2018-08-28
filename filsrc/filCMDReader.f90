subroutine readcmdinput(name,str_len_mid)
	use stateMod
	use stateIndex
	use variables1
	use initCond
	use boundMod
	use parametersBackground
	use filConstants
	use errorMod
	use potential
	use constants
	use wdir
	use dnams
	use fricoptions
	implicit none
integer :: istr
integer :: str_len_mid
integer :: str_len_max=100
character(len=100) :: name
do istr = 2,str_len_mid
select case(name(1:istr))
!**string**
case('s_fileName=','fileName=')
fileName = name(istr+1:str_len_max)
exit
!**string**
case('s_fileDir=','fileDir=')
fileDir = name(istr+1:str_len_max)
exit
!**string**
case('s_tecPlotDir=','tecPlotDir=')
tecPlotDir = name(istr+1:str_len_max)
exit
!**string**
case('s_tecPlotFileName=','tecPlotFileName=')
tecPlotFileName = name(istr+1:str_len_max)
exit
!**logical**
case('b_tecPlot3D=','tecPlot3D=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) tecPlot3D
exit
!**string**
case('s_am03ModelDir=','am03ModelDir=')
am03ModelDir = name(istr+1:str_len_max)
exit
!**string**
case('s_am03Model=','am03Model=')
am03Model = name(istr+1:str_len_max)
exit
!**real**
case('f_am03ModelTime=','am03ModelTime=')
read(name(istr+1:str_len_max),*) am03ModelTime
exit
!**integer**
case('i_upscalingFactor=','upscalingFactor=')
read(name(istr+1:str_len_max),*) upscalingFactor
exit
!**string**
case('s_tsyModelDir=','tsyModelDir=')
tsyModelDir = name(istr+1:str_len_max)
exit
!**string**
case('s_tsyModel=','tsyModel=')
tsyModel = name(istr+1:str_len_max)
exit
!**string**
case('s_fricDir=','fricDir=')
fricDir = name(istr+1:str_len_max)
exit
!**logical**
case('b_useGrid=','useGrid=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useGrid
exit
!**logical**
case('b_buildKubyshData=','buildKubyshData=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) buildKubyshData
exit
!**logical**
case('b_useAM03ModelGrid=','useAM03ModelGrid=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useAM03ModelGrid
exit
!**logical**
case('b_useTsyModelGrid=','useTsyModelGrid=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useTsyModelGrid
exit
!**logical**
case('b_buildTsyData=','buildTsyData=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) buildTsyData
exit
!**logical**
case('b_useAnalyticModelGrid=','useAnalyticModelGrid=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useAnalyticModelGrid
exit
!**logical**
case('b_rcmdata=','rcmdata=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) rcmdata
exit
!**logical**
case('b_write_binary_data=','write_binary_data=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) write_binary_data
exit
!**logical**
case('b_write_full_filament=','write_full_filament=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) write_full_filament
exit
!**logical**
case('b_write_midpoint=','write_midpoint=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) write_midpoint
exit
!**logical**
case('b_write_boundary=','write_boundary=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) write_boundary
exit
!**logical**
case('b_load_external_file=','load_external_file=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) load_external_file
exit
!**logical**
case('b_save_external_file=','save_external_file=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) save_external_file
exit
!**logical**
case('b_save_after_sim=','save_after_sim=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) save_after_sim
exit
!**string**
case('s_external_file=','external_file=')
external_file = name(istr+1:str_len_max)
exit
!**integer**
case('i_method=','method=')
read(name(istr+1:str_len_max),*) method
exit
!**integer**
case('i_algIn=','algIn=')
read(name(istr+1:str_len_max),*) algIn
exit
!**integer**
case('i_drvIn=','drvIn=')
read(name(istr+1:str_len_max),*) drvIn
exit
!**integer**
case('i_cnfIn=','cnfIn=')
read(name(istr+1:str_len_max),*) cnfIn
exit
!**integer**
case('i_dimj=','dimj=')
read(name(istr+1:str_len_max),*) dimj
exit
!**real**
case('f_apex=','apex=')
read(name(istr+1:str_len_max),*) apex
exit
!**real**
case('f_boundary=','boundary=')
read(name(istr+1:str_len_max),*) boundary
exit
!**real**
case('f_pressureScale=','pressureScale=')
read(name(istr+1:str_len_max),*) pressureScale
exit
!**real**
case('f_tau=','tau=')
read(name(istr+1:str_len_max),*) tau
exit
!**real**
case('f_tau_min_limit=','tau_min_limit=')
read(name(istr+1:str_len_max),*) tau_min_limit
exit
!**real**
case('f_tau_scale=','tau_scale=')
read(name(istr+1:str_len_max),*) tau_scale
exit
!**real**
case('f_time=','time=')
read(name(istr+1:str_len_max),*) time
exit
!**real**
case('f_tMax=','tMax=')
read(name(istr+1:str_len_max),*) tMax
exit
!**real**
case('f_tInterv=','tInterv=')
read(name(istr+1:str_len_max),*) tInterv
exit
!**integer**
case('i_n=','n=')
read(name(istr+1:str_len_max),*) n
exit
!**integer**
case('i_m=','m=')
read(name(istr+1:str_len_max),*) m
exit
!**integer**
case('i_nInterv=','nInterv=')
read(name(istr+1:str_len_max),*) nInterv
exit
!**integer**
case('i_nOutInterv=','nOutInterv=')
read(name(istr+1:str_len_max),*) nOutInterv
exit
!**integer**
case('i_run=','run=')
read(name(istr+1:str_len_max),*) run
exit
!**logical**
case('b_plotL=','plotL=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) plotL
exit
!**logical**
case('b_adaptL=','adaptL=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) adaptL
exit
!**logical**
case('b_do_parallel=','do_parallel=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) do_parallel
exit
!**integer**
case('i_threads_to_use=','threads_to_use=')
read(name(istr+1:str_len_max),*) threads_to_use
exit
!**logical**
case('b_bbfrun=','bbfrun=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) bbfrun
exit
!**logical**
case('b_use_drag_force=','use_drag_force=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) use_drag_force
exit
!**logical**
case('b_use_driver_force=','use_driver_force=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) use_driver_force
exit
!**real**
case('f_drive_tmin=','drive_tmin=')
read(name(istr+1:str_len_max),*) drive_tmin
exit
!**real**
case('f_drive_tmax=','drive_tmax=')
read(name(istr+1:str_len_max),*) drive_tmax
exit
!**real**
case('f_drag_tmin=','drag_tmin=')
read(name(istr+1:str_len_max),*) drag_tmin
exit
!**real**
case('f_drag_tmax=','drag_tmax=')
read(name(istr+1:str_len_max),*) drag_tmax
exit
!**real**
case('f_drive_coeff=','drive_coeff=')
read(name(istr+1:str_len_max),*) drive_coeff
exit
!**real**
case('f_drag_coeff=','drag_coeff=')
read(name(istr+1:str_len_max),*) drag_coeff
exit
!**real**
case('f_drive_omega_0=','drive_omega_0=')
read(name(istr+1:str_len_max),*) drive_omega_0
exit
!**real**
case('f_drive_omega_f=','drive_omega_f=')
read(name(istr+1:str_len_max),*) drive_omega_f
exit
!**logical**
case('b_plotPVG=','plotPVG=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) plotPVG
exit
!**real**
case('f_plotPVG_xmin=','plotPVG_xmin=')
read(name(istr+1:str_len_max),*) plotPVG_xmin
exit
!**real**
case('f_plotPVG_xmax=','plotPVG_xmax=')
read(name(istr+1:str_len_max),*) plotPVG_xmax
exit
!**integer**
case('i_plotPVG_grid=','plotPVG_grid=')
read(name(istr+1:str_len_max),*) plotPVG_grid
exit
!**logical**
case('b_useAnalyticDipole=','useAnalyticDipole=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useAnalyticDipole
exit
!**logical**
case('b_forceModel=','forceModel=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) forceModel
exit
!**real**
case('f_fmsection=','fmsection=')
read(name(istr+1:str_len_max),*) fmsection
exit
!**integer**
case('i_interpType=','interpType=')
read(name(istr+1:str_len_max),*) interpType
exit
!**logical**
case('b_use_dPdK=','use_dPdK=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) use_dPdK
exit
!**logical**
case('b_dPdK_fixed=','dPdK_fixed=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) dPdK_fixed
exit
!**logical**
case('b_useKpk=','useKpk=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useKpk
exit
!**real**
case('f_Kpk=','Kpk=')
read(name(istr+1:str_len_max),*) Kpk
exit
!**real**
case('f_ampli=','ampli=')
read(name(istr+1:str_len_max),*) ampli
exit
!**logical**
case('b_trace_scale_step=','trace_scale_step=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) trace_scale_step
exit
!**real**
case('f_tBack=','tBack=')
read(name(istr+1:str_len_max),*) tBack
exit
!**real**
case('f_dBack=','dBack=')
read(name(istr+1:str_len_max),*) dBack
exit
!**logical**
case('b_useBackgroundDensity=','useBackgroundDensity=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) useBackgroundDensity
exit
!**logical**
case('b_use_gallagherDensity=','use_gallagherDensity=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) use_gallagherDensity
exit
!**logical**
case('b_use_tm2003_density=','use_tm2003_density=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) use_tm2003_density
exit
!**logical**
case('b_use_bao_density=','use_bao_density=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) use_bao_density
exit
!**logical**
case('b_wallBound=','wallBound=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) wallBound
exit
!**logical**
case('b_radialDamping=','radialDamping=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) radialDamping
exit
!**logical**
case('b_dontRunSim=','dontRunSim=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) dontRunSim
exit
!**real**
case('f_sigmaP=','sigmaP=')
read(name(istr+1:str_len_max),*) sigmaP
exit
!**integer**
case('i_iopen=','iopen=')
read(name(istr+1:str_len_max),*) iopen
exit
!**real**
case('f_fac=','fac=')
read(name(istr+1:str_len_max),*) fac
exit
!**logical**
case('b_do_sweep=','do_sweep=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) do_sweep
exit
!**string**
case('s_sw_var=','sw_var=')
sw_var = name(istr+1:str_len_max)
exit
!**string**
case('s_sw_min=','sw_min=')
sw_min = name(istr+1:str_len_max)
exit
!**string**
case('s_sw_max=','sw_max=')
sw_max = name(istr+1:str_len_max)
exit
!**integer**
case('i_sw_points=','sw_points=')
read(name(istr+1:str_len_max),*) sw_points
exit
!**logical**
case('b_rebuild_background=','rebuild_background=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) rebuild_background
exit
!**real**
case('f_A0=','A0=')
read(name(istr+1:str_len_max),*) A0
exit
!**real**
case('f_ht=','ht=')
read(name(istr+1:str_len_max),*) ht
exit
!**real**
case('f_alpha=','alpha=')
read(name(istr+1:str_len_max),*) alpha
exit
!**real**
case('f_alpha0=','alpha0=')
read(name(istr+1:str_len_max),*) alpha0
exit
!**real**
case('f_txmax=','txmax=')
read(name(istr+1:str_len_max),*) txmax
exit
!**real**
case('f_txmin=','txmin=')
read(name(istr+1:str_len_max),*) txmin
exit
!**real**
case('f_tzmax=','tzmax=')
read(name(istr+1:str_len_max),*) tzmax
exit
!**real**
case('f_tzmin=','tzmin=')
read(name(istr+1:str_len_max),*) tzmin
exit
!**real**
case('f_tyslice=','tyslice=')
read(name(istr+1:str_len_max),*) tyslice
exit
!**integer**
case('i_tgridx=','tgridx=')
read(name(istr+1:str_len_max),*) tgridx
exit
!**integer**
case('i_tgridz=','tgridz=')
read(name(istr+1:str_len_max),*) tgridz
exit
!**real**
case('f_swpress=','swpress=')
read(name(istr+1:str_len_max),*) swpress
exit
!**real**
case('f_dst=','dst=')
read(name(istr+1:str_len_max),*) dst
exit
!**real**
case('f_byimf=','byimf=')
read(name(istr+1:str_len_max),*) byimf
exit
!**real**
case('f_bzimf=','bzimf=')
read(name(istr+1:str_len_max),*) bzimf
exit
!**real**
case('f_tsyg1=','tsyg1=')
read(name(istr+1:str_len_max),*) tsyg1
exit
!**real**
case('f_tsyg2=','tsyg2=')
read(name(istr+1:str_len_max),*) tsyg2
exit
!**real**
case('f_tilt=','tilt=')
read(name(istr+1:str_len_max),*) tilt
exit
!**real**
case('f_sw_density=','sw_density=')
read(name(istr+1:str_len_max),*) sw_density
exit
!**real**
case('f_sw_velocity=','sw_velocity=')
read(name(istr+1:str_len_max),*) sw_velocity
exit
!**integer**
case('i_kpindex=','kpindex=')
read(name(istr+1:str_len_max),*) kpindex
exit
!**string**
case('s_tsyVersion=','tsyVersion=')
tsyVersion = name(istr+1:str_len_max)
exit
!**real**
case('f_pressureCap=','pressureCap=')
read(name(istr+1:str_len_max),*) pressureCap
exit
!**integer**
case('i_initialization_option=','initialization_option=')
read(name(istr+1:str_len_max),*) initialization_option
exit
!**logical**
case('b_NSsymmetric_option=','NSsymmetric_option=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) NSsymmetric_option
exit
!**logical**
case('b_fric_use_inner_bc=','fric_use_inner_bc=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) fric_use_inner_bc
exit
!**real**
case('f_fric_vis_mom=','fric_vis_mom=')
read(name(istr+1:str_len_max),*) fric_vis_mom
exit
!**integer**
case('i_fric_limiter=','fric_limiter=')
read(name(istr+1:str_len_max),*) fric_limiter
exit
!**integer**
case('i_fric_pvg_corr_inter=','fric_pvg_corr_inter=')
read(name(istr+1:str_len_max),*) fric_pvg_corr_inter
exit
!**integer**
case('i_fric_push_rho_opt=','fric_push_rho_opt=')
read(name(istr+1:str_len_max),*) fric_push_rho_opt
exit
!**integer**
case('i_fric_sheathrho=','fric_sheathrho=')
read(name(istr+1:str_len_max),*) fric_sheathrho
exit
!**logical**
case('b_fric_stopsign=','fric_stopsign=')
call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))
read(name(istr+1:str_len_max),*) fric_stopsign
exit
!**real**
case('f_fric_tgamma=','fric_tgamma=')
read(name(istr+1:str_len_max),*) fric_tgamma
exit
!**real**
case('f_fric_gammainp=','fric_gammainp=')
read(name(istr+1:str_len_max),*) fric_gammainp
exit
!**real**
case('f_fric_rhomin=','fric_rhomin=')
read(name(istr+1:str_len_max),*) fric_rhomin
exit
!**real**
case('f_fric_pressmin=','fric_pressmin=')
read(name(istr+1:str_len_max),*) fric_pressmin
exit
!**real**
case('f_fric_cfl=','fric_cfl=')
read(name(istr+1:str_len_max),*) fric_cfl
exit
!**integer**
case('i_fric_force_norm=','fric_force_norm=')
read(name(istr+1:str_len_max),*) fric_force_norm
exit
!**real**
case('f_fric_max_factor=','fric_max_factor=')
read(name(istr+1:str_len_max),*) fric_max_factor
exit
!**real**
case('f_fric_min_factor=','fric_min_factor=')
read(name(istr+1:str_len_max),*) fric_min_factor
exit
!**integer**
case('i_fric_algorithm=','fric_algorithm=')
read(name(istr+1:str_len_max),*) fric_algorithm
exit
!**real**
case('f_fric_pressure_factor=','fric_pressure_factor=')
read(name(istr+1:str_len_max),*) fric_pressure_factor
exit
!**integer**
case('i_fric_print_diagnostic=','fric_print_diagnostic=')
read(name(istr+1:str_len_max),*) fric_print_diagnostic
exit
!**real**
case('f_fric_min_avg_force=','fric_min_avg_force=')
read(name(istr+1:str_len_max),*) fric_min_avg_force
exit
!**integer**
case('i_fric_max_steps=','fric_max_steps=')
read(name(istr+1:str_len_max),*) fric_max_steps
exit
!**integer**
case('i_fric_interval=','fric_interval=')
read(name(istr+1:str_len_max),*) fric_interval
exit
!**real**
case('f_fric_alpha=','fric_alpha=')
read(name(istr+1:str_len_max),*) fric_alpha
exit

case default

end select
enddo
end subroutine readcmdinput