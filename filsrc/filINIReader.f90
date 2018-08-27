subroutine readini(name,iflag)
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
	character(len=100) :: name
	integer :: iflag,ifound
	character(len=100) :: iniStr

	if (iflag==1) then
	call setIniFilename(ininam)

	else
	call setIniFilename(trim(name))
	endif

	!********** section data *************

	!**string**
	call getValue(ifound,iflag,'data','s_fileName',iniStr,'data')
	if(ifound==1 .or. iflag==1)  fileName= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'data','s_fileDir',iniStr,'fildata')
	if(ifound==1 .or. iflag==1)  fileDir= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'data','s_tecPlotDir',iniStr,'rcme')
	if(ifound==1 .or. iflag==1)  tecPlotDir= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'data','s_tecPlotFileName',iniStr,'tec.dat')
	if(ifound==1 .or. iflag==1)  tecPlotFileName= trim(adjustl(iniStr))

	!**logical**
	call getValue(ifound,iflag,'data','b_tecPlot3D',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tecPlot3D

	!**string**
	call getValue(ifound,iflag,'data','s_am03ModelDir',iniStr,'am03')
	if(ifound==1 .or. iflag==1)  am03ModelDir= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'data','s_am03Model',iniStr,'marina_data.dat')
	if(ifound==1 .or. iflag==1)  am03Model= trim(adjustl(iniStr))

	!**real**
	call getValue(ifound,iflag,'data','f_am03ModelTime',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) am03ModelTime

	!**integer**
	call getValue(ifound,iflag,'data','i_upscalingFactor',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) upscalingFactor

	!**string**
	call getValue(ifound,iflag,'data','s_tsyModelDir',iniStr,'tsydata')
	if(ifound==1 .or. iflag==1)  tsyModelDir= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'data','s_tsyModel',iniStr,'tsymodel.dat')
	if(ifound==1 .or. iflag==1)  tsyModel= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'data','s_fricDir',iniStr,'fric')
	if(ifound==1 .or. iflag==1)  fricDir= trim(adjustl(iniStr))

	!**logical**
	call getValue(ifound,iflag,'data','b_useGrid',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useGrid

	!**logical**
	call getValue(ifound,iflag,'data','b_buildKubyshData',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) buildKubyshData

	!**logical**
	call getValue(ifound,iflag,'data','b_useAM03ModelGrid',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useAM03ModelGrid

	!**logical**
	call getValue(ifound,iflag,'data','b_useTsyModelGrid',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useTsyModelGrid

	!**logical**
	call getValue(ifound,iflag,'data','b_buildTsyData',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) buildTsyData

	!**logical**
	call getValue(ifound,iflag,'data','b_useAnalyticModelGrid',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useAnalyticModelGrid

	!**logical**
	call getValue(ifound,iflag,'data','b_rcmdata',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) rcmdata

	!**logical**
	call getValue(ifound,iflag,'data','b_write_binary_data',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) write_binary_data

	!**logical**
	call getValue(ifound,iflag,'data','b_write_full_filament',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) write_full_filament

	!**logical**
	call getValue(ifound,iflag,'data','b_write_midpoint',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) write_midpoint

	!**logical**
	call getValue(ifound,iflag,'data','b_write_boundary',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) write_boundary

	!**logical**
	call getValue(ifound,iflag,'data','b_load_external_file',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) load_external_file

	!**logical**
	call getValue(ifound,iflag,'data','b_save_external_file',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) save_external_file

	!**logical**
	call getValue(ifound,iflag,'data','b_save_after_sim',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) save_after_sim

	!********** section main *************

	!**integer**
	call getValue(ifound,iflag,'main','i_method',iniStr,'5')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) method

	!**integer**
	call getValue(ifound,iflag,'main','i_algIn',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) algIn

	!**integer**
	call getValue(ifound,iflag,'main','i_drvIn',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drvIn

	!**integer**
	call getValue(ifound,iflag,'main','i_cnfIn',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) cnfIn

	!**integer**
	call getValue(ifound,iflag,'main','i_dimj',iniStr,'499')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) dimj

	!**real**
	call getValue(ifound,iflag,'main','f_apex',iniStr,'-12.5')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) apex

	!**real**
	call getValue(ifound,iflag,'main','f_boundary',iniStr,'5.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) boundary

	!**real**
	call getValue(ifound,iflag,'main','f_pressureScale',iniStr,'1.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) pressureScale

	!**real**
	call getValue(ifound,iflag,'main','f_tau',iniStr,'0.01')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tau

	!**real**
	call getValue(ifound,iflag,'main','f_tau_min_limit',iniStr,'1.0e-5')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tau_min_limit

	!**real**
	call getValue(ifound,iflag,'main','f_tau_scale',iniStr,'0.1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tau_scale

	!**real**
	call getValue(ifound,iflag,'main','f_time',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) time

	!**real**
	call getValue(ifound,iflag,'main','f_tMax',iniStr,'3600.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tMax

	!**real**
	call getValue(ifound,iflag,'main','f_tInterv',iniStr,'10.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tInterv

	!**integer**
	call getValue(ifound,iflag,'main','i_n',iniStr,'0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) n

	!**integer**
	call getValue(ifound,iflag,'main','i_m',iniStr,'100')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) m

	!**integer**
	call getValue(ifound,iflag,'main','i_nInterv',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) nInterv

	!**integer**
	call getValue(ifound,iflag,'main','i_nOutInterv',iniStr,'1000')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) nOutInterv

	!**integer**
	call getValue(ifound,iflag,'main','i_run',iniStr,'0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) run

	!**logical**
	call getValue(ifound,iflag,'main','b_plotL',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) plotL

	!**logical**
	call getValue(ifound,iflag,'main','b_adaptL',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) adaptL

	!**logical**
	call getValue(ifound,iflag,'main','b_do_parallel',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) do_parallel

	!**integer**
	call getValue(ifound,iflag,'main','i_threads_to_use',iniStr,'0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) threads_to_use

	!**logical**
	call getValue(ifound,iflag,'main','b_bbfrun',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) bbfrun

	!********** section forces *************

	!**logical**
	call getValue(ifound,iflag,'forces','b_use_drag_force',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) use_drag_force

	!**logical**
	call getValue(ifound,iflag,'forces','b_use_driver_force',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) use_driver_force

	!**real**
	call getValue(ifound,iflag,'forces','f_drive_tmin',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drive_tmin

	!**real**
	call getValue(ifound,iflag,'forces','f_drive_tmax',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drive_tmax

	!**real**
	call getValue(ifound,iflag,'forces','f_drag_tmin',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drag_tmin

	!**real**
	call getValue(ifound,iflag,'forces','f_drag_tmax',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drag_tmax

	!**real**
	call getValue(ifound,iflag,'forces','f_drive_coeff',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drive_coeff

	!**real**
	call getValue(ifound,iflag,'forces','f_drag_coeff',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drag_coeff

	!**real**
	call getValue(ifound,iflag,'forces','f_drive_omega_0',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drive_omega_0

	!**real**
	call getValue(ifound,iflag,'forces','f_drive_omega_f',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) drive_omega_f

	!********** section initCond *************

	!**logical**
	call getValue(ifound,iflag,'initCond','b_plotPVG',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) plotPVG

	!**real**
	call getValue(ifound,iflag,'initCond','f_plotPVG_xmin',iniStr,'-14.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) plotPVG_xmin

	!**real**
	call getValue(ifound,iflag,'initCond','f_plotPVG_xmax',iniStr,'-6.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) plotPVG_xmax

	!**integer**
	call getValue(ifound,iflag,'initCond','i_plotPVG_grid',iniStr,'100')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) plotPVG_grid

	!**logical**
	call getValue(ifound,iflag,'initCond','b_useAnalyticDipole',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useAnalyticDipole

	!**logical**
	call getValue(ifound,iflag,'initCond','b_forceModel',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) forceModel

	!**real**
	call getValue(ifound,iflag,'initCond','f_fmsection',iniStr,'0.3')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fmsection

	!**integer**
	call getValue(ifound,iflag,'initCond','i_interpType',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) interpType

	!**logical**
	call getValue(ifound,iflag,'initCond','b_use_dPdK',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) use_dPdK

	!**logical**
	call getValue(ifound,iflag,'initCond','b_dPdK_fixed',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) dPdK_fixed

	!**logical**
	call getValue(ifound,iflag,'initCond','b_useKpk',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useKpk

	!**real**
	call getValue(ifound,iflag,'initCond','f_Kpk',iniStr,'-22.44')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) Kpk

	!**real**
	call getValue(ifound,iflag,'initCond','f_ampli',iniStr,'0.5')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) ampli

	!**logical**
	call getValue(ifound,iflag,'initCond','b_trace_scale_step',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) trace_scale_step

	!**real**
	call getValue(ifound,iflag,'initCond','f_tBack',iniStr,'4.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tBack

	!**real**
	call getValue(ifound,iflag,'initCond','f_dBack',iniStr,'2.80102')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) dBack

	!**logical**
	call getValue(ifound,iflag,'initCond','b_useBackgroundDensity',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) useBackgroundDensity

	!**logical**
	call getValue(ifound,iflag,'initCond','b_use_gallagherDensity',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) use_gallagherDensity

	!**logical**
	call getValue(ifound,iflag,'initCond','b_use_tm2003_density',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) use_tm2003_density

	!**logical**
	call getValue(ifound,iflag,'initCond','b_use_bao_density',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) use_bao_density

	!**logical**
	call getValue(ifound,iflag,'initCond','b_wallBound',iniStr,'1')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) wallBound

	!**logical**
	call getValue(ifound,iflag,'initCond','b_radialDamping',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) radialDamping

	!**logical**
	call getValue(ifound,iflag,'initCond','b_dontRunSim',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) dontRunSim

	!**real**
	call getValue(ifound,iflag,'initCond','f_sigmaP',iniStr,'4.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) sigmaP

	!**integer**
	call getValue(ifound,iflag,'initCond','i_iopen',iniStr,'0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) iopen

	!**real**
	call getValue(ifound,iflag,'initCond','f_fac',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fac

	!********** section sweep *************

	!**logical**
	call getValue(ifound,iflag,'sweep','b_do_sweep',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) do_sweep

	!**string**
	call getValue(ifound,iflag,'sweep','s_sw_var',iniStr,'apex')
	if(ifound==1 .or. iflag==1)  sw_var= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'sweep','s_sw_min',iniStr,'-15.0')
	if(ifound==1 .or. iflag==1)  sw_min= trim(adjustl(iniStr))

	!**string**
	call getValue(ifound,iflag,'sweep','s_sw_max',iniStr,'-11.0')
	if(ifound==1 .or. iflag==1)  sw_max= trim(adjustl(iniStr))

	!**integer**
	call getValue(ifound,iflag,'sweep','i_sw_points',iniStr,'10')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) sw_points

	!**logical**
	call getValue(ifound,iflag,'sweep','b_rebuild_background',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) rebuild_background

	!********** section cw99 *************

	!**real**
	call getValue(ifound,iflag,'cw99','f_A0',iniStr,'-530.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) A0

	!**real**
	call getValue(ifound,iflag,'cw99','f_ht',iniStr,'4.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) ht

	!**real**
	call getValue(ifound,iflag,'cw99','f_alpha',iniStr,'-0.0785')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) alpha

	!**real**
	call getValue(ifound,iflag,'cw99','f_alpha0',iniStr,'-0.0785')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) alpha0

	!********** section tsyganenko *************

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_txmax',iniStr,'5.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) txmax

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_txmin',iniStr,'-35.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) txmin

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_tzmax',iniStr,'10.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tzmax

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_tzmin',iniStr,'-10.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tzmin

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_tyslice',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tyslice

	!**integer**
	call getValue(ifound,iflag,'tsyganenko','i_tgridx',iniStr,'40')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tgridx

	!**integer**
	call getValue(ifound,iflag,'tsyganenko','i_tgridz',iniStr,'40')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tgridz

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_swpress',iniStr,'1.5')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) swpress

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_dst',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) dst

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_byimf',iniStr,'1.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) byimf

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_bzimf',iniStr,'1.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) bzimf

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_tsyg1',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tsyg1

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_tsyg2',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tsyg2

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_tilt',iniStr,'-0.1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) tilt

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_sw_density',iniStr,'6.6')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) sw_density

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_sw_velocity',iniStr,'450.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) sw_velocity

	!**integer**
	call getValue(ifound,iflag,'tsyganenko','i_kpindex',iniStr,'2')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) kpindex

	!**string**
	call getValue(ifound,iflag,'tsyganenko','s_tsyVersion',iniStr,'t89')
	if(ifound==1 .or. iflag==1)  tsyVersion= trim(adjustl(iniStr))

	!**real**
	call getValue(ifound,iflag,'tsyganenko','f_pressureCap',iniStr,'20.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) pressureCap

	!********** section friction *************

	!**integer**
	call getValue(ifound,iflag,'friction','i_initialization_option',iniStr,'0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) initialization_option

	!**logical**
	call getValue(ifound,iflag,'friction','b_NSsymmetric_option',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) NSsymmetric_option

	!**logical**
	call getValue(ifound,iflag,'friction','b_fric_use_inner_bc',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_use_inner_bc

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_vis_mom',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_vis_mom

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_limiter',iniStr,'2')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_limiter

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_pvg_corr_inter',iniStr,'-1000')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_pvg_corr_inter

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_push_rho_opt',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_push_rho_opt

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_sheathrho',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_sheathrho

	!**logical**
	call getValue(ifound,iflag,'friction','b_fric_stopsign',iniStr,'0')
	if(ifound==1 .or. iflag==1)  call intStrToLStr(iniStr,iniStr)
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_stopsign

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_tgamma',iniStr,'0.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_tgamma

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_gammainp',iniStr,'1.66667')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_gammainp

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_rhomin',iniStr,'5.0e-2')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_rhomin

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_pressmin',iniStr,'1.0e-12')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_pressmin

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_cfl',iniStr,'0.25')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_cfl

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_force_norm',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_force_norm

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_max_factor',iniStr,'0.5')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_max_factor

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_min_factor',iniStr,'0.05')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_min_factor

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_algorithm',iniStr,'1')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_algorithm

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_pressure_factor',iniStr,'1.0')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_pressure_factor

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_print_diagnostic',iniStr,'10')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_print_diagnostic

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_min_avg_force',iniStr,'0.00001')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_min_avg_force

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_max_steps',iniStr,'30000')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_max_steps

	!**integer**
	call getValue(ifound,iflag,'friction','i_fric_interval',iniStr,'1000')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_interval

	!**real**
	call getValue(ifound,iflag,'friction','f_fric_alpha',iniStr,'0.05')
	if(ifound==1 .or. iflag==1)  read(iniStr,*) fric_alpha

	
end subroutine readini