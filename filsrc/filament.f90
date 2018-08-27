!***********************************************************************************************************************!
! filament main program
! Simulation of thin filament in a background environment
! execute ./fil -h for more information
!***********************************************************************************************************************!
 
program filament
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
use omp_lib
use wdir
use dnams
use rkVectors
implicit none
external :: euler,backEuler,trapezoidal,rk2,rk4
external :: derivs1,algebra1
external :: eulerp,backEulerp,trapezoidalp,rk2p,rk4p
external :: derivs1p,algebra1p
real :: localTime
integer :: i,localN,istr
logical :: notMax
real :: cpuT1,cpuT2,seconds
!**********************************************
integer :: narg,cptArg ! number of arg & counter of arg
character(len=100) :: name,hud
character(len=20) :: str
logical :: outName = .false.
logical :: inputArgs = .false.
logical :: outDir = .false.
logical :: ini_input = .false.
real :: testVal
real :: filK
real :: alphaT
integer :: id, nthreads
real :: maxval,minval
integer :: i_thr
integer :: n_pnt
integer :: str_len_max = 100
integer :: str_len_mid = 40
real :: fil_KinNRG
real,allocatable,dimension(:) :: thr_val
!**********************************************
call cpu_time(cpuT1)
!$ seconds = omp_get_wtime()

call readini(ininam,1)
!**********************************************
!Check if any arguments are found
narg=command_argument_count()
!Loop over the arguments
if(narg>0)then
	!loop across options
	do cptArg=1,narg
		if(outName .or. outDir)then
			outName = .false.
			outDir = .false.
			cycle
		end if
		if(inputArgs)then
			inputArgs = .false.
			cycle
		end if
		if(ini_input)then
			ini_input = .false.
			cycle
		endif
		call get_command_argument(cptArg,name)
		!**********************************************
		select case(adjustl(name))
			case("--help","-h")
				write(*,*)"***Thin filament Code (TFC)***"
				write(*,*)"    Simulation of thin filament in a background environment"
				write(*,*)"    "
				write(*,*)"->Choose data output naming: --output, -o"
				write(*,*)"    Accepts a string e.g. -o myfildata"
				write(*,*)"    "
				write(*,*)"->Choose data output path: --output_directory, -dir"
				write(*,*)"    Accepts a string e.g. -dir myfildatapath"
				write(*,*)"    "
				write(*,*)"->Specify input variable: --input, -i"
				write(*,*)"    Accepts a .ini value in the same format as *.ini"
				write(*,*)"    Any variable listed in default.ini can be changed"
				write(*,*)"    from the command line, with or without the type prefix e.g.:"
				write(*,*)"     ./"//exenam//" -i apex=-12.5  or  ./"//exenam//" -i f_apex=-20.0"
				write(*,*)"    would change the equatorial crossing of the filament."
				write(*,*)"    No spaces allowed in the input statement."
				write(*,*)"    "
				write(*,*)"->Specify ini file to read from:  --inifile -ini"
				write(*,*)"    This will read the specified file after the "//ininam
				write(*,*)"    is loaded so changes in the new file will overwrite"
				write(*,*)"    the old one. e.g. ./"//exenam//" -ini testvals.ini"
				write(*,*)"    will load "//ininam//" first, then any value"
				write(*,*)"    that appears in testvals.ini will be updated."
				write(*,*)"    This is really useful for scheduling a series"
				write(*,*)"    that varies more than one variable"
				write(*,*)"    "
				write(*,*)"->ini usage:"
				write(*,*)"    This executable looks for "//ininam//" and loads"
				write(*,*)"    options from that file.  *.ini files are a loose"
				write(*,*)"    file format for runtime initialization of binary"
				write(*,*)"    files.  There are nice packages for syntax highlighting"
				write(*,*)"    and they support comments (no trailing comments)."
				write(*,*)"    "
				write(*,*)"    default.ini lists all the default values that are"
				write(*,*)"    hardcoded into the executable.  If a value is not"
				write(*,*)"    found in "//ininam//" then the default is loaded"
				write(*,*)"    from the binary file, independent of default.ini."
				write(*,*)"    This allows you to comment out values of "//ininam
				write(*,*)"    or have multiple values saved for different tasks"
				write(*,*)"    in a given file."
				write(*,*)"    "
				write(*,*)"    If the values of default.ini change, i.e. you want to"
				write(*,*)"    change the default initial values, be sure to run"
				write(*,*)"    make runiniscript so that the changes are written"
				write(*,*)"    to the fil source files. "
				write(*,*)"    !!! DO NOT REMOVE ENTRIES FROM default.ini !!!"
				write(*,*)"    This will interfere with make runiniscript."
				write(*,*)"    Any removed entries will not be initialized properly"
				write(*,*)"    if they are still being used by the program."
				write(*,*)"    "
				stop
			case("--output","-o")
				outName = .true.
			case("--output_directory","-dir")
				outDir = .true.
			case("--input","-i")
				inputArgs = .true.
			case("--inifile","-ini")
				ini_input = .true.
			case default
				write(*,*)"Option ",adjustl(name),"unknown"
		end select
		!**********************************************
		if (outName .and. cptArg<narg) then
			call get_command_argument(cptArg+1,name)
			fileName = adjustl(name)
		endif

		if (outDir .and. cptArg<narg) then
			call get_command_argument(cptArg+1,name)
			fileDir = adjustl(name)
		endif

		if (inputArgs .and. cptArg<narg) then
			call get_command_argument(cptArg+1,name)
			name = adjustl(name)
			!**********************************
			call readcmdinput(name,str_len_mid)
			!**********************************
		end if

		if (ini_input .and. cptArg<narg) then
			call get_command_argument(cptArg+1,name)
			name = adjustl(name)
			!**********************************
			call readini(name,0)
			!**********************************
		end if
	end do
end if


!**********************************************
if (do_sweep) then
	n_pnt = sw_points
	allocate(thr_val(n_pnt))
	read(sw_min,*) minval
	read(sw_max,*) maxval
	if (n_pnt > 1) then
		do i = 1,n_pnt
					 thr_val(i) = minval+(maxval-minval)*((i-1.0)/(n_pnt-1.0))
		enddo
	endif
else
	n_pnt = 1
endif
!**********************************************




if (mod(dimj,2)<1) then
	dimj = dimj+1
endif
!**********************************************
allocate(state(dimi,dimj))
allocate(s_0(dimi,dimj))
allocate(stateTemp(dimi,dimj))
allocate(newState(dimi,dimj))
allocate(oldState(dimi,dimj))
allocate(recallPos(dimj,2))
allocate(s1(dimi,dimj));allocate(k1(dimi,dimj))
allocate(s2(dimi,dimj));allocate(k2(dimi,dimj))
allocate(s3(dimi,dimj));allocate(k3(dimi,dimj))
allocate(sa(dimi,dimj));allocate(k4(dimi,dimj))

!leave this untouched for no openMP, single thread
nthreads = 1
id = 0


!$ if (do_parallel) then
!$ 	if (threads_to_use>0) then
!$ 		call omp_set_num_threads(threads_to_use)
!$ 	endif
!$ else
!$ 	call omp_set_num_threads(1)
!$ endif


!$OMP parallel default(shared) private(id)
!$ 		id = omp_get_thread_num()
!$OMP barrier
!$OMP master
!$ 		nthreads = omp_get_num_threads()
!$ 		write(*,*) 'Thread Count = ', nthreads
!$ 		write(*,*) 'Master thread:',id 
!$OMP end master


!$OMP master
allocate(thr_domain(nthreads,2))
do i = 1,nthreads
	if ( i == 1 ) then
		i_thr = 1
		thr_domain(i,1) = i_thr
		i_thr = i_thr+dimj/nthreads-1
		if (mod(dimj,nthreads) /= 0) i_thr = i_thr+1
		thr_domain(i,2) = i_thr
	elseif ( i<=mod(dimj,nthreads) ) then
		i_thr = i_thr+1
		thr_domain(i,1) = i_thr
		i_thr = i_thr+dimj/nthreads
		thr_domain(i,2) = i_thr
	else
		i_thr = i_thr+1
		thr_domain(i,1) = i_thr
		i_thr = i_thr+dimj/nthreads-1
		thr_domain(i,2) = i_thr
	endif
enddo

call system('mkdir -p ' // adjustl(trim(fricDir)))
call system('mkdir -p ' // adjustl(trim(tsyModelDir)))

!$OMP end master

!$OMP barrier

do i_thr = 1,n_pnt !*********begin sweep sequence ****

	!$OMP master
	if(do_sweep .and. n_pnt>1)then
		write(name,*) thr_val(i_thr)
		name = adjustl(trim(sw_var)//'='//trim(adjustl(name)))
		write(*,*) name
		call readcmdinput(name,str_len_mid)
		write (name,'(I3.3)') i_thr
		call system('mkdir -p ' // adjustl(trim(fileDir)))
		datapath = trim(fileDir)//'/'//trim(fileName)//'.'//trim(name)
		write(*,*) 'writing data to:'
		write(*,*) datapath
	else
		call system('mkdir -p ' // adjustl(trim(fileDir)))
		datapath = trim(fileDir)//'/'//trim(fileName)
		write(*,*) 'writing data to:'
		write(*,*) datapath
	endif

	!**********************************************
	localN = 0
	localTime = 0.0
	!**********************************************
	recallPos = -1
	time = 0.0
	n = 0
	if (cnfIn == 1) then
		call initialConditions1
	else
		write(6,*) 'Error: invalid cnfIn'
		stop
	end if
	!**********************************************
	if (write_binary_data) then
		open(unit=111, file=trim(datapath)//".bin",form="unformatted", status="replace",access="stream")
		close(111)
		open(unit=111, file=trim(datapath)//".bin",form="unformatted",position="append", status="old",access="stream")
		write(111) dimi
		write(111) dimj
		write(111) dimjp
		write(111) time
		write(111) n
		write(111) state
		if (dimjp /= 0) then
			write(111) statePlas
		end if
		close(111)
	endif

	if(write_full_filament) then
		open(unit=111, file=trim(datapath), status="replace")
		close(111)
		open(unit=111, file=trim(datapath),position="append", status="old")
		write(111,*) dimi
		write(111,*) dimj + dimjp
		write(111,*) time
		write(111,*) n
		do i = 1,dimj
			write(111,*) state(:,i)
		end do
		if (dimjp /= 0) then
			do i = 1,dimjp
				write(111,*) statePlas(:,i)
			end do
		end if
		close(111)
	endif

	if(write_full_filament) then
		open(unit=111, file=trim(datapath), status="replace")
		close(111)
		open(unit=111, file=trim(datapath),position="append", status="old")
		write(111,*) dimi
		write(111,*) dimj + dimjp
		write(111,*) time
		write(111,*) n
		do i = 1,dimj
			write(111,*) state(:,i)
		end do
		if (dimjp /= 0) then
			do i = 1,dimjp
				write(111,*) statePlas(:,i)
			end do
		end if
		close(111)
	endif

	if (write_midpoint) then
		open(unit=111, file=trim(datapath)//'.mid', status="replace")
		close(111)
		open(unit=111, file=trim(datapath)//'.mid',position="append", status="old")
		write(111,*) time,state(x,dimj/2+1),state(z,dimj/2+1),state(vx,dimj/2+1),state(vz,dimj/2+1), &
		& state(d,dimj/2+1),state(p,dimj/2+1),state(b,dimj/2+1),filK(state)
		close(111)
	endif

	if (write_boundary) then
		open(unit=111, file=trim(datapath)//'.bnd', status="replace")
		close(111)
		open(unit=111, file=trim(datapath)//'.bnd',position="append", status="old")
		write(111,*) time,state(x,1),state(z,1),state(x,dimj),state(z,dimj)
		close(111)
	endif

	if (plotL) then
		notMax = n <= m
	else
		notMax = time < tMax
	end if
	!$OMP end master

	!$OMP barrier

	do while (notMax)
		!$OMP master
		if (mod(n,nOutInterv) == 0) then
			write(6,*) 'n =', n,'t =',time,'tau =',tau
			write(6,*) 'xe=',state(x,dimj/2+1),'KE=',fil_KinNRG(state)
			!call testState
		end if
		if (tau < tau_min_limit) then
			write(*,*) 'time step too small'
			stop
		end if
		!$OMP end master

		!$OMP barrier
		select case(method)
			case(1)
				!method = euler
				if (adaptL) then 
					call adaptive(euler,derivs1,algebra1,eulerp,derivs1p,algebra1p,id)
				else
					call fixedtau(euler,derivs1,algebra1,eulerp,derivs1p,algebra1p,id)
				end if
			case(2)
				!methodP = backEuler
				if (adaptL) then 
					call adaptive(backEuler,derivs1,algebra1,backEulerp,derivs1p,algebra1p,id)
				else
					call fixedtau(backEuler,derivs1,algebra1,backEulerp,derivs1p,algebra1p,id)
				end if
			case(3)
				!methodP = trapezoidal
				if (adaptL) then 
					call adaptive(trapezoidal,derivs1,algebra1,trapezoidalp,derivs1p,algebra1p,id)
				else
					call fixedtau(trapezoidal,derivs1,algebra1,trapezoidalp,derivs1p,algebra1p,id)
				end if
			case(4)
				!methodP = rk2
				if (adaptL) then 
					call adaptive(rk2,derivs1,algebra1,rk2p,derivs1p,algebra1p,id)
				else
					call fixedtau(rk2,derivs1,algebra1,rk2p,derivs1p,algebra1p,id)
				end if
			case(5)
				!methodP = rk4
				if (adaptL) then 
					call adaptive(rk4,derivs1,algebra1,rk4p,derivs1p,algebra1p,id)
				else
					call fixedtau(rk4,derivs1,algebra1,rk4p,derivs1p,algebra1p,id)
				end if
			case default
				write(6,*) 'Error: invalid method'
				stop
		end select
		!$OMP barrier

		!$OMP master
		oldState = state
		state = newState
		if(dimjp/=0)then
			oldStatePlas = statePlas
			statePlas = newStatePlas
		endif
		!if (.not. errorFlag2) then
		!call reconnect2
		!end if
		localTime = localTime + tau
		localN = localN + 1
		time = time + tau
		n = n + 1
		!alpha = alphaT(time)
		if(bbfrun)then
			if (state(vx,dimj/2+1)<0.0 .or. abs(state(x,dimj/2+1))<boundary .and. time>tInterv) then
				notMax = .false. !exit when tailward motion starts
				localTime = tInterv+1.0 !write last data
			endif
		endif

		if (plotL) then
			if (localN >= nInterv) then
				localN = 0
				if (write_binary_data) then
					open(unit=111, file=trim(datapath)//".bin",form="unformatted",position="append", status="old")
					write(111,*) dimi
					write(111,*) dimj
					write(111,*) dimjp
					write(111,*) time
					write(111,*) n
					write(111,*) state
					if (dimjp /= 0) then
						write(111,*) statePlas
					end if
					close(111)
				endif

				if(write_full_filament) then
					open(unit=111, file=trim(datapath),position="append", status="old")
					write(111,*) dimi
					write(111,*) dimj + dimjp
					write(111,*) time
					write(111,*) n
					do i = 1,dimj
						write(111,*) state(:,i)
					end do
					if (dimjp /= 0) then
						do i = 1,dimjp
							write(111,*) statePlas(:,i)
						end do
					end if
					close(111)
				endif

				if (write_midpoint) then
					open(unit=111, file=trim(datapath)//'.mid',position="append", status="old")
					write(111,*) time,state(x,dimj/2+1),state(z,dimj/2+1),state(vx,dimj/2+1),state(vz,dimj/2+1)
					close(111)
				endif

				if (write_boundary) then
					open(unit=111, file=trim(datapath)//'.bnd',position="append", status="old")
					write(111,*) time,state(x,1),state(z,1),state(x,dimj),state(z,dimj)
					close(111)
				endif
			end if
		else
			if (localTime > tInterv) then
				localTime = 0.0
				if (write_binary_data) then
					open(unit=111, file=trim(datapath)//".bin",form="unformatted",&
					position="append", status="old",access="stream")
					write(111) dimi
					write(111) dimj
					write(111) dimjp
					write(111) time
					write(111) n
					write(111) state
					if (dimjp /= 0) then
						write(111) statePlas
					end if
					close(111)
				endif

				if(write_full_filament) then
					open(unit=111, file=trim(datapath),position="append", status="old")
					write(111,*) dimi
					write(111,*) dimj + dimjp
					write(111,*) time
					write(111,*) n
					do i = 1,dimj
						write(111,*) state(:,i)
					end do
					if (dimjp /= 0) then
						do i = 1,dimjp
							write(111,*) statePlas(:,i)
						end do
					end if
					close(111)
				endif

				if (write_midpoint) then
					open(unit=111, file=trim(datapath)//'.mid',position="append", status="old")
					write(111,*) time,state(x,dimj/2+1),state(z,dimj/2+1),state(vx,dimj/2+1),state(vz,dimj/2+1), &
					& state(d,dimj/2+1),state(p,dimj/2+1),state(b,dimj/2+1),filK(state)
					close(111)
				endif

				if (write_boundary) then
					open(unit=111, file=trim(datapath)//'.bnd',position="append", status="old")
					write(111,*) time,state(x,1),state(z,1),state(x,dimj),state(z,dimj)
					close(111)
				endif
			end if
		end if

		if (plotL .and. notMax) then
			notMax = n <= m
		endif
		if (.not. plotL .and. notMax) then
			notMax = time < tMax
		endif
		!$OMP end master
		!$OMP barrier	 
	enddo
enddo!*********end thr_val sequence ****
!$OMP end parallel
run = 0
call cpu_time(cpuT2)
write(*,*) 'Elapsed CPU time =',cpuT2-cpuT1
!$ seconds = omp_get_wtime() - seconds
!$ write(*,*) 'Wallclock time =',seconds


end program
