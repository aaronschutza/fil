! Routines for calculating the position of the contours of a general bacgkround potential

subroutine contour(x,y,z,ni,nj,c,xin,yin) ! find contour of potential grid for initial conditions
	use cont
	implicit none
	real :: c,x0,y0,x1,y1,xin,yin,z0,z1,t
	integer :: i,j,k,l,m,n,dimc,dir,ni,nj
	logical :: l0,l1,l2,l3
	real,dimension(ni,nj) :: z
	real,dimension(ni) :: x
	real,dimension(nj) :: y
	integer,dimension(ni-1,nj-1) :: v
	
	
	n=0
	m=0
	
	write(*,*) xin,yin
	v = 0
	do i = 1,ni-1
		if ( x(i) < xin .and. x(i+1) >= xin  ) then
			n = i
			write(*,*) n
		end if
	end do
	do j = 1,nj-1
		if ( y(j) < yin .and. y(j+1) >= yin  ) then
			m = j
			write(*,*) m
		end if
	end do
	dimc = 10000
	allocate(coni(dimc))
	allocate(conj(dimc))
	dir = 1
	k = 1
	coni(k) = xin
	conj(k) = yin
	do k = 2,dimc
	  write(6,*) 'k=',k,'n=',n,'m=',m
	  if (v(n,m) /= 0) then
		write(6,*) 'cell visited'
		call truncCon(k,dimc)
		return
	  end if
	  v(n,m) = 1
	  l0 =  z(n+1,m) > c .and. z(n+1,m+1) > c &
	  &.or. z(n+1,m) < c .and. z(n+1,m+1) < c
	  l1 =  z(n,m+1) > c .and. z(n+1,m+1) > c &
	  &.or. z(n,m+1) < c .and. z(n+1,m+1) < c
	  l2 =    z(n,m) > c .and. z(n,m+1) > c &
	  &.or.   z(n,m) < c .and. z(n,m+1) < c
	  l3 =    z(n,m) > c .and. z(n+1,m) > c &
	  &.or.   z(n,m) < c .and. z(n+1,m) < c
	  if (l0 .and. l1 .and. l2 .and. l3) then
		write(6,*) 'error, no contour found'
		return
	  end if
	  select case (dir)
	    case (2)
		if (.not. l1) then
			z0 = z(n,m+1)
			z1 = z(n+1,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m+1)
			x1 = x(n+1)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 1
		else if (.not. l2) then
			z0 = z(n,m)
			z1 = z(n,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m)
			x1 = x(n)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 2
		else if (.not. l3) then
			z0 = z(n,m)
			z1 = z(n+1,m)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m)
			x1 = x(n+1)
			y1 = y(m)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 3
		else
			write(6,*) 'error0',k
			return
		end if
	    case (3)
		if (.not. l0) then
			z0 = z(n+1,m)
			z1 = z(n+1,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n+1)
			y0 = y(m)
			x1 = x(n+1)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 0
		else if (.not. l2) then
			z0 = z(n,m)
			z1 = z(n,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m)
			x1 = x(n)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 2
		else if (.not. l3) then
			z0 = z(n,m)
			z1 = z(n+1,m)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m)
			x1 = x(n+1)
			y1 = y(m)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 3
		else
			write(6,*) 'error1',k
			return
		end if
	    case (0)
		if (.not. l0) then
			z0 = z(n+1,m)
			z1 = z(n+1,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n+1)
			y0 = y(m)
			x1 = x(n+1)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 0
		else if (.not. l1) then
			z0 = z(n,m+1)
			z1 = z(n+1,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m+1)
			x1 = x(n+1)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 1
		else if (.not. l3) then
			z0 = z(n,m)
			z1 = z(n+1,m)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m)
			x1 = x(n+1)
			y1 = y(m)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 3
		else
			write(6,*) 'error'
			return
		end if
	    case (1)
		if (.not. l0) then
			z0 = z(n+1,m)
			z1 = z(n+1,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n+1)
			y0 = y(m)
			x1 = x(n+1)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 0
		else if (.not. l1) then
			z0 = z(n,m+1)
			z1 = z(n+1,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m+1)
			x1 = x(n+1)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 1
		else if (.not. l2) then
			z0 = z(n,m)
			z1 = z(n,m+1)
			t = (c-z0)/(z1-z0)
			x0 = x(n)
			y0 = y(m)
			x1 = x(n)
			y1 = y(m+1)
			coni(k) = x0+t*(x1-x0)
			conj(k) = y0+t*(y1-y0)
			dir = 2
		else
			write(6,*) 'error'
			return
		end if
	    case default
		write(6,*) 'error, dir \= 0,1,2,3'
	  end select
	  select case(dir)
	    case (0)
		n = n+1
	    case (1)
		m = m+1
	    case (2)
		n = n-1
	    case (3)
		m = m-1
	    case default
		write(6,*) 'error'
	  end select
	  if (n>=ni-1 .or. n <= 0) then
		write(6,*) 'error, n out of bounds'
		return
	  end if
	  if (m>=nj-1 .or. m <= 0) then
		write(6,*) 'error, m out of bounds'
		return
	  end if
	  if ( sqrt(x(n)**2+y(m)**2) < 0.75 ) then
		write(6,*) '0.75 Re reached'
		call truncCon(k,dimc)
		return
	  end if
	end do
end subroutine contour

subroutine truncCon(k,dimc) ! gets rid of empty cells from coni and conj
	use cont
	implicit none
	integer :: i,k
	integer :: dimc
	real,dimension(dimc) :: coniTemp
	real,dimension(dimc) :: conjTemp
	coniTemp = coni
	conjTemp = conj
	deallocate (coni)
	deallocate (conj)
	allocate (coni(k))
	allocate (conj(k))
	Do i = 1,k
		coni(i) = coniTemp(i)
	end do
	do i = 1,k
		conj(i) = conjTemp(i)
	end do
end subroutine truncCon
