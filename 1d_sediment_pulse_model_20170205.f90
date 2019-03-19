program model_1d
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
!                                                                              !
!                  ONE-DIMENSIONAL HYDRO-MORPHODYNAMIC MODEL                   !
!                                                                              !
!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Written by: Jacob A. Morgan
!			  Colorado State University
!			  Fort Collins, Colorado
!
! Last Edit: 05 February 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 character (LEN=13) :: fname
 character (LEN=2) :: RUN
 character (LEN=4) :: TIMEM
 character (LEN=3) :: ii
integer :: &
	i, &! Downstream spatial index
	j, &! Grain size distiribution size class index
	k, &! Vertical grain sorting index
	f, &! Index for output file units
	o, &
	tt
integer, parameter :: &
	M = 229, &! Number of downstream nodes
	nw = 6, &! Number of wavelenths of width variation
	npp = 29, &! Number of grain size classes in CDF distributions
	np = npp-1, &! Number of grain size classes in PDF distributions
	N = 200 ! Maximum number of vertical nodes for storing stratigraphy
integer, dimension(1:M) :: &
	Nn ! Number of vertical nodes with stored stratigraphy data
real, parameter :: &
	! Physical Constants
	pi = 3.14159265359, &! Constant, pi
	g = 9.81, &! Gravitational acceleration (m/s^2)
	rho = 1000., &! Density of water (kg/m^3)
	Rr = 1.65, &! Submerged specific gravity of sediment
	lps = 0.4, &! Sediment porosity
	! Geometric Constants
	Ac = 0.1676, &! Dimensionless amplitude of width variations
	Lc = 1.08, &! Dimensional wavelength of width variations (m)
	enL = 1.58, &! Length of flume entrance reach (m)
	exL = 1.08, &! Length of flume exit reach (m)
	enB = 0.108, &! Half-width of flume entrance reach (m)
	exB = 0.108, &! Half-width of flume exit reach (m)
	B0 = 0.0925, &! Mean half-width of width variations (m)
	Sb = 0.007, &! Initial constant bed slope
	etad = 0.048, &! Fixed elevation of downstream-most node (m)
	L = enL+nw*Lc+exL, &! Total length of flume (m)
	Ls = 0.0025, &! Thickness of substrate layers (m)
	! Numerical Constants
	nk = 2., &! Factor multiplied by surface D90 for roughness height, ks
	nactive = 2., &! Factor multiplied by surface D90 for active layer thickness
	ar = 8.1, &! Coefficient in Manning-Strickler resistance relation
	au = 0.75, &! Upwinding coefficient for load derivatives in Exner equation
	atrans = 0.7, &! Coefficient for material transferred to substrate for agg.
	dx = L/(M-1), &! Spatial step length (m)
	dt = 0.1, &! Time step increment (s)
	T = (13.+23.6+4.5+28.6+711./60.)*60.*60., &! Total run time (s)
	! Hydraulic and Morphodynamic Constants
	Q = 0.00091, &! Water discharge (m^3/s)
	Qbf = 150./1000./2650./60. ! Sediment feed rate (m^3/s)
real, parameter, dimension(1:npp) :: &
	di = (/ 8., 6.73, 5.66, 4.76, 4., 3.36, &
		2.83, 2.38, 2., 1.68, 1.41, 1.19, 1., 0.841, 0.707, 0.595, 0.5, 0.420, &
		0.354, 0.297, 0.25, 0.210, 0.177, 0.149, 0.125, 0.105, 0.0884, 0.0743, &
		0.0625 /), &! CDF grain size classes (mm)
	cdf_surfi = (/ 100., 99.96, 99.93, 99.9, 99.8, 96.865, 93.93, 90.995, &
		88.06, 80.51, 72.96, 65.41, 57.86, 49.27, 42.95, 34.675, 26.4, 20.683, &
		14.967, 9.25, 5.17, 4.045, 2.92, 2.17, 1.42, 0.9467, 0.4733, 0., 0. &
		/), &! CDF percent finer of initial surface sediment
	cdf_subi = (/ 100., 99.96, 99.93, 99.9, 99.8, 96.865, 93.93, 90.995, &
		88.06, 80.51, 72.96, 65.41, 57.86, 49.27, 42.95, 34.675, 26.4, 20.683, &
		14.967, 9.25, 5.17, 4.045, 2.92, 2.17, 1.42, 0.9467, 0.4733, 0., 0. &
		/), &! CDF percent finer of initial substrate sediment
	cdf_feed = (/ 100., 99.96, 99.93, 99.9, 99.8, 96.865, 93.93, 90.995, &
		88.06, 80.51, 72.96, 65.41, 57.86, 49.27, 42.95, 34.675, 26.4, 20.683, &
		14.967, 9.25, 5.17, 4.045, 2.92, 2.17, 1.42, 0.9467, 0.4733, 0., 0. &
		/) ! CDF percent finer of sediment feed
real, parameter, dimension(1:np) :: &
	ds = (/ (sqrt(di(j)*di(j-1)), j = 2, npp) /), &! PDF grain size classes (mm)
	psi = log10(ds)/log10(2.) ! PDF phi scale grain size classes
real, dimension(1:np) :: pdf_surfi, pdf_subi, pdf_feed
real, dimension(1:M) :: &
	x, &! Downstream node positions (m)
	Bx, &! Local channel half-width (m)
	B, &! Local channel width (m)
	eta, &! Bed elevation at each surface node (m)
	qw, &! Volumetric discharge per unit width (m^2/s)
	dgsurf, &! GSD geometric mean at each surface node (mm)
	sigsurf, &! GSD geometric standard deviation at each surface node
	d90surf, &! D90 at each surface node (mm)
	d50surf, &! D50 at each surface node (mm)
	H, &! Flow depth (m)
	tausg, &! Boundary Shields parameter referenced to GSD geometric mean
	qbT, &! Volumetric bedload transport rate (m^2/s)
	La, &! Thickness of active layer (m)
	Laold ! Thickness of active layer at previous time step (m)
real, dimension(M,npp) :: cdf_surf
real, dimension(M,np) :: pdf_surf, pdf_load
real, dimension(M,N) :: eta_strat, d50strat, d90strat, dgstrat, sigstrat
real, dimension(M,N,np) :: pdf_strat
real, dimension(M,N,npp) :: cdf_strat
real :: dgsurfi, dgsubi, sigsurfi, sigsubi, d50surfi, d50subi, d90surfi, &
	d90subi, time, qbTf, start, finish
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Set up Flume Geometry
!------------------------------------------------------------------------------!
do i = 1, M
	x(i) = (i-1)*dx
	if (x(i) <= enL) then
		Bx(i) = enB
	else if (x(i) >= (enL+nw*Lc)) then
		Bx(i) = exB
	else
		Bx(i) = B0*(1+Ac*cos(2*pi*(x(i)-enL)/Lc))
	end if
end do
B = 2*Bx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Set up Boundary Conditions
!------------------------------------------------------------------------------!
qw = Q/B
	call cdf2pdf(npp, cdf_feed, pdf_feed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Set up Initial Conditions
!------------------------------------------------------------------------------!
do i = 1, M
	eta(i) = Sb*(L-x(i))+etad
end do

!open(unit=1,file="geometry.csv",form="formatted",status="replace")
!write(1,'(3(A10,:,","))') 'x(m)', 'B(m)', 'eta(m)'
!do i = 1, M
!	write(1,'(3(F10.4,:,","))') x(i), B(i), eta(i)
!end do
! close(unit=1)

! Initial grain size distribution statistics
	! Grain size distribution PDFs
	call cdf2pdf(npp, cdf_surfi, pdf_surfi)
	call cdf2pdf(npp, cdf_subi, pdf_subi)
	! Grain size distribution geometric means
	call geomean(np, psi, pdf_surfi, dgsurfi, sigsurfi)
	call geomean(np, psi, pdf_subi, dgsubi, sigsubi)
	! Grain size distribution D90
	call finer(npp, di, cdf_surfi, 90., d90surfi)
	call finer(npp, di, cdf_surfi, 90., d90subi)
	! Grain size distribution D50
	call finer(npp, di, cdf_surfi, 50., d50surfi)
	call finer(npp, di, cdf_surfi, 50., d50subi)
! Assign initial grain size distributions and statistics to all surface nodes
do i = 1, M
	cdf_surf(i,:) = cdf_surfi
	pdf_surf(i,:) = pdf_surfi
	dgsurf(i) = dgsurfi
	sigsurf(i) = sigsurfi
	d90surf(i) = d90surfi
	d50surf(i) = d50surfi
end do
! Assign initial grain size distributions and statistics to all substrate nodes
do i = 1, M
	Nn(i) = nint((eta(i)-(nactive*d90surf(i)/1000.))/Ls)+ 2
	do k = 1, N
		if (k == 1) then
			eta_strat(i,k) = 0.
		else if (k < Nn(i)) then
			eta_strat(i,k) = eta_strat(i,k-1)+Ls
		else if (k == Nn(i)) then
			eta_strat(i,k) = eta(i)-(nactive*d90surf(i)/1000.)
		else if (k > Nn(i)) then
			eta_strat(i,k) = 0./o
		end if
		do j = 1, np
			if (k > Nn(i)) then
				pdf_strat(i,k,j) = 0./o
			else
				pdf_strat(i,k,j) = pdf_subi(j)
			end if
		end do
		if (eta_strat(i,k) > (eta(i)-nactive*d90surf(i)/1000.)) then
			eta_strat(i,k) = 0./o
		end if
	end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Main model
!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Laold = 0.0
time = 0.0
f = 2
tt = 0
 call cpu_time(start)
do while (time <= T)
	time = tt*dt
	if (time <= 13*60*60) then
		RUN = "V1"
		qbTf = Qbf/B(1)
	elseif (time <= (13+23.6)*60*60) then
		RUN = "V2"
		qbTf = 0
	elseif (time <= (13+23.6+4.5)*60*60) then
		RUN = "V3"
		qbTf = Qbf/B(1)
	elseif (time <= (13+23.6+4.5+28.6)*60*60) then
		RUN = "V4"
		qbTf = 0
	else
		RUN = "V5"
		if (time <= (13+23.6+4.5+28.6+27./60.)*60*60) then
			qbTf = 4.*Qbf/B(1)
			call cdf2pdf(npp, (/ 100., 100., 100., 100., 100., 100., 100., &
				100., 100., 100., 100., 100., 100., 0., 0., 0., 0., 0., 0., &
				0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /), pdf_feed)
		else
			qbTf = 0
		end if
	end if
	! Print current time step to disply ---------------------------------------!
	if (time == 0.0 .or. (time)/int(time) == 1.0) then
		!print*, "time = ", nint(time), "s"
		write(*,'("time = ",I6.6," seconds")') int(time)
	end if
	! Calculate backwater and bedload -----------------------------------------!
	call backwater(qw, eta, dx, M, g, Rr, ar, nk, d90surf, dgsurf, H, tausg, &
		time)
	call bedload(tausg,dgsurf,pdf_surf,ds,M,np,Rr,g,qbT,pdf_load,time)
	! Write output to file ----------------------------------------------------!
!	if (time == 0 .or. &
!			(time/60.0)/int((time/60.0)) == 1 .or. &
!			time == (13*60*60) .or. &
!			time == ((13+23.6)*60*60) .or. &
!			time == ((13+23.6+4.5)*60*60) .or. &
!			time == ((13+23.6+4.5+28.6)*60*60) .or. &
!			time == ((13+23.6+4.5+28.6+27./60.)*60*60) .or. &
!			time == ((13+23.6+4.5+28.6+711./60.)*60*60)) then
	write(TIMEM,"(I4.4)") int(time/60.0)
	if (time == 0 .or. (time/60.0)/int(time/60.0) == 1.0) then
		!
		do i = 1, M
			do k = 1, Nn(i)
				call pdf2cdf(npp, pdf_strat(i,k,:), cdf_strat(i,k,:))
			end do
			call pdf2cdf(npp, pdf_surf(i,:), cdf_surf(i,:))
		end do
		!
		fname = "run_"//TIMEM//"m.csv"
		open(unit=f,file="OUTPUT/"//fname,form="formatted",&
			status="replace")
		write(f,'("#","time = ",A4," minutes")') TIMEM
		write(f,'("#",36(A10,:,","))') "x(m)", "z(m)", "B(m)", "eta(m)", &
			"H(m)", "tausg(-)", "QbT(g/s)", "8.00 mm", "6.73 mm", "5.66 mm", &
			"4.76 mm", "4.00 mm", "3.36 mm", "2.83 mm", "2.38 mm", "2.00 mm", &
			"1.68 mm", "1.41 mm", "1.19 mm", "1.00 mm", "0.84 mm", "0.71 mm", &
			"0.59 mm", "0.50 mm", "0.42 mm", "0.35 mm", "0.30 mm", "0.25 mm", &
			"0.21 mm", "0.18 mm", "0.15 mm", "0.13 mm", "0.11 mm", "0.09 mm", &
			"0.07 mm", "0.06 mm"
		do i = 1, M
			do k = 1, Nn(i)
				if (eta_strat(i,k)==eta_strat(i,k)) then
					write(f,'(36(F10.5,:,","))') x(i), eta_strat(i,k), B(i), &
						eta(i), H(i), tausg(i), qbT(i)*B(i)*2650.*1000., &
						cdf_strat(i,k,:)
				end if
			end do
			write(f,'(36(F10.5,:,","))') x(i), eta(i), B(i), &
				eta(i), H(i), tausg(i), qbT(i)*B(i)*2650.*1000., cdf_surf(i,:)
		end do
		close(unit=f)
	end if
	! Calculate bed evolution -------------------------------------------------!
	call bedevolv(qbT, qbTf, eta, M, N, Nn, npp, np, di, psi, dx, time, Ls, &
		dt, x, pdf_feed, pdf_strat, pdf_load, pdf_surf, d90surf, d50surf, &
		dgsurf, sigsurf, au, nactive, atrans, lps, Laold, eta_strat, o)
	! Update time -------------------------------------------------------------!
	! time = time + dt
	f = f + 1
	tt = tt + 1
end do

 call cpu_time(finish)
 print '("Total time elapsed= ",f7.3," minutes")',(finish-start)/60.0

stop
end program model_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Subroutines
!------------------------------------------------------------------------------!
subroutine renormalize(xp, np)
	implicit none
	! Renormalizes PDF gran size distribution, eliminating negative values
	integer :: j
	integer, intent(in) :: np
	real :: s
	real, intent(inout), dimension(1:np) :: xp
	do j = 1, np
		if (xp(j) < 0) then
			xp(j) = 0.0
		end if
	end do
	s = sum(xp)
	if (s > 0) then
		do j = 1, np
			xp(j) = xp(j)/s
		end do
	else
		do j = 1, np
			xp(j) = 0.0
		end do
	end if
end subroutine
!------------------------------------------------------------------------------!
subroutine cdf2pdf(npp, xpf, xp)
	implicit none
	! Converts sediment grain size distribution from a CDF in percent values to
	! a PDF in fraction values
	integer :: j
	integer, intent(in) :: npp
	real, intent(in), dimension(1:npp) :: xpf
	real, intent(out), dimension(1:npp-1) :: xp
	xp = (/ ((xpf(j-1)-xpf(j))/100., j=2, npp) /)
end subroutine
!------------------------------------------------------------------------------!
subroutine pdf2cdf(npp, xp, xpf)
	implicit none
	! Converts sediment grain size distribution froma  PDF in fraction values to
	! a CDF in percent values
	integer :: j
	integer, intent(in) :: npp
	real, intent(in), dimension(1:npp-1) :: xp
	real, intent(out), dimension(1:npp) :: xpf
	xpf(1) = 100.0
	do j = 2, npp
		xpf(j) = xpf(j-1)-100.0*xp(j-1)
	end do
	xpf(npp) = 0.0
end subroutine
!------------------------------------------------------------------------------!
subroutine geomean(np, psi, xp, xdm, xsg)
	implicit none
	! Calculates the geometric mean of a grain size distribution
	integer :: j
	real :: psibar, s2
	integer, intent(in) :: np
	real, intent(in), dimension(1:np) :: psi, xp
	real, intent(out) :: xdm, xsg
	psibar = sum((/ (psi(j)*xp(j), j=1, np) /))
	xdm = 2.**psibar
	s2 = sum((/ ((psi(j)-psibar)**2.*xp(j), j=1, np) /))
	xsg = 2.**s2
end subroutine
!------------------------------------------------------------------------------!
subroutine finer(npp, di, xpf, xty, xtdy)
	implicit none
	! Determines grain diameter for which xty percent of the grain size 
	! distribution is finer by weight
	integer :: j, z
	integer, intent(in) :: npp
	real, intent(in), dimension(1:npp) :: di, xpf
	real, intent(in) :: xty
	real, intent(out) :: xtdy
	do j = 1, npp-1
		if (xpf(j) >= xty) then
			if (xpf(j+1) < xty) then
				xtdy = exp(log(di(j+1))+(log(di(j))-log(di(j+1)))/ &
					(xpf(j)-xpf(j+1))*(xty-xpf(j+1)))
			end if
		end if
	end do
end subroutine
!------------------------------------------------------------------------------!
subroutine backwater(qw, eta, dx, M, g, Rr, ar, nk, d90surf, dgsurf, H, tausg, &
	time)
	implicit none
	! Calculates water depths and associated bed Shields stresses using the
	! standard step backwater procedure
	integer :: i, z
	integer, intent(in) :: M
	real, dimension(1:M) :: ks, dfhf, Hn, Hc, S0, U, hv, Fr2, E, Cf, Sf, Sfb, hf
	real, dimension(1:M) :: B
	real, intent(in) :: dx, g, ar, nk, Rr
	real, intent(in), dimension(1:M) :: qw, eta, d90surf, dgsurf
	real, intent(out), dimension(1:M) :: H, tausg
	real, intent(inout) :: time
	real :: H_sub, H_sup, SpFc_sub, SpFc_sup
	B = 0.035*(0.3048**3.)/qw
	S0(1) = (eta(1)-eta(2))/dx
	S0(2:M-1) = (/ ((eta(i-1)-eta(i+1))/2./dx, i=2, M-1) /)
	S0(M) = (eta(M-1)-eta(M))/dx
	ks = nk*d90surf/1000.
	Hc = ((qw**2.)/g)**(1./3.)
	Hn = ((ks**(1./3.)*qw**2.)/(ar**2.*g*S0))**(3./10.)
	dfhf(:) = 1.
	do i = M, 1, -1
		z = 0
		if (i == M) then
			H(i) = Hn(i) ! Downstream boundary condition, normal depth (m)
		else
			do while (abs(dfhf(i))>10.**(-5.))
				if (z == 0) then
					H(i) = H(i+1) ! Initial water depth guess
				else
					H(i) = H(i)-sign(0.00001,dfhf(i)) ! Adjust water depth guess
				end if
				U(i) = qw(i)/H(i) ! Mean flow velocity (m/s)
				hv(i) = U(i)**2./2./g ! Velocity head (m)
				Fr2(i) = U(i)**2./g/H(i) ! Froude number
				E(i) = eta(i)+H(i)+hv(i) ! Total energy (m)
				Cf(i) = (ar*(H(i)/ks(i))**(1./6.))**(-2.) ! Flow resistance
				Sf(i) = Cf(i)*Fr2(i) ! Friction slope
				Sfb(i) = 0.5*(Sf(i)+Sf(i+1))
				hf(i) = Sfb(i)*dx
				dfhf(i) = E(i)-(E(i+1)+hf(i))
				z = z+1
				if (z > 1000) then
					H(i) = 0.95*Hc(i) ! Critical depth after 1000 iterations
					dfhf(i) = 0.0
				end if
			end do
		end if
		U(i) = qw(i)/H(i) ! Mean flow velocity (m/s)
		hv(i) = U(i)**2./2./g ! Velocity head (m)
		Fr2(i) = U(i)**2./g/H(i) ! Froude number
		E(i) = eta(i)+H(i)+hv(i) ! Total energy (m)
		Cf(i) = (ar*(H(i)/ks(i))**(1./6.))**(-2.) ! Flow resistance
		Sf(i) = Cf(i)*Fr2(i) ! Friction slope
		tausg(i) = Cf(i)*U(i)**2./Rr/g/(dgsurf(i)/1000.) ! Shields stress
	end do
	do i = 1, M ! Determine supercritical flow depths
		if (Fr2(i) >= 0.8) then ! .or. H(i) <= 1.06*Hc(i)) then
			if (i == 1) then
				H(i) = Hn(i)
			else
				SpFc_sub = (qw(i-1)*B(i-1))**2.0/(B(i-1)*H(i-1))/g+ &
					0.5*H(i-1)**2.0*B(i-1)
				H_sup = 10.0**(-10.0)
				H_sub = 10.0
				z = 0
				dfhf(i) = 1.0
				do while (z < 1)
					H_sup = H_sup+10.0**(-7.0)
					SpFc_sup = (qw(i)*B(i))**2.0/(B(i)*H_sup)/g+ &
						0.5*H_sup**2.0*B(i)
					dfhf(i) = abs(SpFc_sub-SpFc_sup)
					if (dfhf(i) > H_sub) then
						z = 1
					end if
					if (H_sup >= Hc(i)) then
						write(*,*) 'ERROR: No solution found for flow depth!'
						H_sup = 1.05*Hc(i)
						z = 1
					end if
					H_sub = dfhf(i)
				end do
				H(i) = H_sup
			end if
		end if
		U(i) = qw(i)/H(i) ! Mean flow velocity (m/s)
		hv(i) = U(i)**2./2./g ! Velocity head (m)
		Fr2(i) = U(i)**2./g/H(i) ! Froude number
		E(i) = eta(i)+H(i)+hv(i) ! Total energy (m)
		Cf(i) = (ar*(H(i)/ks(i))**(1./6.))**(-2.) ! Flow resistance
		Sf(i) = Cf(i)*Fr2(i) ! Friction slope
		tausg(i) = Cf(i)*U(i)**2./Rr/g/(dgsurf(i)/1000.) ! Shields stress
		if (H(i) /= H(i)) then
			write(*,*) 'ERROR: NaN value for flow depth!'
			time = 10.0**6.0
		end if
	end do
end subroutine
!------------------------------------------------------------------------------!
subroutine bedload(tausg,dgsurf,pdf_surf,ds,M,np,Rr,g,qbT,pdf_load,time)
	implicit none
	! Calculates the volumetric bedload per unit width and corresponding grain
	! grain sizes being transported using a modified version of the Ashida and
	! Michiue bedload relation
	integer :: i, j, M, np
	real, parameter :: &
		taussrg = 0.043, &! Critical Shields stress for surface geometric mean
		alpha = 0.27 !
	real :: npar, gNAM, tauci, taubi, g, Rr
	real, intent(inout) :: time
	real, dimension(1:np) :: qbi
	real, dimension(1:np), intent(in) :: ds
	real, dimension(M,np), intent(in) :: pdf_surf
	real, dimension(1:M), intent(in) :: tausg, dgsurf
	real, dimension(1:M), intent(out) :: qbT
	real, dimension(M,np), intent(out) :: pdf_load
	do i = 1, M
		do j = 1, np
			npar = ds(j)/dgsurf(i)
			if (npar > 1) then
				gNAM = (1./npar)**0.68
			else
				gNAM = (1./npar)**0.98
			end if
			tauci = taussrg*gNAM
			taubi = tausg(i)*dgsurf(i)/ds(j)
			if (taubi <= tauci) then
				qbi(j) = 0.
			else
				qbi(j) = 17.*alpha*(taubi-tauci)*(sqrt(taubi)-sqrt(tauci))* &
					sqrt(Rr*g*ds(j)/1000.)*ds(j)/1000.*pdf_surf(i,j)
			endif
		end do
		qbT(i) = sum(qbi)
		if (qbT(i) > 0) then
			do j = 1, np
				pdf_load(i,j) = qbi(j)/qbT(i)
			end do
		else
			do j = 1, np
				pdf_load(i,j) = 0.0
			end do
		end if
		call renormalize(pdf_load(i,:),np)
		do j = 1, np
			if (pdf_load(i,j) /= pdf_load(i,j)) then
				write(*,'("ERROR: NaN value for bedload grain size!")')
				time = 10.0**6.0
			end if
		end do
	end do
end subroutine
!------------------------------------------------------------------------------!
subroutine bedevolv(qbT, qbTf, eta, M, N, Nn, npp, np, di, psi, dx, time, Ls, &
	dt, x, pdf_feed, pdf_strat, pdf_load, pdf_surf, d90surf, d50surf, dgsurf, &
	sigsurf, au, nactive, atrans, lps, Laold, eta_strat, o)
	implicit none
	integer :: i, M, N, np, k, npp, j, Nnew, o, ii
	integer, dimension(1:M), intent(inout) :: Nn
	real :: qbTdev1, qbTdev2, qbTdev, dgsurft, sigsurft, d90surft, d50surft, &
		Ls, des, qe, qw, A1, A2, Xc, Yc, etaf1, etaf2, time
	real, dimension(1:M) :: La, etanew, de
	real, dimension(1:np) :: qjj1dev, qjj2dev, pdf_feed, Ft, psi
	real, dimension(1:npp) :: Fft
	real, dimension(1:npp), intent(in) :: di
	real, intent(in) :: dx, qbTf, au, nactive, atrans, lps, dt
	real, dimension(1:M), intent(inout) :: Laold, eta, d90surf, d50surf, &
		dgsurf, sigsurf
	real, dimension(1:M), intent(in) :: qbT, x
	real, dimension(M,np) :: Fnew, FIexc
	real, dimension(M,np), intent(in) :: pdf_load
	real, dimension(M,np), intent(inout) :: pdf_surf
	real, dimension(M,N), intent(inout) :: eta_strat
	real, dimension(M,N,np), intent(inout) :: pdf_strat
	do i = 1, M
		La(i) = nactive*d90surf(i)/1000. ! Active layer thickness (m)
		if (i == 1) then ! First node ------------------------------------------
			qbTdev = au*(qbT(i)-qbTf)/dx+(1.-au)*(qbT(i+1)-qbT(i))/dx
			!qe = qbT(i)+(1./8.)*(qbT(i)-qbTf)+(3./8.)*(qbT(i+1)-qbT(i))
			!qw = qbTf+(3./8.)*(qbT(i)-qbTf)
			!qbTdev = (qe-qw)/dx
			if (qbTdev>0) then ! For degradation
				FIexc(i,:) = pdf_strat(i,Nn(i),:)
			else ! For aggradation
				FIexc(i,:) = atrans*pdf_surf(i,:)+(1.-atrans)*pdf_load(i,:)
			end if
			qjj1dev = au*(qbT(i)*pdf_load(i,:)-qbTf*pdf_feed)/dx+ &
				(1.-au)*(qbT(i+1)*pdf_load(i+1,:)-qbT(i)*pdf_load(i,:))/dx
			qjj2dev = FIexc(i,:)*qbTdev
		!else if (i == 2) then ! Second node------------------------------------
		!	qe = qbT(i)+(1./8.)*(qbT(i)-qbT(i-1))+(3./8.)*(qbT(i+1)-qbT(i))
		!	qw = qbT(i-1)+(1./8.)*(qbT(i-1)-qbTf)+(3./8.)*(qbT(i)-qbT(i-1))
		!	qbTdev = (qe-qw)/dx
		!	if (qbTdev>0) then ! For degradation
		!		FIexc(i,:) = pdf_strat(i,Nn(i),:)
		!	else ! For aggradation
		!		FIexc(i,:) = atrans*pdf_surf(i,:)+(1.0-atrans)*pdf_load(i,:)
		!	end if
		!	qjj1dev = au*(qbT(i)*pdf_load(i,:)-qbT(i-1)*pdf_load(i-1,:))/ &
		!		dx + (1.-au)*(qbT(i+1)*pdf_load(i+1,:)-qbT(i)* &
		!		pdf_load(i,:))/dx
		!	qjj2dev = FIexc(i,:)*qbTdev
		else
			if (i == M) then ! Last node----------------------------------------
				qbTdev = (qbT(M)-qbT(M-1))/dx
				if (qbTdev>0) then ! For degradation
					FIexc(i,:) = pdf_strat(i,Nn(i),:)
				else ! For aggradation
					FIexc(i,:) = atrans*pdf_surf(i,:)+(1.-atrans)*pdf_load(i,:)
				end if
				qjj1dev = (qbT(i)*pdf_load(i,:)-qbT(i-1)*pdf_load(i-1,:))/dx
				qjj2dev = FIexc(i,:)*qbTdev
			else ! In between nodes---------------------------------------------
				qbTdev1 = au*(qbT(i)-qbT(i-1))/dx
				qbTdev2 = (1.-au)*(qbT(i+1)-qbT(i))/dx
				qbTdev = qbTdev1+qbTdev2
				!qe = qbT(i)+(1./8.)*(qbT(i)-qbT(i-1))+(3./8.)*(qbT(i+1)-qbT(i))
				!qw = qbT(i-1)+(1./8.)*(qbT(i-1)-qbT(i-2))+(3./8.)*&
				!	(qbT(i)-qbT(i-1))
				!qbTdev = (qe-qw)/dx
				if (qbTdev>0) then ! For degradation
					FIexc(i,:) = pdf_strat(i,Nn(i),:)
				else ! For aggradation
					FIexc(i,:) = atrans*pdf_surf(i,:)+(1.0-atrans)*pdf_load(i,:)
				end if
				qjj1dev = au*(qbT(i)*pdf_load(i,:)-qbT(i-1)*pdf_load(i-1,:))/ &
					dx + (1.-au)*(qbT(i+1)*pdf_load(i+1,:)-qbT(i)* &
					pdf_load(i,:))/dx
				qjj2dev = FIexc(i,:)*qbTdev
			end if !------------------------------------------------------------
		end if
		etanew(i) = eta(i)+dt*(-qbTdev/(1.0-lps))
		etanew(M) = eta(M)
		do j = 1, np
			Fnew(i,j) = pdf_surf(i,j)+dt*(-qjj1dev(j)+qjj2dev(j))/(1.-lps)/La(i)
			if (time > 0 .and. La(i) /= Laold(i)) then
				Fnew(i,j) = Fnew(i,j)+(FIexc(i,j)-pdf_surf(i,j))/&
					(La(i)/(La(i)-Laold(i)))
			end if
		end do
	end do
	i = 1
	do while (i <= M)
		if (etanew(i) /= etanew(i)) then
			write(*,'("ERROR: NaN value for bed elevation!")')
			i = M
			time = 10.0**6.0
		end if
		j = 1
		do while (j <= np)
			Ft(j) = Fnew(i,j)
			if (Ft(j) /= Ft(j)) then
				write(*,'("ERROR: NaN value for surface grain size!")')
				j = np
				i = M
				time = 10.0**6.0
			!else if (Ft(j) < 0) then
			!	Ft(j) = 0.0
			end if
			j = j + 1
		end do
		Laold(i) = La(i)
		call renormalize(Ft, np)
		call geomean(np, psi, Ft, dgsurft, sigsurft)
		dgsurf(i) = dgsurft
		sigsurf(i) = sigsurft
		call pdf2cdf(npp, Ft, Fft)
		do j = 1, npp
			if (Fft(j) < 0) then
				Fft(j) = 0.0
			end if
		end do
		call finer(npp, di, Fft, 90., d90surft)
		call finer(npp, di, Fft, 50., d50surft)
		d90surf(i) = d90surft
		d50surf(i) = d50surft
		call cdf2pdf(npp, Fft, pdf_surf(i,:))
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! STRATIGRAPHY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Nnew = nint((etanew(i)-(nactive*d90surf(i)/1000.))/Ls)+ 2
		de(i) = etanew(i)-eta(i)
		des = (etanew(i)-(nactive*d90surf(i)/1000.))-(eta(i)-La(i))
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (de(i) > 0) then! AGGRADATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!if (des < 0) then
			!	des = 0.0
			!end if
			if (Nnew > Nn(i)) then ! More stratigraphy nodes
				pdf_strat(i,Nn(i),:) = (pdf_strat(i,Nn(i)- 1,:)* &
					(eta_strat(i,Nn(i))-(Nn(i)- 2)*Ls))/Ls+ &
					(FIexc(i,:)*((Nn(i)- 1 )*Ls-eta_strat(i,Nn(i))))/Ls
				do k = 1, Nnew
					if (k > Nn(i)) then
						pdf_strat(i,k,:) = FIexc(i,:)
					end if
				end do
			else if (Nnew == Nn(i)) then ! Same number of stratigraphy nodes
				if (des > 0) then
					pdf_strat(i,Nnew,:) = (pdf_strat(i,Nn(i),:)* &
						(eta_strat(i,Nn(i))-(Nn(i)- 1)*Ls)+des*FIexc(i,:))/ &
						((eta_strat(i,Nn(i))-(Nn(i)- 1)*Ls)+des)
				end if
			else if (Nnew < Nn(i)) then ! Less stratigraphy nodes
				do k = 1, Nn(i)
					if (k > Nnew) then
						do j = 1, np
							pdf_strat(i,k,j) = 0./o
						end do
					end if
				end do
			end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else if (de(i) < 0) then! DEGRADATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (Nnew < Nn(i)) then ! Less stratigraphy nodes
				do k = 1, Nn(i)
					if (k > Nnew) then
						do j = 1, np
							pdf_strat(i,k,j) = 0./o
						end do
					end if
				end do
			else if (Nnew >= Nn(i)) then ! More stratigraphy nodes
				!if (des > Ls) then
					do k = 1, Nnew
						if (k > Nn(i)) then
							pdf_strat(i,k,:) = pdf_surf(i,:)
						end if
					end do
				!end if
			end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		end if
		Nn(i) = Nnew
		k = 1
		do while (k <= N)
			if (k == 1) then
				eta_strat(i,k) = 0.
			else if (k < Nn(i)) then
				eta_strat(i,k) = eta_strat(i,k-1)+Ls
			else if (k == Nn(i)) then
				eta_strat(i,k) = etanew(i)-(nactive*d90surf(i)/1000.)
			else if (k > Nn(i)) then
				eta_strat(i,k) = 0./o
				do j = 1, np
					pdf_strat(i,k,j) = 0./o
				end do
			end if
			if (eta_strat(i,k) > (etanew(i)-(nactive*d90surf(i)/1000.))) then
				eta_strat(i,k) = 0./o
			end if
			call renormalize(pdf_strat(i,k,:), np)
			if (k <= Nn(i)) then
				j = 1
				do while (j <= np)
					if (pdf_strat(i,k,j) /= pdf_strat(i,k,j)) then
						write(*,'("ERROR: NaN for stratigraphy grain size")')
						time = 10.0**6.0
						j = np
						i = M
						k = N
					end if
					j = j + 1
				end do
			end if
			k = k + 1
		end do
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		eta(i) = etanew(i)
		i = i + 1
	end do
end subroutine

