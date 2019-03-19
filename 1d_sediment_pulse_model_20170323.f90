program model_1d
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------!
!                                                                      !
!              ONE-DIMENSIONAL HYDRO-MORPHODYNAMIC MODEL               !
!                                                                      !
!----------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Written by: Jacob A. Morgan
!			  Colorado State University
!			  Fort Collins, Colorado
!
! Creation: 31 July 2015
! Last Edit: 23 March 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 character (LEN=15) :: fname
 character (LEN=5) :: TIMEM
 character (LEN=3) :: ii
 character (LEN=2) :: RUN
 character (LEN=1) :: junk
logical :: dir_e
integer :: &
	i, &! Downstream spatial index
	j, &! Grain size distiribution size class index
	k, &! Vertical grain sorting index
	f, &! Index for output file units
	o, &
	tt, &
	NR, ios, v1
real, allocatable :: dat(:,:)
integer, parameter :: &
	! M = 229, &! Number of downstream nodes
	nw = 10, &! Number of wavelenths of width variation
	npp = 29, &! Number of grain size classes in CDF distributions
	np = npp-1, &! Number of grain size classes in PDF distributions
	N = 200, & ! Maximum number of vertical nodes for storing stratigraphy
	RUN_NO = 10 ! Run number
real, parameter, dimension(1:55) :: &
	MP = (/ 4.05, 4.05, 8.10, & ! Runs 01-03
 		4.05, 8.10, 16.20, 32.40, 64.80, & ! Runs 04-08
 		4.05, 8.10, 16.20, 32.40, 64.80, & ! Runs 09-13
 		4.05, 8.10, 16.20, 32.40, 64.80, & ! Runs 14-18
 		4.05, 8.10, 16.20, 32.40, 64.80, & ! Runs 19-23
 		4.05, 8.10, 16.20, 32.40, 64.80, & ! Runs 24-28
 		32.40, 64.80, 64.80, & ! Runs 29-31
 		16.20, 16.20, 16.20, 16.20, & ! Runs 32-35
 		16.20, 16.20, 16.20, 16.20, & ! Runs 36-39
 		16.20, 16.20, 16.20, 16.20, & ! Runs 40-43
 		16.20, 16.20, 16.20, 16.20, & ! Runs 44-47
 		16.20, 16.20, 16.20, 16.20, & ! Runs 48-51
 		16.20, 16.20, 16.20, 16.20 /), & ! Runs 52-55
	TP = (/ 108.0, 54.0, 108.0, & ! Runs 01-03
		27.0, 54.0, 108.0, 216.0, 432.0, & ! Runs 04-08
		13.5, 27.0, 54.0, 108.0, 216.0, & ! Runs 09-13
		6.75, 13.5, 27.0, 54.0, 108.0, & ! Runs 14-18
		3.375, 6.75, 13.5, 27.0, 54.0, & ! Runs 19-23
		1.6875, 3.375, 6.75, 13.5, 27.0, & ! Runs 24-28
		6.75, 13.5, 6.75, & ! Runs 29-31
		27.0, 27.0, 27.0, 27.0, & ! Runs 32-35
		27.0, 27.0, 27.0, 27.0, & ! Runs 36-39
		27.0, 27.0, 27.0, 27.0, & ! Runs 40-43
		27.0, 27.0, 27.0, 27.0, & ! Runs 44-47
		27.0, 27.0, 27.0, 27.0, & ! Runs 48-51
		27.0, 27.0, 27.0, 27.0 /) ! Runs 52-55
real, parameter :: &
	! Physical Constants
	pi = 3.14159265359, &! Constant, pi
	g = 9.81, &! Gravitational acceleration (m/s^2)
	rho = 1000., &! Density of water (kg/m^3)
	Rr = 1.65, &! Submerged specific gravity of sediment
	lps = 0.4, &! Sediment porosity
	! Geometric Constants
	B0 = 0.0925, &! Mean half-width of width variations (m)
	Ac = 0.10, &!0.1676, &! Dimensionless amplitude of width variations
	lambdac = 0.2, &! Dimensionless wavenumber of width variations
	Lc = 2.0*pi*B0/lambdac, &!1.08, &! Wavelength of width variations (m)
	enL = 1.00, &!1.58, &! Length of flume entrance reach (m)
	exL = 1.00, &!1.08, &! Length of flume exit reach (m)
	enB = B0, &! Half-width of flume entrance reach (m)
	exB = B0, &! Half-width of flume exit reach (m)
	Sb = 0.005, &! Initial constant bed slope
	etad = 0.048, &! Fixed elevation of downstream-most node (m)
	L = 40, &!enL+nw*Lc+exL, &! Total length of flume (m)
	Ls = 0.0025, &! Thickness of substrate layers (m)
	! Numerical Constants
	nk = 2.0, &! Factor multiplied by surface D90 for roughness height, ks
	nactive = 2.0, &! Factor multiplied by surface D90 for active layer
	ar = 8.1, &! Coefficient in Manning-Strickler resistance relation
	au = 0.75, &! Upwinding coefficient for load derivatives in Exner equation
	atrans = 0.7, &! Coefficient for material transferred to substrate for agg.
	dx = 0.05, &!L/(M-1), &! Spatial step length (m)
	dt = 0.0625, &! Time step increment (s)
	T = 240.0*60.0*60.0, &!(13.+23.6+4.5+28.6+711./60.)*60.*60., &! Total run time (s)
	! Hydraulic and Morphodynamic Constants
	Q = 0.00091, &! Water discharge (m^3/s)
	!Qbf = 150./1000./2650./60., &! Sediment feed rate (m^3/s)
	!
	!
	!
	Qbf = MP(RUN_NO)/TP(RUN_NO)/2650.0/60.0
integer, parameter :: M = L/dx + 1
integer, dimension(1:M) :: &
	Nn ! Number of vertical nodes with stored stratigraphy data
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
	psi = log10(ds)/log10(2.) ! PDF psi scale grain size classes
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
write(RUN,"(I2.2)") RUN_NO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Set up Flume Geometry
!------------------------------------------------------------------------------!
do i = 1, M
	x(i) = (i-1)*dx
!	if (x(i) <= enL) then
!		Bx(i) = enB
!	else if (x(i) >= (enL+nw*Lc)) then
!		Bx(i) = exB
!	else
!		Bx(i) = B0*(1+Ac*sin(2*pi*(x(i)-enL)/Lc))
!	end if
	Bx(i) = B0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! START FROM PREVIOUS RUN !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read file for elevation and grain size values for all nodes from previous run
NR = 0
open(unit=1,&
    file="run00_12000.dat")
do i = 1, 1000000
    read(1,*,iostat=ios) junk
    if (ios /= 0) then
        exit
    else if (i == 1000000) then
        write(*,*) "Error: Maximum number of records exceeded..."
        write(*,*) "Exiting program now..."
        STOP
    end if
    NR = NR + 1
end do
rewind(1)
allocate(dat(NR,36))
do i = 1, NR
    read(1,*) dat(i,:)
end do
close(1)

! write(*,*) "HERE 1!"

! 
i = 1
k = 1
do v1 = 1, size(dat,1)
    if (dat(v1,2) /= dat(v1,4)) then
        !!!!!!! x(i) = dat(v1,1)
        eta_strat(i,k) = dat(v1,2)
        eta(i) = dat(v1,4)
        cdf_strat(i,k,:) = dat(v1,8:36)
        k = k + 1
    else
        !!!!!!! x(i) = dat(v1,1)
        eta(i) = dat(v1,4)
        H(i) = dat(v1,5)
        tausg(i) = dat(v1,6)
        qbT(i) = dat(v1,7)/B(i)/2650.0/1000.0
        cdf_surf(i,:) = dat(v1,8:36)
        i = i + 1
        k = 1
    end if
end do

!!!!!!!!do while (i <= M)
!!!!!!!!    Nn(i) =  nint((eta(i)-(nactive*d90surf(i)/1000.))/Ls)+ 2
!!!!!!!!    do while (k <= Nn(i))
!!!!!!!!        if (cdf_strat(i,k,1) == 0.0) then
!!!!!!!!            call pdf2cdf(npp, pdf_subi(:), cdf_strat(i,k,:))
!!!!!!!!        end if
!!!!!!!!    end do
!!!!!!!!    i = i + 1
!!!!!!!!end do
!
!!
i = 1
do while (i <= M)
    call finer(npp, di, cdf_surf(i,:), 90., d90surf(i))
    call finer(npp, di, cdf_surf(i,:), 50., d50surf(i))
    call cdf2pdf(npp, cdf_surf(i,:), pdf_surf(i,:))
    call geomean(np, psi, pdf_surf(i,:), dgsurf(i), sigsurf(i))
    Nn(i) =  nint((eta(i)-(nactive*d90surf(i)/1000.))/Ls)+ 2
    k = 1
    do while (k <= N)
        if (cdf_strat(i,k,1) == 0.0) then
            call pdf2cdf(npp, pdf_subi(:), cdf_strat(i,k,:))
        end if
        call finer(npp, di, cdf_strat(i,k,:), 90., d90strat(i,k))
        call finer(npp, di, cdf_strat(i,k,:), 50., d50strat(i,k))
        call cdf2pdf(npp, cdf_strat(i,k,:), pdf_strat(i,k,:))
        call geomean(np, psi, pdf_strat(i,k,:), dgstrat(i,k), sigstrat(i,k))
        if (k == 1) then
            eta_strat(i,k) = 0.0
        else if (k < Nn(i)) then
            eta_strat(i,k) = eta_strat(i,k-1)+Ls
        else if (k == Nn(i)) then
            eta_strat(i,k) = eta(i)-(nactive*d90surf(i)/1000.0)
        else if (k > Nn(i)) then
            eta_strat(i,k) = 0.0/o
            do j = 1, np
                pdf_strat(i,k,j) = 0.0/o
            end do
        end if
        if (eta_strat(i,k) > (eta(i)-(nactive*d90surf(i)/1000.))) then
            eta_strat(i,k) = 0.0/o
        end if
        call renormalize(pdf_strat(i,k,:), np)
        k = k + 1
    end do
    !!!!!!!!!! pdf_strat(i,Nn(i),:) = pdf_subi !!!!!!!!!!
    i = i + 1
end do

! Read D90 surface grain sizes from previous time step (for Laold)
NR = 0
open(unit=10001,&
	file="Z:/JAM/model/RUNS/run00_11990.dat",&
	form="formatted",&
	status="old",&
	action="read")

do i = 1, 1000000
    read(10001,*,iostat=ios) junk
    if (ios /= 0) then
        exit
    else if (i == 1000000) then
        write(*,*) "Error: Maximum number of records exceeded..."
        write(*,*) "Exiting program now..."
        STOP
    end if
    NR = NR + 1
end do
rewind(10001)

i = 1
do while (i <= M)
    read(10001,*) Laold(i)
    Laold(i) = nactive*Laold(i)/1000.0
    i = i + 1
end do
 close(unit=10001)

! write(*,*) "HERE 2!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tt = 720000/dt
time = tt*dt
write(TIMEM,"(I4.4)") int(time/60.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		fname = "run_"//TIMEM//"m.csv"
!		open(unit=1,file="OUTPUT_00/"//fname,form="formatted",&
!			status="replace")
!		write(1,'("#","time = ",A4," minutes")') TIMEM
!		write(1,'("#",36(A10,:,","))') "x(m)", "z(m)", "B(m)", "eta(m)", &
!			"H(m)", "tausg(-)", "QbT(g/s)", "8.00 mm", "6.73 mm", "5.66 mm", &
!			"4.76 mm", "4.00 mm", "3.36 mm", "2.83 mm", "2.38 mm", "2.00 mm", &
!			"1.68 mm", "1.41 mm", "1.19 mm", "1.00 mm", "0.84 mm", "0.71 mm", &
!			"0.59 mm", "0.50 mm", "0.42 mm", "0.35 mm", "0.30 mm", "0.25 mm", &
!			"0.21 mm", "0.18 mm", "0.15 mm", "0.13 mm", "0.11 mm", "0.09 mm", &
!			"0.07 mm", "0.06 mm"
!		do i = 1, M
!			do k = 1, Nn(i)
!				if (eta_strat(i,k)==eta_strat(i,k)) then
!					write(1,'(36(F10.5,:,","))') x(i), eta_strat(i,k), B(i), &
!						eta(i), H(i), tausg(i), qbT(i)*B(i)*2650.0*1000.0, &
!						cdf_strat(i,k,:)
!				end if
!			end do
!			write(1,'(36(F10.5,:,","))') x(i), eta(i), B(i), &
!				eta(i), H(i), tausg(i), qbT(i)*B(i)*2650.0*1000.0, cdf_surf(i,:)
!		end do
!		close(unit=1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! qbTf = Qbf/B(1)
 call cdf2pdf(npp, (/ 100., 100., 100., 100., 100., 100., 100., &
	100., 100., 100., 100., 100., 100., 0., 0., 0., 0., 0., 0., &
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /), pdf_feed)
 call bedevolv(qbT, qbTf, eta, M, N, Nn, npp, np, di, psi, dx, time, Ls, &
	dt, x, pdf_feed, pdf_strat, pdf_load, pdf_surf, d90surf, d50surf, &
	dgsurf, sigsurf, au, nactive, atrans, lps, Laold, eta_strat, o)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
! Main model
!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Laold = 0.0
! time = 0.0
f = 2
tt = tt + 1
 call cpu_time(start)
do while (time <= T)
	time = tt*dt
	if (time < 720000.0 + 60.0*TP(RUN_NO)) then
		qbTf = Qbf/B(1)
	else
		qbTf = 0.0
	end if
	! Print current time step to disply ---------------------------------------!
	if (time == 0.0 .or. (time)/int(time) == 1.0) then
		!print*, "time = ", nint(time), "s"
		write(*,'("    time = ",I7.7," seconds")') int(time)
		if (time > 0.0) then
			write(*,'("max Cr # = ",F5.3,"")') (maxval(qw/H)*(dt/dx))
		end if
	end if
	! Calculate backwater and bedload -----------------------------------------!
	call backwater(qw, eta, dx, M, g, Rr, ar, nk, d90surf, dgsurf, H, tausg, &
		time)
	call AandMmod(tausg,dgsurf,pdf_surf,ds,M,np,Rr,g,qbT,pdf_load,time)
	! Write output to file ----------------------------------------------------!
	write(TIMEM,"(I5.5)") int(time/60.0)
	if (time == 0 .or. (time/60.0)/int(time/60.0) == 1.0) then
		do i = 1, M
			do k = 1, Nn(i)
				call pdf2cdf(npp, pdf_strat(i,k,:), cdf_strat(i,k,:))
			end do
			call pdf2cdf(npp, pdf_surf(i,:), cdf_surf(i,:))
		end do
		! Check if output directy exists and open output file
		!INQUIRE(DIRECT="OUTPUT_"//RUN, EXIST=dir_e)
		fname = "run"//RUN//"_"//TIMEM//".csv"
		call system("mkdir OUTPUT_"//RUN//"\")
		open(unit=f,file="OUTPUT_"//RUN//"\"//fname,form="formatted",&
			status="replace")
		!if (dir_e) then
		!	open(unit=f,file="OUTPUT_"//RUN//"/"//fname,form="formatted",&
		!		status="replace")
		!else
		!	call system("mkdir OUTPUT_"//RUN//"/")
		!	open(unit=f,file="OUTPUT_"//RUN//"/"//fname,form="formatted",&
		!		status="replace")
		!end if
		! Write data to file
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
						eta(i), H(i), tausg(i), qbT(i)*B(i)*2650.0*1000.0, &
						cdf_strat(i,k,:)
				end if
			end do
			write(f,'(36(F10.5,:,","))') x(i), eta(i), B(i), &
				eta(i), H(i), tausg(i), qbT(i)*B(i)*2650.0*1000.0, cdf_surf(i,:)
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
	B = 0.0321363*(0.3048**3.)/qw
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------!
!-                          BEDLOAD RELATIONS                         -!
!----------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AandMmod( & ! Modified Ashida and Michiue formula
	tausg, &
	dgsurf, &
	pdf_surf, &
	ds, &
	M, &
	np, &
	Rr, &
	g, &
	qbT, &
	pdf_load, &
	time)
	implicit none
	! Calculates the volumetric bedload transport rate per unit channel 
	! width and corresponding grain sizes in transport using the Ashida
	! and Michiue bedload (1972) relation as modified by Viparelli et
	! al. (2010).
	!
	! The Viparellie et al. (2010) modification of the Ashida and
	! Michiue (1972) relation is specified as:
	!  
	!                 q*bi = 17α(τ*i-τ*ci)*(√τ*i-√τ*ci)
	!
	! where, q*bi is the dimensionless bedload rate, or Einstein
	! parameter, of grain size class i, α is a calibration coefficient,
	! τ*i is the dimensionless shear stress, or Shields parameter,
	! referenced to grain size class i, and τ*ci is the critical Shields
	! parameter corresponding to the threshold of motion for grain size
	! class i.
	!
	! The critical Shields number for each grain size class (τ*ci), is
	! determined using the hiding function of Viparelli et al. (2010):
	!
	!                    | (Di/Dsg)^-0.98    for    (Di/Dsg) ≤ 1
	!       τ*ci/τ*scg = |
	!                    | (Di/Dsg)^-0.68    for    (Di/Dsg) > 1
	!
	! where τ*scg is the critical Shields stress referenced to the
	! geometric mean diameter of bed surface material, Di is the 
	! diameter of grain in grain size class i, and Dsg is the geometric
	! mean diameter of bed surface material.
	!
	! Viparelli et al. (2010) used linear regression of measured bedload
	! data to determine a value for coefficient value of α = 0.270 and
	! visually estimated τ*scg = 0.043. Using the flume experiments of
	! Nelson et al. (2015), we calibrated both α and τ*scg using
	! elevation profiles and bedload transport measured at the outlet of
	! the flume and obtained values of α = 1.0 and τ*scg = 0.03. These
	! parameters can be calibrated for specific situations (e.g., model
	! domain dimensions and grain sizes present).
	!
	! The dimensional bedload transport (qbi) is calculated from the
	! Einstein parameter as:
	!
	! qbi = q*bi×√(R×g×di)×di×
	!
	! where R is the submerged specific gravity of the sediment (assumed
	! to be 1.65), g is gravitational acceleration (assumed to be 9.81 
	! m/s^2).
	!
	! REFERENCES
	! Ashida, K. and M. Michiue (1972), Study on hydraulic resistance
	!     and bed-load transport rate in alluvial streams, Proceedings
	!     of the Japan Society of Civil Engineers, 206: 59-69.
	! Nelson, P.A., A.K. Brew, and J.A. Morgan (2015), Morphodynamic
	!     response of a variable-width channel to changes in sediment
	!     supply, Water Resources Research, 51(7): 5717-5734,
	!     doi: 10.1002/2014WR016806.
	! Viparelli, E., R. Haydel, M. Salvaro, P. R. Wilcock, and G. Parker
	!     (2010a), River morphodynamics with creation/consumption of
	!     grain size stratigraphy 1: laboratory experiments, Journal of
	!     Hydraulic Research, 48(6): 715-726,
	!     doi: 10.1080/00221686.2010.515383.
	!
	integer :: i, j, M, np
	real, parameter :: &
		taussrg = 0.029, &! Critical Shields stress for surface geometric mean
		alpha = 1.0 ! Coefficient
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
		dgsurf, sigsurf, qbT
	real, dimension(1:M), intent(in) :: x
	real, dimension(M,np) :: Fnew, FIexc
	real, dimension(M,np), intent(inout) :: pdf_surf, pdf_load
	real, dimension(M,N), intent(inout) :: eta_strat
	real, dimension(M,N,np), intent(inout) :: pdf_strat
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!	pdf_load(130,:) = pdf_load(130,:)*qbT(130)/(qbT(130)+0.048358605*qbTf) + &
!!!!		pdf_feed*0.048358605*qbTf/(qbT(130)+0.048358605*qbTf)
!!!!	call renormalize(pdf_load(130,:),np)
!!!!	qbT(130) = qbT(130) + 0.048358605*qbTf
!!!!	
!!!!	pdf_load(131,:) = pdf_load(131,:)*qbT(131)/(qbT(131)+0.074814208*qbTf) + &
!!!!		pdf_feed*0.074814208*qbTf/(qbT(131)+0.074814208*qbTf)
!!!!	call renormalize(pdf_load(131,:),np)
!!!!	qbT(131) = qbT(131) + 0.074814208*qbTf

!	pdf_load(130,:) = pdf_load(130,:)*qbT(130)/(qbT(130)+qbTf) + &
!		pdf_feed*qbTf/(qbT(130)+qbTf)
!	call renormalize(pdf_load(130,:),np)
!	qbT(130) = qbT(130) + qbTf

!	pdf_load(130,:) = pdf_load(130,:)*qbT(130)/(qbT(130)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(130)+0.1*qbTf)
!	call renormalize(pdf_load(130,:),np)
!	qbT(130) = qbT(130) + 0.1*qbTf
!	
!	pdf_load(131,:) = pdf_load(131,:)*qbT(131)/(qbT(131)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(131)+0.1*qbTf)
!	call renormalize(pdf_load(131,:),np)
!	qbT(131) = qbT(131) + 0.1*qbTf
!
!	pdf_load(132,:) = pdf_load(132,:)*qbT(132)/(qbT(132)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(132)+0.1*qbTf)
!	call renormalize(pdf_load(132,:),np)
!	qbT(133) = qbT(132) + 0.1*qbTf
!	
!	pdf_load(133,:) = pdf_load(133,:)*qbT(133)/(qbT(133)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(133)+0.1*qbTf)
!	call renormalize(pdf_load(133,:),np)
!	qbT(133) = qbT(133) + 0.1*qbTf
!	
!	pdf_load(134,:) = pdf_load(134,:)*qbT(134)/(qbT(134)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(134)+0.1*qbTf)
!	call renormalize(pdf_load(134,:),np)
!	qbT(134) = qbT(134) + 0.1*qbTf
!	
!	pdf_load(135,:) = pdf_load(135,:)*qbT(135)/(qbT(135)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(135)+0.1*qbTf)
!	call renormalize(pdf_load(135,:),np)
!	qbT(135) = qbT(135) + 0.1*qbTf
!	
!	pdf_load(136,:) = pdf_load(136,:)*qbT(136)/(qbT(136)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(136)+0.1*qbTf)
!	call renormalize(pdf_load(136,:),np)
!	qbT(136) = qbT(136) + 0.1*qbTf
!	
!	pdf_load(137,:) = pdf_load(137,:)*qbT(137)/(qbT(137)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(137)+0.1*qbTf)
!	call renormalize(pdf_load(137,:),np)
!	qbT(137) = qbT(137) + 0.103780754*qbTf
!	
!	pdf_load(138,:) = pdf_load(138,:)*qbT(138)/(qbT(138)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(138)+0.1*qbTf)
!	call renormalize(pdf_load(138,:),np)
!	qbT(138) = qbT(138) + 0.1*qbTf
!	
!	pdf_load(139,:) = pdf_load(139,:)*qbT(139)/(qbT(139)+0.1*qbTf) + &
!		pdf_feed*0.1*qbTf/(qbT(139)+0.1*qbTf)
!	call renormalize(pdf_load(139,:),np)
!	qbT(139) = qbT(139) + 0.1*qbTf
!!!!!!	pdf_load(115,:) = pdf_load(115,:)*qbT(115)/(qbT(115)+qbTf) + &
!!!!!!		pdf_feed*qbTf/(qbT(115)+qbTf)
!!!!!!	call renormalize(pdf_load(115,:),np)
!!!!!!	qbT(115) = qbT(115)+qbTf
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = 1, M
		La(i) = nactive*d90surf(i)/1000.0 ! Active layer thickness (m)
		if (i == 1) then ! First node ------------------------------------------
			qbTdev = au*(qbT(i)-qbTf)/dx+(1.0-au)*(qbT(i+1)-qbT(i))/dx
			if (qbTdev>0.0) then ! For degradation
				FIexc(i,:) = pdf_strat(i,Nn(i),:)
			else ! For aggradation
				FIexc(i,:) = atrans*pdf_surf(i,:)+(1.0-atrans)*pdf_load(i,:)
			end if
			qjj1dev = au*(qbT(i)*pdf_load(i,:)-qbTf*pdf_feed)/dx+ &
				(1.0-au)*(qbT(i+1)*pdf_load(i+1,:)-qbT(i)*pdf_load(i,:))/dx
			qjj2dev = FIexc(i,:)*qbTdev
		else
			if (i == M) then ! Last node----------------------------------------
				qbTdev = (qbT(M)-qbT(M-1))/dx
				if (qbTdev>0.0) then ! For degradation
					FIexc(i,:) = pdf_strat(i,Nn(i),:)
				else ! For aggradation
					FIexc(i,:) = atrans*pdf_surf(i,:)+(1.-atrans)*pdf_load(i,:)
				end if
				qjj1dev = (qbT(i)*pdf_load(i,:)-qbT(i-1)*pdf_load(i-1,:))/dx
				qjj2dev = FIexc(i,:)*qbTdev
			else ! In between nodes---------------------------------------------
				qbTdev = au*(qbT(i)-qbT(i-1))/dx+(1.0-au)*(qbT(i+1)-qbT(i))/dx
				!qbTdev = (qe-qw)/dx
				if (qbTdev>0.0) then ! For degradation
					FIexc(i,:) = pdf_strat(i,Nn(i),:)
				else ! For aggradation
					FIexc(i,:) = atrans*pdf_surf(i,:)+(1.0-atrans)*pdf_load(i,:)
				end if
				qjj1dev = au*(qbT(i)*pdf_load(i,:)-qbT(i-1)*pdf_load(i-1,:))/ &
					dx + (1.0-au)*(qbT(i+1)*pdf_load(i+1,:)-qbT(i)* &
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!		if (i < 116) then
!!!!!!			etanew(i) = eta(i)
!!!!!!			call pdf2cdf(npp, pdf_surf(i,:), Fft)
!!!!!!		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

