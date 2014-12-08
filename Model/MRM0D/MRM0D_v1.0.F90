PROGRAM MAIN
IMPLICIT NONE
!		"model of particle remineralization that includes bacterial 
!       extracellular enzyme dynamics" 
!
!----------------------------------
!  DESCRIPTION OF VARIABLES
!----------------------------------
!
!   	FLUX IN         
!
!				jday:			julian day
!				
!				hhmm:			hour, minute 
!				
!				pom_in_top:		particle flux at upper boundary (mg C m^-2)
!
!
!   	WATER COLUMN PROFILE
!
!				depth:			depth in the water column
!
!				temp:			temperature of seawater (°C)
!
!				s:				salinity of seawater (psu)
!
!				rho:			density of seawater (kg m^-3)
!  
!				sldoc:			semi-labile dissolved organic carbon (mg C m^-3)
!
!               
!       INITIAL STATE VARIABLES
!               
!               hdom:			hydrolysate dissolved organic matter (mg C m^-3)
!           
!               pom:			particulate organic matter (mg C m^-3)
!                       
!               bf:				free living bacteria (mg C m^-3)
!
!               ba:				attached bacteria (mg C m^-3)
!
!               enz:			enzyme in particle (mg C m^-3)
!	
!				h:				hydrolysate in particle (mg C m^-3)
!
!				denz:			inactive enzyme in particle (mg C m^-3)
!
!				edom:			active enzyme in dissolved environment (mg C m^-3)
!
!				dedom:			inactive enzyme in dissolved environment (mg C m^-3)
!
!
!       TEMPERATURE DEPENDENCE
!
!               q10_b:			q10 value for bacterial temperature dependence.
!							
!				q10_enz:		q10 value for enzyme activity rate.
!		
!				q10_hl:			q10 value for enzyme half life.  
!
!
!		WATER COLUMN PROPERTIES
!				
!				dz:				depth interval (m)
!				
!				dt:				time step (h)
!
!
!		PARTICLE PROPERTIES
!
!				w:				sinking speed of the particles (m h^-1)
!
!				r_top:			radius of the particle at the upper boundary (m)
!
!				ba_poc:			ratio of ba to pom at upper boundary (Ghiglione et al. 2009)
!
!
!       BACTERIA PARAMETERS
!
!               up_max_bf:		maximum free-living bacteria uptake rate at 20°C (h^-1)
!
!				up_max_ba:		maximum attached bacteria uptake rate at 20°C (h^-1)
!							
!               kh:				half saturation constant for attached bacterial uptake of hydrolysate (mg C m^-3 particle)
!
!               kdom:			half saturation constant for free-living bacterial uptake of dissolved hydrolysate (mg C m^-3 seawater)
!
!				bge_max:		maximum bacterial growth efficiency (Fig. 3, "Estuary", Del Giorgio and Cole 2000)
!
!				basal:			basal energetic cost at 20°C (h^-1) (based on the 0.01 day^-1 at 0°C from COBALT v1, Stock et al. in prep)
!		
!				quorum:			bacterial biomass at quorum sensing threshold (mg C m^-3 particle) 
!								**converted from 5E7 CFU (colony-forming units = active bacteria cells) mL^-1 (Gram et al. 2002)
!				
!				x:				exponent for slope of the quorum threshold response
!
!				max_epsilon:	maximum enzyme production rate (h^-1)
!
!				m_bf:			mortality free-living bacteria (h^-1)
!
!				m_ba:			mortality of attached bacteria (h^-1)
!	
!				swim:			swimming velocity of bacteria (m h^-1) (Kiorboe et al. 2002)
!
!				run:			time of bacterial swimming interval (h) (Kiorboe et al. 2002)
!
!				angle:			angle between bacteria runs (cosine of angle)
!	
!				max_detach:		maximum detachment rate (h^-1) (Kiorboe et al. 2002, 2003)
!								**observations of antagonistic interactions (Long and Azam 2001)
!
!
!       ENZYME PARAMETERS
!
!				ehl				half life of the enzyme at 7°C (h)	(Steen and Arnosti 2011)
!
!				epom			maximum degradation rate of the particle by extracellular enzyme at 7°C 
!								(mg C particle (mg C extracellular enzyme)^-1) (Steen and Arnosti 2011)
!
!				f_CaCO3			calcium carbonate ballast preservation fraction mg C (mg CaCO3)^-1 (Klaas and Archer 2002)
!
!				f_opal			opal ballast preservation fraction mg C (mg opal)^-1 (Klaas and Archer 2002)
!
!				f_lith			lithogenic material ballast preservation fraction mg C (mg lithogenic)^-1 (Klaas and Archer 2002)
!			
!				flux_CaCO3		flux of calcium carbonate (mg CaCO3 m^-2 d^-1)
!
!				flux_opal		flux of opal (mg opal m^-2 d^-1)
!
!				flux_lith		flux of lithogenic material (mg lithogenic m^-2 d^-1)
!               
!               
!		DIFFUSION PARAMETERS
!
!				bc				Boltzmann constant (m^2 kg s^-2 K^-1)
!
!				mr_enz			molecular radius of a beta-glucosidase molecule (m) (Jeng et al. 2010)
!
!				mr_h			molecular radius of a glucose molecule (m) (Pappenheimer et al. 1951)
!
!
!       POM PARAMETERS
!   
!               epom:		rate of pom degradation through extracellular 
!							enzyme hydrolysis (per hr)
!
!               mr_h:		molecular radius of the hydrolysate
!
!				r:			characteristic radius of a particle in m
!
!				w:			flux rate of POM in m/hr
!
!
!       DOM PARAMETERS
!   
!               tt:		turnover time in years of semi-labile dom into labile dom (hdom)
!							
!
!	
!------------------------------------
!   FORTRAN VARIABLE INITIALIZATION
!------------------------------------

	CHARACTER (150)::	TimeForce
	CHARACTER (150)::	DepthForce
	REAL::	hdom_initial
	REAL::	pom_initial
	REAL::	bf_initial
	REAL::	ba_initial
	REAL::	enz_initial
	REAL::	denz_initial
	REAL::	h_initial
	REAL::	edom_initial
	REAL::	dedom_initial
	REAL::	q10_b
	REAL::	q10_enz
	REAL::	q10_hl
	REAL::	dz
	REAL::	dt
	REAL::	w
	REAL::	r_top
	REAL::	ba_poc
	REAL::	up_max_bf
	REAL::	up_max_ba
	REAL::	kh
	REAL::	kdom
	REAL::	bge_max
	REAL::	basal
	REAL::	quorum
	REAL::	max_epsilon
	REAL::	m_bf
	REAL::	m_ba
	REAL::	swim
	REAL::	run
	REAL::	angle
	REAL::	max_detach
	REAL::	ehl
	REAL::	epom
	REAL::	f_CaCO3
	REAL::	f_opal
	REAL::	f_lith
	REAL::	flux_CaCO3
	REAL::	flux_opal
	REAL::	flux_lith
	REAL::	bc
	REAL::	mr_enz
	REAL::	mr_h
	REAL::	cperb
	REAL::  bdiam
	REAL::  blength
	INTEGER::	OpenStatus1
	INTEGER::   inputstatus1
	INTEGER::	jday
	INTEGER::	hhmm
	REAL::	pom_in_top
	INTEGER::	OpenStatus10
	INTEGER::   inputstatus10
	INTEGER::	depth
	REAL::	temp
	REAL::	s
	REAL::	rho
	INTEGER:: tt
	REAL::	sldoc
	INTEGER:: t
	INTEGER::	MaxDepth = 99999	
	INTEGER, PARAMETER:: max_dim = 8000
	INTEGER, DIMENSION(max_dim) :: list_time
	REAL, DIMENSION(max_dim) :: list_pom = 0
	REAL, DIMENSION(max_dim) :: list_enz = 0
	REAL, DIMENSION(max_dim) :: list_h = 0
	REAL, DIMENSION(max_dim) :: list_ba = 0
	REAL, DIMENSION(max_dim) :: list_hdom = 0
	REAL, DIMENSION(max_dim) :: list_bf = 0 
	REAL, DIMENSION(max_dim) :: list_denz = 0
	REAL, DIMENSION(max_dim) :: list_edom = 0
	REAL, DIMENSION(max_dim) :: list_dedom = 0
	REAL, DIMENSION(max_dim) :: list_flux = 0
	REAL, DIMENSION(max_dim) :: list_attach = 0
	REAL, DIMENSION(max_dim) :: list_detach = 0
	REAL, DIMENSION(max_dim) :: list_ba_grow = 0
	REAL, DIMENSION(max_dim) :: list_bf_grow = 0
	REAL, DIMENSION(max_dim) :: list_epsilon = 0
	REAL, DIMENSION(max_dim) :: list_pnum = 0
	REAL, DIMENSION(max_dim) :: list_pom_dgd = 0
	REAL, DIMENSION(max_dim) :: list_h_hdom = 0
	REAL, DIMENSION(max_dim) :: list_ehl_k = 0
	REAL, DIMENSION(max_dim) :: list_enz_dom = 0
	REAL, DIMENSION(max_dim) :: list_denz_dom = 0
	REAL, DIMENSION(max_dim) :: list_tfacb = 0
	REAL, DIMENSION(max_dim) :: list_tfacenz = 0
	REAL, DIMENSION(max_dim) :: list_tfachl = 0
	REAL, DIMENSION(max_dim) :: list_b_ratio = 0
	REAL, DIMENSION(max_dim) :: list_bp_ratio = 0
	REAL, DIMENSION(max_dim) :: list_up_ba = 0
	REAL, DIMENSION(max_dim) :: list_up_bf = 0
	REAL::  decayK
	REAL::	ldoc
	REAL::	tfac_b
	REAL::  tfac_enz
	REAL::  tfac_hl
	REAL::	E
	REAL::	F
	REAL::	Clv
	REAL::	dv_fw
	REAL::	dv_sw
	REAL::  kv
	REAL::	r
	REAL::	dw
	REAL::  pnum
	REAL::  y
	REAL::  Vtotal
	REAL::  Vsurf
	REAL::	Re
	REAL::	up_bf
	REAL::  hp
	REAL::  up_ba
	REAL::  bapop
	REAL::  k_epsilon
	REAL::	epsilon
	REAL::	grow_bf
	REAL::	grow_ba
	REAL::	pomp
	REAL::	gamma
	REAL::	gammap
	REAL::	availC
	REAL::	pom_dgd
	REAL::	D
	REAL::	Sh
	REAL::	RW
	REAL::	attach
	REAL::  k_detach
	REAL::	detach
	REAL::	diff_h
	REAL::	sc_h
	REAL::	sh_h
	REAL::	h_hdom
	REAL::	diff_enz
	REAL::	sc_enz
	REAL::	sh_enz
	REAL::	enz_edom
	REAL::	denz_dedom
	REAL::	ehl_k
	REAL::	denzdecay
	REAL::  bacover
	REAL::	fpom_in
	REAL::	fba_in
	REAL::	fenz_in
	REAL::	fh_in
	REAL::	fdenz_in
	REAL::	fpom_out
	REAL::	fenz_out
	REAL::	fh_out
	REAL::	fba_out
	REAL::	fdenz_out
	REAL::  jm_ba
	REAL::	jpom_dgd
	REAL::	jh_hdom
	REAL::	jup_ba
	REAL::  hdiff
	REAL::  hplus
	REAL::	jgrow_ba
	REAL::	jattach
	REAL::	jdetach
	REAL::	jgrow_bf
	REAL::	jm_bf
	REAL::	jepsilon
	REAL::	jenz_edom
	REAL::	jm_enz
	REAL::	jup_bf
	REAL::	jdedomdecay
	REAL::	jdenzdecay
	REAL::	jdenz_dedom
	REAL::	jm_edom
	REAL::	d_pom
	REAL::	d_h
	REAL::	d_ba
	REAL::	d_bf
	REAL::	d_enz 
	REAL::	d_hdom
	REAL::	d_denz
	REAL::	d_edom
	REAL::	d_dedom
	REAL::	b_ratio
	REAL::	bp_ratio
	INTEGER:: i
    REAL:: hdom_check
	REAL::  hcheck
	REAL::  hcheck2
	REAL::  qtest
	INTEGER:: finalt
	REAL::  qfrac
    REAL::  rho_Csub
    REAL::  Vcarbon
    REAL::  qk
    REAL::  qthresh
    REAL::  havail
    INTEGER::  scenario
    REAL::  qxfit
    INTEGER:: etype
    REAL::  Ratio_mmC

!----------------------------------
!			NAMELIST
!----------------------------------

NAMELIST /control/ qfrac, max_epsilon

!----------------------------------
!		SCENARIO CHOICE
!----------------------------------
scenario = 4

!scenario = 1  (Interior Arrangement)
!scenario = 3  (Solute Retention)
!scenario = 4  (Hydrolysate Interception)

!----------------------------------
!		ENZYME CHOICE
!----------------------------------
etype = 2

!etype = 1  (Beta-glucosidase)
!etype = 2  (Leucine aminopeptidase)

!----------------------------------
!	  DEFAULT PARAMETER VALUES
!----------------------------------
!		******************************************
		hdom_initial	= 1.5			!mg C m^-3 seawater  
		pom_initial		= 0             !mg C m^-3 seawater
		bf_initial		= 1.5			!mg C m^-3 seawater
		ba_initial		= 0				!mg C m^-3 seawater 
		enz_initial		= 0				!mg C m^-3 seawater
		denz_initial	= 0				!mg C m^-3 seawater
		h_initial		= 0				!mg C m^-3 seawater
		edom_initial	= 0				!mg C m^-3 seawater
		dedom_initial	= 0				!mg C m^-3 seawater
!		******************************************
		q10_b			= 2				!
		q10_enz			= 2				!
		q10_hl			= 2				!
!		******************************************
		dz				= 10			!m
        dt				= 0.00833		!h                          (30 s)
!		******************************************
		w				= 2				!m h^-1						(~50 m/day)
		r_top			= 0.001			!m							(1 mm)
!		******************************************				
		up_max_bf		= 0.25			!h^-1
		up_max_ba		= 0.25			!h^-1
        kh				= 200           !5 mg C m^-3 particle
		kdom			= 30            !20 mg C m^-3 seawater
		bge_max			= 0.4			!
		basal			= 0.002			!h^-1 at 20°C
		max_epsilon		= 0.028         !h^-1
		m_bf			= 0.02			!h^-1
		m_ba			= 0.04			!h^-1
		swim			= 0.09			!m h^-1						(ranges from 0.076-0.12 m h^-1) (Kiorboe et al. 2002)
		run				= 0.002			!h							(ranges from 0.0001-0.004 h) (Kiorboe et al. 2002)
		angle			= -0.67			!							(Kiorboe et al. 2002)
!		******************************************				
		ehl				= 30			!30 at 20°C; 72 at 7°C		(Steen and Arnosti 2011)
		f_CaCO3			= 0.075			!mg C (mg CaCO3)^-1			(Klaas and Archer 2002)
		f_opal			= 0.029			!mg C (mg opal)^-1			(Klaas and Archer 2002)
		f_lith			= 0.052			!mg C (mg lithogenic)^-1	(Klaas and Archer 2002)
        flux_CaCO3      = 21			!mg CaCO3 m^-2 d^-1			(Conte et al. 2001)
		flux_opal       = 4.34			!mg opal m^-2 d^-1			(Conte et al. 2001)
        flux_lith       = 6.77			!mg lith m^-2 d^-1			(Conte et al. 2001)
!		******************************************
		bc				= 1.38E-23		!m^2 kg s^-2 K^-1
!		******************************************				
		cperb			= 5E-12			!mg C (bacteria)^-1
		bdiam			= 5E-7			!m							(Kirchman 2012)
		blength			= 1E-6			!m							(Kirchman 2012)
!		******************************************		
		tt				= 5.			!years						(Abell et al. 2000)
!		******************************************		
		temp			= 20            !degrees C, 19.19 BATS depth 150 m
		s				= 36.632		!psu, BATS depth 150 m
		rho				= 1026.226		!density, BATS depth 150 m
		sldoc			= 174.77		!doc, BATS depth 150 m
		pom_in_top		= 1.4			!mg C m^-3 h^-1
		qfrac			= 1.1
!		******************************************

!----------------------------------
!   SCENARIO SPECIFIC PARAMETERS
!----------------------------------

IF(etype .eq. 1)THEN
!Beta-glucosidase
epom            = 8.6           !at 20°C, q10=2 (Steen and Arnosti 2011)
rho_Csub        = 1.54E9        !mg m^-3 density of glucose
mr_enz			= 6.4E-9		!m glucosidase	(Jeng et al. 2010)
mr_h			= 3.8E-10		!m glucose	(Pappenheimer et al. 1951)
Ratio_mmC       = 180/72        !dimensionless, ratio of molar mass glucose to molar mass of carbon in glucose
ENDIF

IF(etype .eq. 2)THEN
!Leucine aminopeptidase
epom            = 118           !at 20°C, q10=2 (Steen and Arnosti 2011)
rho_Csub        = 1.29E9        !mg m^-3 density of leucine
mr_enz			= 4.9E-9        !Leucine aminopeptidase	(Burley et al. 1990)
mr_h            = 1.5E-10       !m Leucine (Rawn 1989)
Ratio_mmC       = 131/72        !dimensionless, ratio of molar mass Leucine to molar mass of carbon in Leucine
ENDIF

PRINT*, "scenario = ", scenario

!----------------------------------
!		NAMELIST PARAMETERS
!----------------------------------

READ (5, control) !	Read in the namelist
WRITE (6, control)
PRINT *, "------------------------------------"


!----------------------------------
!		MAKE OUTPUT FILES
!----------------------------------

OPEN(UNIT=2, FILE="StateVariables.out", STATUS='REPLACE')
OPEN(UNIT=3, FILE="BacteriaRates.out", STATUS='REPLACE')
OPEN(UNIT=4, FILE="MassTransferRates.out", STATUS='REPLACE')

!----------------------------------
!	  INITIAL STATE VARIABLES
!----------------------------------
	list_hdom(1) = hdom_initial
	list_bf(1) = bf_initial
	!list_enz(1) = enz_initial
	list_denz(1) = denz_initial
	list_h(1) = h_initial
	list_edom(1) = edom_initial
	list_dedom(1) = dedom_initial
	list_time(1) = 0
	
    finalt = 600

!----------------------------
!   WATER COLUMN PROFILE
!----------------------------

! Calculate the amount of semi-labile doc converted to labile doc per time step

	decayK = (1./(tt*24.*365.))

!----------------------------
!  TEMPERATURE DEPENDENCE
!----------------------------	

! Calculate temperature scaling for each rate
	
	tfac_b =	q10_b**((temp-20)/10)

	tfac_enz =	q10_enz**((temp-20)/10)
	
	tfac_hl =	q10_hl**((20-temp)/10)


!----------------------------
!  SEAWATER CHARACTERISTICS
!----------------------------	

! Calculate dynamic viscosity and kinematic viscosity of seawater (equations A1, A7, A8, A9, A10 from Jumars et al. 1993)
	
	E = 0.0001068 + 0.00005185*temp	!A8
	
	F = 0.002591 + 0.000033*temp	!A9
	
	Clv = ((rho/1000)*s)/1.80655	!A10  density converted to g cm**-3
	
	dv_fw = (10**((-157.095 - 3.09695*temp - 0.001827*temp**2)/(89.93+temp))) !A1 
	
	dv_sw = (dv_fw * (1 + E*SQRT(Clv) + F*Clv))*0.1	!A7 dynamic viscosity in kg m^-1 s^-1
	
	kv = (dv_sw/rho)*3600 !kinematic viscosity in m^2 h^-1 


!----------------------------
!  PARTICLE CHARACTERISTICS
!----------------------------

    r = r_top
	
    dw = (17*(r*2*1000)**1.8)  !Iverson et al. 2010 (r*1000 is converting m to mm)

    pnum = 84
	
    pom_initial = pnum*(dw*0.001*0.12)  !0.12 is the ratio of POC to dry weight (dw) from Iverson et al. 2010 (dw*0.001 is converting ug to mg)

    list_pom(1) = pom_initial
	
!Calculate particle volume in order to scale the constants for quorum sensing and hydrolysate up
	
	Vtotal = (4./3.)*3.14*r**3
	
! Reynolds Number of a particle

	Re = (w*r)/kv


! Calculation of Quorum Threshold  (Gram et al. 2002)
	
	quorum = ((4*3.14*r**2.)/(bdiam*blength))*cperb  ! mg C particle^-1  One layer of bacteria covering the particle

	list_ba(1) = quorum*pnum*qfrac
	list_enz(1) = 0  


DO t = 2, finalt
!------------------------
!    RATE CONSTANTS
!------------------------	

! Calculate attached bacteria covering the surface area of a particle

    IF(scenario .eq. 3)THEN

        bacover = (((list_ba(t-1)/pnum)/cperb)*(bdiam*blength))/(4*3.14*r**2.)

        IF (bacover .gt. 1)THEN

            bacover = 1

        ENDIF

    ELSE

        bacover = 0

    ENDIF

! Particle degradation

!	amount of pom degraded which is dependent on the amount of pom and temperature
	
	pomp = list_pom(t-1)/(pnum*Vtotal)
	
	gamma = ((f_CaCO3*flux_CaCO3 + f_opal*flux_opal + f_lith*flux_lith)/24/w)!ballast scaling constant for enzyme degradation (Klaas and Archer 2002)
	
	gammap = gamma/(pnum*Vtotal)
	
	availC = (pomp-gammap)/pomp

	IF (availC .lt. 0)THEN
	
	availC = 0
	
	ENDIF
	
	pom_dgd = tfac_enz*epom*availC

! Calculate Volume of Substrate Carbon Available Per Particle

    Vcarbon = (((list_pom(t-1)*availC)/pnum)*(Ratio_mmC))/rho_Csub

! Bacteria attachment to particles

!	rate constant using equation 6a from Kiorboe et al. 2002

	D = tfac_b*((swim*swim*run)/(6*(1-angle)))!(Diffusivity)
	
	Sh = 1+0.619*(Re**0.412)*(kv/D)**(1./3.) !(Sherwood Number)
			
	RW = 4*3.14*D*r*Sh !Encounter Rate Kernel for Random Walk 

	attach = 0 !RW*pnum
	
! Bacteria detachment from particles
	
	k_detach = (quorum)
	
	detach = 0 !max_detach*(bapop/(k_detach+bapop))

! Hydrolysate and extracellular enzyme leaving sinking pom

!	hydrolysate rate constant using equation 4 from Kiorboe et al. 2001

	diff_h = ((bc*(temp+273.16))/(6*3.14*dv_sw*mr_h))*3600	!equation 6 from Jumars et al. 1993 converted to h^-1
	
	sc_h = kv/diff_h !Schmidt Number
	
	sh_h = 1 + 0.619*(Re**0.412)*(sc_h)**(1./3.) !(Sherwood Number)
	
	h_hdom = (4*3.14*diff_h*r*sh_h) !h^-1 Kiorboe et al. 2001

!	active extracellular enzyme rate constant using equation 4 from Kiorboe et al. 2001
	
	diff_enz = ((bc*(temp+273.16))/(6*3.14*dv_sw*mr_enz))*3600	!equation 6 from Jumars et al. 1993 converted to h^-1

	sc_enz = kv/diff_enz !Schmidt Number
	
	sh_enz = 1 + 0.619*(Re**0.412)*(sc_enz)**(1./3.) !(Sherwood Number)
	
	enz_edom = (4*3.14*diff_enz*r*sh_enz)  !h^-1 Kiorboe et al. 2001

!	inactive extracellular enzyme rate constant using equation 4 from Kiorboe et al. 2001
	
	denz_dedom = (4*3.14*diff_enz*r*sh_enz) !h^-1 Kiorboe et al. 2001

! Bacterial uptake and growth

!	uptake rate constant of free-living bacteria which is dependent on the amount of dom available and temperature

    up_bf = tfac_b*up_max_bf*(list_hdom(t-1)/(kdom+list_hdom(t-1)))


    hp = list_h(t-1)/(pnum*Vtotal)

    up_ba = tfac_b*up_max_ba*(hp/(kh+hp))

!	enzyme production constant dependent on ba biomass per particle - quorum sensing

    bapop = list_ba(t-1)/(pnum)

    epsilon = max_epsilon

!	calculate bacteria growth rate constants

    grow_bf = (bge_max)*up_bf - tfac_b*basal

    grow_ba = (bge_max)*up_ba - tfac_b*basal - epsilon

! Extracellular enzyme halflife

!	enzyme half-life
	ehl_k = (LOG(0.5)/(tfac_hl*ehl))*(-1)

!------------------------
!        FLUXES
!------------------------

!Sinking Rates In
	fpom_in = 0 !pom_in_top/dz
	fba_in = 0 !ba_poc*fpom_in !See Table 1 in Ghighlione et al. 2009, 1E11 is attached bacteria at 150 m and 26.4 is pom at 150 m 
	fenz_in = 0 
	fh_in = 0
	fdenz_in = 0 

!Sinking Rates Out
	fpom_out = 0 !(w*list_pom(t-1))/dz			!pom loss to sinking
	fenz_out = 0 !(w*list_enz(t-1))/dz			!enz loss to sinking
	fh_out = 0 !(w*list_h(t-1))/dz				!h loss to sinking
	fba_out = 0 !(w*list_ba(t-1))/dz				!ba loss to sinking
	fdenz_out = 0 !(w*list_denz(t-1))/dz			!dead enzyme loss to sinking

!------------------------
!        RATES
!------------------------

	jm_ba = m_ba*list_ba(t-1)**2				!mortality of attached bacteria
	jpom_dgd = pom_dgd*list_enz(t-1)			!degradation of pom by extracellular enzyme
	jup_ba = up_ba*list_ba(t-1)					!uptake by attached bacteria

    IF(scenario .eq. 4)THEN

        hdom_check = h_hdom*((list_h(t-1)/(Vtotal*pnum)))*pnum - up_ba*list_ba(t-1)!Hydrolysate interception

        IF (hdom_check .lt. 0)THEN

            jh_hdom = 0

        ELSE

            jh_hdom = hdom_check

        ENDIF

    ELSE

        jh_hdom = (1-bacover)*h_hdom*((list_h(t-1)/(Vtotal*pnum)))*pnum     !hydrolysate flux out

    ENDIF 

	jgrow_ba = grow_ba*list_ba(t-1)			!growth rate of attached bacteria	
	jattach = attach*list_bf(t-1)*pnum		!attachment rate of bacteria
	jdetach = detach*list_ba(t-1)				!detachment rate of bacteria
	jgrow_bf = grow_bf*list_bf(t-1)			!growth rate of free-living bacteria
	jm_bf = m_bf*list_bf(t-1)**2				!mortality rate of free-living bacteria
	
	jepsilon = epsilon*list_ba(t-1)			!enzyme production rate by attached bacteria

    IF(scenario .eq. 3)THEN

        jenz_edom = (1-bacover)*enz_edom*((list_enz(t-1)-list_enz(t-1)*(Vcarbon/Vtotal))/Vtotal) !active enzyme flux

    ELSE

        jenz_edom = enz_edom*((list_enz(t-1)-list_enz(t-1)*(Vcarbon/Vtotal))/Vtotal)

    ENDIF
    
    jm_enz = ehl_k*list_enz(t-1)				!enzyme deactivation rate in the particle
    jup_bf = up_bf*list_bf(t-1)                 !uptake by free-living bacteria
	jdedomdecay = decayK*list_dedom(t-1)		!inactive enzyme decay rate to hdom
	jdenzdecay = decayK*list_denz(t-1)

    IF(scenario .eq. 3)THEN

        jdenz_dedom = (1-bacover)*denz_dedom*((list_denz(t-1)/(Vtotal*pnum)))*pnum !inactive enzyme flux 

    ELSE

        jdenz_dedom = denz_dedom*((list_denz(t-1)/(Vtotal*pnum)))*pnum !inactive enzyme flux 

    ENDIF

	jm_edom = ehl_k*list_edom(t-1)			!enzyme deactivation rate in the dissolved environment

!-----------------------------
!   CHANGE IN CONCENTRATION
!-----------------------------

	d_pom = (fpom_in + jm_ba - jpom_dgd - fpom_out)*dt
	
	d_h = (fh_in + jpom_dgd + jdenzdecay - jh_hdom - jup_ba - fh_out)*dt

	d_ba = (fba_in + jgrow_ba + jattach - jdetach - jm_ba - fba_out)*dt
	
	d_bf = (jgrow_bf - jattach + jdetach - jm_bf)*dt
			
	d_enz = (fenz_in + jepsilon - jenz_edom - jm_enz - fenz_out)*dt

	d_hdom = (decayK*sldoc + jh_hdom + jdedomdecay + jm_bf - jup_bf)*dt 
	
	d_denz = (fdenz_in + jm_enz - jdenz_dedom - jdenzdecay - fdenz_out)*dt
	
	d_edom = (jenz_edom - jm_edom)*dt
	
	d_dedom = (jm_edom + jdenz_dedom - jdedomdecay)*dt


!---------------------------
!     NEW CONCENTRATIONS
!---------------------------

	list_time(t) = (t-1)*60
	list_pom(t) =	list_pom(t-1) + d_pom
	list_enz(t) =	list_enz(t-1) + d_enz
	list_h(t)	=	list_h(t-1) + d_h
	list_ba(t)	=	list_ba(t-1) + d_ba
	list_hdom(t) =	list_hdom(t-1) + d_hdom
	list_bf(t)	=	list_bf(t-1) + d_bf
	list_denz(t) = list_denz(t-1) + d_denz
	list_edom(t) = list_edom(t-1) + d_edom
	list_dedom(t) = list_dedom(t-1) + d_dedom 


!---------------------------
!      SAVE RATES
!---------------------------
	
	!Bacteria
	list_attach(t) = attach 
	list_detach(t) = detach
	list_ba_grow(t) = grow_ba
	list_bf_grow(t) = grow_bf
	list_up_ba(t) = up_ba
	list_up_bf(t) = up_bf
	list_epsilon(t) = epsilon
	
	!Mass Transfer
	list_pnum(t) = pnum 
	list_pom_dgd(t) = pom_dgd 
	list_h_hdom(t) = h_hdom
	list_ehl_k(t) = ehl_k
	list_enz_dom(t) = enz_edom 
	list_denz_dom(t) = denz_dedom
	

!---------------------------
!      ERROR CHECKING
!---------------------------

    IF(list_h(t) .lt. 0)THEN

        PRINT *, "***ERROR Negative Hydrolysate***"
        PRINT *, "CountTime =", t

    STOP
    END IF

!-----------------------------------
!   PRINT FINAL OUTPUT TO FILES
!-----------------------------------

	IF (t == finalt) THEN
		
			!PRINT *, "Number of Time Steps = ", finalt
			
			DO i=1, finalt

				WRITE(2,5) list_time(i), list_pom(i), list_enz(i), list_h(i), list_ba(i), list_hdom(i), list_bf(i),&
				 list_edom(i), list_denz(i), list_dedom(i)
				5 FORMAT(I8,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5)

				WRITE(3,7) list_time(i), list_attach(i), list_detach(i), list_ba_grow(i), list_bf_grow(i), list_up_ba(i), list_up_bf(i),&
				list_epsilon(i)
				7 FORMAT(I8,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5)
		
				WRITE(4,8) list_time(i), list_pnum(i), list_pom_dgd(i), list_h_hdom(i), list_ehl_k(i), list_enz_dom(i), list_denz_dom(i) 
				8 FORMAT(I8,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5)
				
		
			END DO
			
			CLOSE(1)
			STOP
	
	ENDIF

!-------------------------------
!   PRINT PROGRESS TO TERMINAL
!-------------------------------

!PRINT *, "Time Iteration Completed = ", t
!PRINT *, "------------------------------------"

END DO

END PROGRAM MAIN