PROGRAM MAIN
IMPLICIT NONE
!		"1D model of of sinking particle remineralization in the ocean water column that includes bacterial dynamics"
!
!----------------------------------
!  DESCRIPTION OF VARIABLES
!----------------------------------
!
!   SCENARIO
!
!       scenario:       choices are:
!                       (1) Interior - bacterial uptake is separate from diffusive flux from particles
!                       (2) Interception - bacterial uptake is from the diffusive flux from particles
!                       (3) Retention - exoenzyme and hydrolysate flux is stopped by particle-attached bacteria
!
!   FLUX IN
!
!		jday:			julian day
!				
!		hhmm:			hour, minute
!				
!		pom_in_top:		particle flux at upper boundary (mg C m^-2)
!
!
!   WATER COLUMN PROFILE
!
!		depth:			depth in the water column
!
!		temp:			temperature of seawater (°C)
!
!		s:				salinity of seawater (psu)
!
!		rho:			density of seawater (kg m^-3)
!  
!		sldoc:			semi-labile dissolved organic carbon (mg C m^-3)
!
!
!	TIME AND SPACE
!
!		dz:				depth interval (m)
!
!       k:              depth iteration
!
!		dt:				time step (h)
!
!       t:              time iteration
!
!               
!   STATE VARIABLES
!              
!       hdom:			hydrolysate dissolved organic matter (mg C m^-3)
!           
!       pom:			particulate organic matter (mg C m^-3)
!                       
!       bf:				free living bacteria (mg C m^-3)
!
!       ba:				attached bacteria (mg C m^-3)
!
!       enz:			exoenzyme in particle (mg C m^-3)
!	
!       h:				hydrolysate in particle (mg C m^-3)
!
!       denz:			inactive exoenzyme in particle (mg C m^-3)
!
!       edom:			active exoenzyme in dissolved environment (mg C m^-3)
!
!       dedom:			inactive exoenzyme in dissolved environment (mg C m^-3)
!
!
!   TEMPERATURE DEPENDENCE
!
!       q10_b:			q10 value for bacterial temperature dependence.
!							
!       q10_enz:		q10 value for exoenzyme activity rate.
!		
!       q10_hl:			q10 value for exoenzyme half-life.
!
!
!   PARTICULATE ORGANIC MATTER PARAMETERS
!
!       w:				sinking speed of the particles (m h^-1)
!
!       r_top:			radius of the particle at the upper boundary (m)
!
!
!   BACTERIA PARAMETERS
!
!       up_max_bf:		maximum free-living bacteria uptake rate at 20°C (h^-1)
!
!       up_max_ba:		maximum attached bacteria uptake rate at 20°C (h^-1)
!							
!       kh:				half-saturation constant for attached bacterial uptake of hydrolysate (mg C m^-3 particle)
!
!       kdom:			half-saturation constant for free-living bacterial uptake of dissolved hydrolysate (mg C m^-3 seawater)
!
!       bge_max:		maximum bacterial growth efficiency (Fig. 3, Del Giorgio and Cole 2000)
!
!       basal:			basal energetic cost at 20°C (h^-1) (based on the 0.01 day^-1 at 0°C from COBALT v1, Stock et al. in prep)
!		
!       quorum:			bacterial biomass at quorum sensing threshold (mg C m^-3 particle)
!		
!       qfrac:          fraction of the surface area of the particle covered by particle-attached bacteria
!
!       x:				exponent for slope of the quorum threshold response
!
!       max_epsilon:	maximum exoenzyme production rate (h^-1)
!
!       m_bf:			mortality free-living bacteria (h^-1)
!
!       m_ba:			mortality of attached bacteria (h^-1)
!	
!       cperb:          mg carbon per bacteria (mg C (bacteria)^-1)
!
!       bdiam:          diameter of a rod-shaped bacteria (m)
!
!       blength:        length of a rod-shaped bacteria (m)
!
!       swim:			swimming velocity of bacteria (m h^-1) (Kiorboe et al. 2002)
!
!       run:			time of bacterial swimming interval (h) (Kiorboe et al. 2002)
!
!       angle:			angle between bacteria runs (cosine of angle) (Kiorboe et al. 2002)
!
!
!   EXOENZYME PARAMETERS
!
!       ehl:            half-life of exoenzyme activity (h)	(Steen and Arnosti 2011)
!
!       epom:           maximum exoenzyme activity rate (mg C particle (mg C exoenzyme)^-1) (Steen and Arnosti 2011)
!
!       f_CaCO3:        calcium carbonate ballast preservation fraction mg C (mg CaCO3)^-1 (Klaas and Archer 2002)
!
!       f_opal:         opal ballast preservation fraction mg C (mg opal)^-1 (Klaas and Archer 2002)
!
!       f_lith:         lithogenic material ballast preservation fraction mg C (mg lithogenic)^-1 (Klaas and Archer 2002)
!
!       flux_CaCO3:     flux of calcium carbonate (mg CaCO3 m^-2 d^-1)
!
!		flux_opal:      flux of opal (mg opal m^-2 d^-1)
!
!		flux_lith:      flux of lithogenic material (mg lithogenic m^-2 d^-1)
!               
!               
!   DIFFUSION PARAMETERS
!
!		bc:             Boltzmann constant (m^2 kg s^-2 K^-1)
!
!		mr_enz:         molecular radius of a exoenzyme leucine aminopeptidase molecule (m)
!
!		mr_h:           molecular radius of a hydrolysate leucine molecule (m)
!
!       rho_Csub:       density of leucine (mg m^-3)
!
!       Ratio_mmC:      ratio of molar mass Leucine to molar mass of carbon in Leucine (dimensionless)
!
!
!   DOM PARAMETERS
!   
!       tt:             turnover time in years of semi-labile dom into labile dom (hdom)
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
    INTEGER::	ss
	REAL::	pom_in_top
	INTEGER::	OpenStatus10
	INTEGER::   inputstatus10
	INTEGER::	depth
	REAL::	temp
	REAL::	s
	REAL::	rho
	INTEGER:: tt
	REAL::	sldoc
	INTEGER::	CountTime
	INTEGER::	k
	INTEGER::	MaxDepth = 99999	
	INTEGER, PARAMETER:: max_dim = 500
	INTEGER, DIMENSION(max_dim) :: list_depth
	REAL, DIMENSION(max_dim) :: list_pom = 0
	REAL, DIMENSION(max_dim) :: list_enz = 0
	REAL, DIMENSION(max_dim) :: list_h = 0
	REAL, DIMENSION(max_dim) :: list_ba = 0
	REAL, DIMENSION(max_dim) :: list_hdom = 0
	REAL, DIMENSION(max_dim) :: list_bf = 0 
	REAL, DIMENSION(max_dim) :: list_denz = 0
	REAL, DIMENSION(max_dim) :: list_edom = 0
	REAL, DIMENSION(max_dim) :: list_dedom = 0
	REAL, DIMENSION(max_dim) :: list_pomflux = 0
	REAL, DIMENSION(max_dim) :: list_baflux = 0
    REAL, DIMENSION(max_dim) :: list_enzflux = 0
    REAL, DIMENSION(max_dim) :: list_hflux = 0
    REAL, DIMENSION(max_dim) :: list_denzflux = 0
	REAL, DIMENSION(max_dim) :: list_attach = 0
	REAL, DIMENSION(max_dim) :: list_detach = 0
	REAL, DIMENSION(max_dim) :: list_ba_grow = 0
	REAL, DIMENSION(max_dim) :: list_bf_grow = 0
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
	REAL, DIMENSION(max_dim) :: list_bacover = 0
    INTEGER, DIMENSION(100) :: RetentionQfrac = 0
    REAL, DIMENSION(100) :: RetentionEpsilon = 0
    INTEGER:: QfracIndex = 0
    REAL :: EpsilonOpt = 0
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
    REAL::  hdom_check
	REAL::  hcheck
	REAL::  hcheck2
	REAL::  flusht
	REAL::  qfrac
	REAL::  tt_dedom
	REAL::  decayK2
    REAL::  rho_Csub
    REAL::  Vcarbon
    REAL::  qk
    REAL::  qthresh
    REAL::  qxfit
    INTEGER::  scenario
    REAL::  havail
    REAL::  Ratio_mmC
    INTEGER:: QfracFind
    INTEGER:: Index
    REAL:: Days

!----------------------------------
!			NAMELIST
!----------------------------------

NAMELIST /control/ w, r_top, qfrac

!----------------------------------
!		SCENARIO CHOICE
!----------------------------------

scenario = 2

!scenario = 1  (Interior)
!scenario = 2  (Interception)
!scenario = 3  (Retention)


!----------------------------------
!		PARAMETER VALUES
!----------------------------------
!******************************************
hdom_initial	= 1.5			!mg C m^-3 seawater
pom_initial		= 0.025			!mg C m^-3 seawater
bf_initial		= 1.5			!mg C m^-3 seawater
ba_initial		= 0.000001		!mg C m^-3 seawater
denz_initial	= 0.			!mg C m^-3 seawater
enz_initial		= 0.			!mg C m^-3 seawater
h_initial		= 0.			!mg C m^-3 seawater
edom_initial	= 0.			!mg C m^-3 seawater
dedom_initial	= 0.			!mg C m^-3 seawater
!******************************************
q10_b			= 2.			!
q10_enz			= 2.			!
q10_hl			= 2.			!
!******************************************
dz				= 10.			!m
dt				= 0.00833		!h (30 s)
!******************************************
up_max_bf		= 0.25			!h^-1
up_max_ba		= 0.25			!h^-1
kh				= 200.			!mg C m^-3 particle
kdom			= 30.			!mg C m^-3 seawater
bge_max			= 0.4			!
basal			= 0.002			!h^-1 at 20°C
m_bf			= 0.02			!h^-1
m_ba			= 0.04			!h^-1
swim			= 0.09			!m h^-1	(0.076-0.12 Kiorboe et al. 2002)
run				= 0.002			!h (0.0001-0.004 Kiorboe et al. 2002)
angle			= -0.67			!(Kiorboe et al. 2002)
!******************************************
ehl				= 30.			!72	(Steen and Arnosti 2011)
epom            = 118.          !at 20°C, q10=2 (Steen and Arnosti 2011)
mr_enz			= 4.9E-9        !Leucine aminopeptidase	(Burley et al. 1990)
rho_Csub        = 1.29E9        !mg m^-3 density of leucine
mr_h            = 1.5E-10       !m Leucine (Rawn 1989)
Ratio_mmC       = 131./72.      !dimensionless, ratio of molar mass Leucine to molar mass of carbon in Leucine
!******************************************
f_CaCO3			= 0.075			!mg C (mg CaCO3)^-1	(Klaas and Archer 2002)
f_opal			= 0.029			!mg C (mg opal)^-1	(Klaas and Archer 2002)
f_lith          = 0.052			!mg C (mg lithogenic)^-1 (Klaas and Archer 2002)
flux_CaCO3		= 21.			!mg CaCO3 m^-2 d^-1	(Conte et al. 2001)
flux_opal		= 4.34			!mg opal m^-2 d^-1 (Conte et al. 2001)
flux_lith		= 6.77			!mg lith m^-2 d^-1 (Conte et al. 2001)
!******************************************
bc				= 1.38E-23		!m^2 kg s^-2 K^-1
!******************************************
cperb			= 5.0E-12		!mg C (bacteria)^-1
bdiam			= 5.0E-7		!m (Kirchman 2012)
blength			= 1.0E-6		!m (Kirchman 2012)
!******************************************
tt				= 5.			!years (Abell et al. 2000)
tt_dedom		= 5.			!years
!******************************************
!w				= 2.			!set in namelist
!r_top			= 0.001			!set in namelist
!qfrac			= 0.50			!set in namelist
!******************************************


!----------------------------------
!   SCENARIO SPECIFIC PARAMETERS
!----------------------------------

!**Exoenzyme production parameters for scenario 1 (Interior)**

IF(scenario .eq. 1)THEN
max_epsilon     = 0.027
qk              = 0.52          !fractional particle coverage
qthresh         = 0.37          !fractional particle coverage
qxfit           = 1.42          !exponent for quorum sensing equation
ENDIF


!**Exoenzyme production parameters for scenario 2 (Interception)**

IF(scenario .eq. 2)THEN
max_epsilon     = 0.025
qk              = 0.43          !fractional particle coverage
qthresh         = 0.31          !fractional particle coverage
qxfit           = 1.44          !exponent for quorum sensing equation
ENDIF


!**Exoenzyme production from lookup table for scenario 3 (Retention)**

IF(scenario .eq. 3)THEN
OPEN(60, FILE = 'RetentionEpsilon.txt')
DO i=1, 100
READ(60, *) QfracIndex, EpsilonOpt
RetentionQfrac(i) = QfracIndex
RetentionEpsilon(i) = EpsilonOpt
ENDDO
CLOSE(60)
ENDIF


!----------------------------------
!		NAMELIST PARAMETERS
!----------------------------------
PRINT *, "------------------------------------"
READ (5, control) !	Read in the namelist
WRITE (6, control)
PRINT *, "------------------------------------"
PRINT*, "Scenario choice = ", scenario


!----------------------------------
!     READ IN FORCING FILES
!----------------------------------
	
OPEN(58, FILE = 'filename1.txt')
READ(58, '(A)') TimeForce
CLOSE(58)
OPEN(59, FILE = 'filename2.txt')
READ(59, '(A)') DepthForce
CLOSE(59)
OPEN(UNIT=1, FILE=TimeForce, STATUS='OLD', IOSTAT=OpenStatus1)
IF (OpenStatus1 > 0) STOP "*** Cannot open the file ***"


!----------------------------------
!       MAKE OUTPUT FILES
!----------------------------------

OPEN(UNIT=2, FILE="StateVariables.out", STATUS='REPLACE')
OPEN(UNIT=3, FILE="BacteriaRates.out", STATUS='REPLACE')
OPEN(UNIT=4, FILE="MassTransferRates.out", STATUS='REPLACE')

	
!----------------------------------
!  INITIAL STATE VARIABLES
!----------------------------------
	list_hdom = hdom_initial
	list_pom = pom_initial
	list_bf = bf_initial
	list_ba = ba_initial
	list_enz = enz_initial
	list_denz = denz_initial
	list_h = h_initial
	list_edom = edom_initial
	list_dedom = dedom_initial


!------------------
!  START TIME LOOP
!------------------

CountTime = 1

DO

!----------------------------------
! PARTICLE FLUX AT UPPER BOUNDARY
!----------------------------------

READ(1, *, IOSTAT = inputstatus1) jday, hhmm, ss, pom_in_top

IF (inputstatus1 > 0) STOP "*** Input error ***"
IF (inputstatus1 < 0) EXIT 

	
OPEN(UNIT=10, FILE=DepthForce, STATUS='OLD', IOSTAT=OpenStatus10)
IF (OpenStatus10 > 0) STOP "*** Cannot open the file ***"


!------------------
!  START DEPTH LOOP
!------------------

k = 1

DO 


!----------------------------
!  WATER COLUMN PROFILE
!----------------------------

READ(10, *, IOSTAT = inputstatus10) depth, temp, s, rho, sldoc 

IF (inputstatus10 > 0) STOP "*** Input error ***"
IF (inputstatus10 < 0) EXIT 


!**Decay of semi-labile to labile DOC per time step**

	decayK = (1./(tt*24.*365.))
	
	decayK2 = (1./(tt_dedom*24.*365.))


!----------------------------
!  TEMPERATURE DEPENDENCE
!----------------------------
	
    tfac_b =	q10_b**((temp-20.)/10.)

    tfac_enz =	q10_enz**((temp-20.)/10.)

    tfac_hl =	q10_hl**((20.-temp)/10.)

!----------------------------
!  SEAWATER CHARACTERISTICS
!----------------------------	

!**Calculate dynamic viscosity**
	
	E = 0.0001068 + 0.00005185*temp	!(A8, Jumars et al. 1993)
	
	F = 0.002591 + 0.000033*temp	!(A9, Jumars et al. 1993)
	
	Clv = ((rho/1000.)*s)/1.80655	!(A10, Jumars et al. 1993)
        !density converted to g cm**-3
	
	dv_fw = (10.**((-157.095 - 3.09695*temp - 0.001827*temp**2.)/(89.93+temp))) !(A1, Jumars et al. 1993)
	
	dv_sw = (dv_fw * (1. + E*SQRT(Clv) + F*Clv))*0.1	!(A7, Jumars et al. 1993)


!**Calculate kinematic viscosity**

	kv = (dv_sw/rho)*3600. !kinematic viscosity in m^2 h^-1


!----------------------------
!  PARTICLE CHARACTERISTICS
!----------------------------	

IF (k == 1)THEN

!**Particle radius**

	r = r_top


!**Particle dry weight**
	
	dw = (17.*(r*2.*1000.)**1.8)  !Iverson et al. 2010 (r*1000 is converting m to mm)


!**Initial particulate number**
	
    pnum = list_pom(k)/(dw*0.001*0.12) !(Iverson et al. 2010)
        !0.12 is the ratio of POC to dry weight (dw)
        !dw*0.001 is converting ug to mg

	
END IF
	
!**Particle volume**

	Vtotal = (4./3.)*3.14*r**3.


!**Reynolds Number of a particle**

	Re = (w*r)/kv


!**Quorum-sensing threshold**

    quorum = ((4.*3.14*r**2.)/(bdiam*blength))*cperb  ! mg C particle^-1
        !One layer of bacteria covering the surface area spherical particle


!------------------------
!    RATE CONSTANTS
!------------------------		

!**Fraction of particle surface area covered by attached bacteria**

    IF(scenario .eq. 3)THEN

        bacover = (((list_ba(k)/pnum)/cperb)*(bdiam*blength))/(4.*3.14*r**2.)

            IF (bacover .gt. 1.)THEN

                bacover = 1.

            ENDIF
    ELSE

        bacover = 0.

    ENDIF


!**Particle degradation**
	
	pomp = list_pom(k)/(pnum*Vtotal)

	gamma = ((f_CaCO3*flux_CaCO3 + f_opal*flux_opal + f_lith*flux_lith)/24./w)!(Klaas and Archer 2002)

	gammap = gamma/(pnum*Vtotal)
	
	availC = (pomp-gammap)/pomp
	
	IF (availC .lt. 0.)THEN
	
	availC = 0.
	
	ENDIF
	
	pom_dgd = tfac_enz*epom*availC


!**Volume of substrate carbon available per particle**

    Vcarbon = (((list_pom(k)*availC)/pnum)*Ratio_mmC)/rho_Csub


!**Bacterial attachment to particles**

	D = tfac_b*((swim*swim*run)/(6.*(1.-angle)))!(Diffusivity)
	
	Sh = 1.+0.619*(Re**0.412)*(kv/D)**(1./3.) !Sherwood Number
			
	RW = 4.*3.14*D*r*Sh !(Equation 6a, Kiorboe et al. 2002)

	attach = RW*pnum


!**Mass transfer of hydrolysate from particles**

	diff_h = ((bc*(temp+273.16))/(6.*3.14*dv_sw*mr_h))*3600. !h^-1 (Equation 6, Jumars et al. 1993)
	
	sc_h = kv/diff_h !Schmidt Number
	
	sh_h = 1. + 0.619*(Re**0.412)*(sc_h)**(1./3.) !Sherwood Number
	
	h_hdom = (4.*3.14*diff_h*r*sh_h) !h^-1 (Equation 4, Kiorboe et al. 2001)


!**Mass transfer of exoenzymes (active and inactive) from particles**
	
	diff_enz = ((bc*(temp+273.16))/(6.*3.14*dv_sw*mr_enz))*3600.  !h^-1 (Equation 6, Jumars et al. 1993)

	sc_enz = kv/diff_enz !Schmidt Number
	
	sh_enz = 1. + 0.619*(Re**0.412)*(sc_enz)**(1./3.) !Sherwood Number
	
	enz_edom = (4.*3.14*diff_enz*r*sh_enz)  !h^-1 (Equation 4, Kiorboe et al. 2001)
	
	denz_dedom = (4.*3.14*diff_enz*r*sh_enz) !h^-1 (Equation 4, Kiorboe et al. 2001)


!**Bacterial uptake**

    up_bf = tfac_b*up_max_bf*(list_hdom(k)/(kdom+list_hdom(k))) !free-living bacterial uptake

    hp = list_h(k)/(pnum*Vtotal)

    IF(scenario .eq. 3 .and. bacover .ge. 1)THEN

        kh = 500.  !higher kh because no hydrolysate is fluxing out of the particle

    ENDIF

    up_ba = tfac_b*up_max_ba*(hp/(kh+hp)) !particle-attached bacterial uptake


!**Exoenzyme production constant dependent on attached bacteria per particle**

    bapop = list_ba(k)/(pnum)

    IF(scenario .eq. 1 .or. scenario .eq. 2)THEN !Exoenzyme production from equation

        IF((bapop/quorum) .le. qthresh)THEN

            epsilon = 0.

        ELSE

            epsilon = max_epsilon*((bapop/quorum) - qthresh)**qxfit/((qk-qthresh)**qxfit + ((bapop/quorum) - qthresh)**qxfit)

        ENDIF

    ELSEIF(scenario .eq. 3)THEN  !Exoenzyme production from lookup table

        QfracFind = NINT((bapop/quorum)*100.)

        IF(QfracFind .gt. 100.)THEN

            epsilon = 0.0065

        ELSE

            Index = 0

            DO i = 1, 100

                IF(RetentionQfrac(i) .eq. QfracFind)THEN

                    Index = i

                    EXIT

                ENDIF

            ENDDO

            IF(Index .ne. 0)THEN

                epsilon = RetentionEpsilon(Index)

            ENDIF

            IF(Index .eq. 0)THEN

                epsilon = 0.

            ENDIF

        ENDIF

    ENDIF


!**Bacterial growth**

grow_bf = (bge_max)*up_bf - tfac_b*basal  !free-living bacterial growth

grow_ba = (bge_max)*up_ba - tfac_b*basal - epsilon !particle-attached bacterial growth


!**Bacterial detachment from particles**

    IF(grow_ba .gt. 0.)THEN

        IF(bapop .gt. quorum)THEN  !Maximum detachment rate dependent on density (i.e. competition)

            max_detach = grow_ba + (attach*list_bf(k))/(list_ba(k))

        ELSE

            max_detach = grow_ba

        END IF

    ELSE

        max_detach = 0.

    END IF

    k_detach = quorum

    detach = max_detach*(bapop/(k_detach+bapop))


!**Exoenzyme half-life**

	ehl_k = (LOG(0.5)/(tfac_hl*ehl))*(-1.)


!------------------------
!        FLUXES
!------------------------

!**Sinking Rates In**

IF(k == 1) THEN

	fpom_in = (pom_in_top)/dz
	fba_in = ((fpom_in)/(dw*0.001*0.12))*(quorum*qfrac)	
	fenz_in = 0.
	fh_in = 0.

ELSE
	
	fpom_in = w*list_pom(k-1)/dz			!pom fluxing in from above
	fenz_in = w*list_enz(k-1)/dz			!exoenzyme fluxing in from above
	fh_in = w*list_h(k-1)/dz				!hydrolysate fluxing in from above
	fba_in = w*list_ba(k-1)/dz				!attached bacteria fluxing in from above
	fdenz_in = w*list_denz(k-1)/dz			!dead exoenzyme fluxing in from above
	
END IF


!**Sinking Rates Out**

	fpom_out = (w*list_pom(k))/dz			!pom loss to sinking
	fenz_out = (w*list_enz(k))/dz			!exoenzyme loss to sinking
	fh_out = (w*list_h(k))/dz				!hydrolysate loss to sinking
	fba_out = (w*list_ba(k))/dz				!attached bacteria to sinking
	fdenz_out = (w*list_denz(k))/dz			!dead exoenzyme loss to sinking


!------------------------
!        RATES
!------------------------

	jm_ba = m_ba*list_ba(k)**2.				!mortality of attached bacteria
	jpom_dgd = pom_dgd*list_enz(k)			!degradation of pom by extracellular exoenzyme
    jup_ba = up_ba*list_ba(k)				!uptake by attached bacteria

    IF(scenario .eq. 2)THEN

        hdom_check = h_hdom*((list_h(k)/(Vtotal*pnum))-list_hdom(k))*pnum - up_ba*list_ba(k) !hydrolysate interception

        IF (hdom_check .lt. 0.)THEN

            jh_hdom = 0. !hydrolysate flux rate out

        ELSE

            jh_hdom = hdom_check !hydrolysate flux rate out

        ENDIF

    ELSE

        jh_hdom = (1.-bacover)*h_hdom*((list_h(k)/(Vtotal*pnum))-list_hdom(k))*pnum !hydrolysate flux rate out

    ENDIF

	jgrow_ba = grow_ba*list_ba(k)			!growth rate of attached bacteria
	jattach = attach*list_bf(k)				!attachment rate of bacteria
	jdetach = detach*list_ba(k)				!detachment rate of bacteria
	jgrow_bf = grow_bf*list_bf(k)			!growth rate of free-living bacteria
	jm_bf = m_bf*list_bf(k)**2.				!mortality rate of free-living bacteria
    jepsilon = epsilon*list_ba(k)			!exoenzyme production rate by attached bacteria

    IF(scenario .eq. 3)THEN

    jenz_edom = (1.-bacover)*enz_edom*(((list_enz(k)-list_enz(k)*(Vcarbon/Vtotal))/(Vtotal*pnum))-list_edom(k))*pnum  !active exoenzyme flux

    ELSE

    jenz_edom = enz_edom*(((list_enz(k)-list_enz(k)*(Vcarbon/Vtotal))/(Vtotal*pnum))-list_edom(k))*pnum  !active exoenzyme flux

    ENDIF

    jm_enz = ehl_k*list_enz(k)				!particle exoenzyme deactivation rate
	jup_bf = up_bf*list_bf(k)				!uptake by free-living bacteria
	jdedomdecay = decayK2*list_dedom(k)		!dissolved inactive exoenzyme decay rate
	jdenzdecay = decayK2*list_denz(k)       !particle inactive exoenzyme decay rate

    IF(scenario .eq. 3)THEN

	jdenz_dedom = (1.-bacover)*denz_dedom*((list_denz(k)/(Vtotal*pnum))-list_dedom(k))*pnum !inactive exoenzyme flux

    ELSE

    jdenz_dedom = denz_dedom*((list_denz(k)/(Vtotal*pnum))-list_dedom(k))*pnum !inactive exoenzyme flux

    ENDIF

    jm_edom = ehl_k*list_edom(k)  !dissolved exoenzyme deactivation rate


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


!-----------------------------------
!   BACTERIA EFFICIENCY CORRECTION
!-----------------------------------	

	hcheck = list_h(k) + d_h
	
	IF(hcheck .lt. 0.)THEN

	PRINT*, "efficiency used"
    PRINT *, "depth = ", depth
    PRINT *, "Time Iteration Completed = ", CountTime
    PRINT*, "-----------------------------"

	jup_ba = list_h(k)/dt + fh_in + jpom_dgd + jdenzdecay - jh_hdom - fh_out

    up_ba = jup_ba/list_ba(k)
	
	jgrow_ba = (bge_max*up_ba - tfac_b*basal - epsilon)*list_ba(k)

	d_h = (fh_in + jpom_dgd + jdenzdecay - jh_hdom - jup_ba - fh_out)*dt
	
	d_ba = (fba_in + jgrow_ba + jattach - jdetach - jm_ba - fba_out)*dt
	
	hcheck2 = list_h(k) + d_h  !In case there is a rounding error

    PRiNT*, "hcheck2 = ", hcheck2
	
	IF(hcheck2 .lt. 0.)THEN
	
	d_h = -list_h(k)  !Fix the rounding error
	
	END IF
		
	END IF
		

!---------------------------
!     NEW CONCENTRATIONS
!---------------------------

	list_depth(k) = depth
	list_pom(k) =	list_pom(k) + d_pom
	list_enz(k) =	list_enz(k) + d_enz
	list_h(k)	=	list_h(k) + d_h
	list_ba(k)	=	list_ba(k) + d_ba
	list_hdom(k) =	list_hdom(k) + d_hdom
	list_bf(k)	=	list_bf(k) + d_bf
	list_denz(k) = list_denz(k) + d_denz
	list_edom(k) = list_edom(k) + d_edom
	list_dedom(k) = list_dedom(k) + d_dedom 


!---------------------------
!       METRICS
!---------------------------

!**Ratio of attached bacteria to free-living bacteria**

    b_ratio = list_ba(k)/list_bf(k)


!**Ratio of attached bacterial production to free-living bacterial production**

    bp_ratio = ((grow_ba*list_ba(k))/(list_ba(k)/cperb))/((grow_bf*list_bf(k))/(list_bf(k)/cperb))



!**Flux of carbon from all sources**

IF (k == 1) THEN	

	list_pomflux(k) = pom_in_top
	list_baflux(k) = ((pom_in_top)/(dw*0.001*0.12))*(quorum*qfrac)
    list_enzflux(k) = 0.
    list_hflux(k) = 0.
    list_denzflux(k) = 0.

ELSE
	
	list_pomflux(k) = w*list_pom(k-1)
	list_baflux(k) = w*list_ba(k-1)
    list_enzflux(k) = w*list_enz(k-1)
    list_hflux(k) = w*list_h(k-1)
    list_denzflux(k) = w*list_denz(k-1)

END IF


!---------------------------
!      SAVE RATES
!---------------------------

	
	!**Bacteria**
	list_attach(k) = attach 
	list_detach(k) = detach
	list_ba_grow(k) = grow_ba
	list_bf_grow(k) = grow_bf
	list_b_ratio(k) = b_ratio 
	list_bp_ratio(k) = bp_ratio
	list_up_ba(k) = up_ba
	list_up_bf(k) = up_bf
	
	!**Mass transfer**
	list_pnum(k) = pnum 
	list_pom_dgd(k) = pom_dgd 
	list_h_hdom(k) = h_hdom
	list_ehl_k(k) = ehl_k
	list_enz_dom(k) = enz_edom 
	list_denz_dom(k) = denz_dedom
	list_bacover(k) = bacover
	
	!**Temperature scaling**
	list_tfacb(k) = tfac_b
	list_tfacenz(k) = tfac_enz
	list_tfachl(k) = tfac_hl
	
	
!---------------------------
!      ERROR CHECKING
!---------------------------

	IF(list_pom(k) .lt. 0.)THEN
	 
		PRINT *, "***ERROR Negative POM***"
		PRINT *, "CountTime =", CountTime
		PRINT *, "Depth = ", depth
	
		STOP
	END IF
	
	IF(list_h(k) .lt. 0.)THEN
	 
		PRINT *, "***ERROR Negative Hydrolysate***"
		PRINT *, "CountTime =", CountTime
		PRINT *, "Depth = ", depth
        PRINT *, "bacover = ", bacover
	
		STOP
	END IF
	
	
	IF(list_ba(k) .lt. 0.)THEN
	 
		PRINT *, "***ERROR Negative Attached Bacteria***"
		PRINT*, "list_ba(k) = ", list_ba(k)
		PRINT *, "CountTime =", CountTime
		PRINT *, "Depth = ", depth
	
		STOP
	END IF


    IF(Vcarbon/Vtotal .gt. 1.)THEN

        PRINT *, "***ERROR Carbon Volume Larger Than Particle Volume***"
        PRINT*, "Vcarbon/Vtotal = ", Vcarbon/Vtotal
        PRINT *, "CountTime =", CountTime
        PRINT *, "Depth = ", depth

        STOP
    END IF


!-----------------------------------
!   PRINT FINAL OUTPUT TO FILES
!-----------------------------------

	IF (CountTime == 43200 .AND. k == MaxDepth) THEN  !518400 (180 days); 259200 (90 days); 172800 (60 days); 86400 (30 days); 43200 (15 days)

            Days = (CountTime*30.)/60./60./24.

            PRINT*, "Days of Simulation = ", Days

            IF(Days .lt. 180.) THEN

                PRINT*, "**WARNING: 180 days needed for steady state**"
                PRINT*, "Change simulation time in the PRINT FINAL OUTPUT TO FILES section of MRM1D_v1.0.F90"

            END IF

			DO i=1, MaxDepth

				WRITE(2,5) list_depth(i), list_pom(i), list_enz(i), list_h(i), list_ba(i), list_hdom(i), list_bf(i),&
				 list_edom(i), list_denz(i), list_dedom(i), list_pomflux(i), list_baflux(i), list_enzflux(i), list_hflux(i),&
                list_denzflux(i)
				5 FORMAT(I4,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,&
1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5)
		
				WRITE(3,7) list_depth(i), list_attach(i), list_detach(i), list_ba_grow(i), list_bf_grow(i), list_b_ratio(i), list_bp_ratio(i),&
				list_up_ba(i), list_up_bf(i)
				7 FORMAT(I4,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5)
	
				WRITE(4,8) list_depth(i), list_pnum(i), list_pom_dgd(i), list_h_hdom(i), list_ehl_k(i), list_enz_dom(i), list_denz_dom(i),&
				list_bacover(i) 
				8 FORMAT(I4,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5,1X,1PE12.5)
				
		
			END DO
		
			STOP

	ENDIF


k = k + 1  !Depth Counter

!------------------
!  END DEPTH LOOP
!------------------

END DO

MaxDepth = k-1


!-------------------------------
!   PRINT PROGRESS TO TERMINAL
!-------------------------------

IF(CountTime .eq. 10000 .or. CountTime .eq. 20000 .or. CountTime .eq. 30000 .or. CountTime .eq. 40000 .or.&
 CountTime .eq. 50000 .or. CountTime .eq. 60000 .or. CountTime .eq. 70000 .or. CountTime .eq. 80000 .or.&
 CountTime .eq. 90000 .or. CountTime .eq. 100000 .or. CountTime .eq. 110000 .or. CountTime .eq. 120000 .or.&
 CountTime .eq. 130000 .or. CountTime .eq. 140000 .or. CountTime .eq. 150000 .or. CountTime .eq. 160000 .or.&
 CountTime .eq. 170000 .or. CountTime .eq. 180000 .or. CountTime .eq. 190000 .or. CountTime .eq. 200000 .or.&
 CountTime .eq. 210000 .or. CountTime .eq. 220000 .or. CountTime .eq. 230000 .or. CountTime .eq. 240000 .or.&
 CountTime .eq. 250000 .or. CountTime .eq. 260000 .or. CountTime .eq. 270000 .or. CountTime .eq. 280000 .or.&
 CountTime .eq. 290000 .or. CountTime .eq. 300000 .or. CountTime .eq. 310000 .or. CountTime .eq. 320000 .or.&
 CountTime .eq. 330000 .or. CountTime .eq. 340000 .or. CountTime .eq. 350000 .or. CountTime .eq. 360000 .or.&
 CountTime .eq. 370000 .or. CountTime .eq. 380000 .or. CountTime .eq. 390000 .or. CountTime .eq. 400000 .or.&
 CountTime .eq. 410000 .or. CountTime .eq. 420000 .or. CountTime .eq. 430000 .or. CountTime .eq. 440000 .or.&
 CountTime .eq. 450000 .or. CountTime .eq. 460000 .or. CountTime .eq. 470000 .or. CountTime .eq. 480000 .or.&
 CountTime .eq. 490000 .or. CountTime .eq. 500000 .or. CountTime .eq. 510000 .or. CountTime .eq. 520000)THEN
PRINT *, "Time Iteration Completed = ", CountTime
END IF



CLOSE(10)

CountTime = CountTime + 1 !Time Counter

!------------------
!  END TIME LOOP
!------------------

END DO

END PROGRAM MAIN