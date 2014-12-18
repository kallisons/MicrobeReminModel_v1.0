PROGRAM MAIN
IMPLICIT NONE
!		"0D model of bacterial growth on particles that includes bacterial exoenzyme dynamics"
!
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
!
!   TIME
!
!       dt:				time step (h)
!
!       t:              time iteration
!
!
!   ENVIRONMENTAL PARAMETERS
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
!		dedom:			inactive exoenzyme in dissolved environment (mg C m^-3)
!
!
!   TEMPERATURE DEPENDENCE
!
!       q10_b:			q10 value for bacterial temperature dependence.
!							
!       q10_enz:		q10 value for exoenzyme activity rate.
!		
!       q10_hl:			q10 value for exoenzyme half life.
!
!
!   PARTICULATE ORGANIC MATTER PARAMETERS
!
!       r_top:			radius of the particle (m)
!
!		w:              flux rate of pom in m h^-1 (In 0D, used to calculate advective flux from particles)
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
!		bge_max:		maximum bacterial growth efficiency (Fig. 3, Del Giorgio and Cole 2000)
!
!		basal:			basal energetic cost at 20°C (h^-1) (based on the 0.01 day^-1 at 0°C from COBALT v1, Stock et al. 2013)
!		
!		quorum:			bacterial biomass at quorum sensing threshold (mg C m^-3 particle)
!
!       qfrac:          fraction of the surface area of the particle covered by particle-attached bacteria
!
!       max_epsilon:    maximum exoenzyme production rate (h^-1)
!
!		m_bf:           mortality free-living bacteria (h^-1)
!
!		m_ba:           mortality of attached bacteria (h^-1)
!
!       cperb:          mg carbon per bacteria (mg C (bacteria)^-1)
! 
!       bdiam:          diameter of a rod-shaped bacteria (m)
!
!       blength:        length of a rod-shaped bacteria (m)
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
!	DIFFUSION PARAMETERS
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
!------------------------------------
!   FORTRAN VARIABLE INITIALIZATION
!------------------------------------

	REAL::	hdom_initial
	REAL::	bf_initial
	REAL::	enz_initial
	REAL::	denz_initial
	REAL::	h_initial
	REAL::	edom_initial
	REAL::	dedom_initial
	REAL::	q10_b
	REAL::	q10_enz
	REAL::	q10_hl
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
	REAL::	temp
	REAL::	s
	REAL::	rho
	INTEGER:: tt
	REAL::	sldoc
	INTEGER:: t
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
    REAL, DIMENSION(max_dim) :: list_up_ba = 0
    REAL, DIMENSION(max_dim) :: list_up_bf = 0
	REAL, DIMENSION(max_dim) :: list_epsilon = 0
	REAL, DIMENSION(max_dim) :: list_pnum = 0
	REAL, DIMENSION(max_dim) :: list_pom_dgd = 0
	REAL, DIMENSION(max_dim) :: list_h_hdom = 0
	REAL, DIMENSION(max_dim) :: list_ehl_k = 0
	REAL, DIMENSION(max_dim) :: list_enz_dom = 0
	REAL, DIMENSION(max_dim) :: list_denz_dom = 0
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
	REAL::  Vtotal
	REAL::	Re
	REAL::	up_bf
	REAL::  hp
	REAL::  up_ba
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
	INTEGER:: i
    REAL:: hdom_check
	INTEGER:: finalt
	REAL::  qfrac
    REAL::  rho_Csub
    REAL::  Vcarbon
    INTEGER::  scenario
    REAL::  Ratio_mmC

!----------------------------------
!			NAMELIST
!----------------------------------

NAMELIST /control/ qfrac, max_epsilon


!----------------------------------
!		SCENARIO CHOICE
!----------------------------------
scenario = 2

!scenario = 1  (Interior)
!scenario = 2  (Interception)
!scenario = 3  (Retention)

PRINT*, "scenario = ", scenario


!----------------------------------
!	  DEFAULT PARAMETER VALUES
!----------------------------------

!******************************************
hdom_initial	= 1.5			!mg C m^-3 seawater
bf_initial		= 1.5			!mg C m^-3 seawater
enz_initial		= 0.			!mg C m^-3 seawater
denz_initial	= 0.			!mg C m^-3 seawater
h_initial		= 0.			!mg C m^-3 seawater
edom_initial	= 0.			!mg C m^-3 seawater
dedom_initial	= 0.			!mg C m^-3 seawater
!******************************************
q10_b			= 2.			!
q10_enz			= 2.			!
q10_hl			= 2.			!
!******************************************
dt				= 0.00833		!h  (30 s)
!******************************************
w				= 2.			!m h^-1	(~50 m/day)
r_top			= 0.001			!m (1 mm)
!******************************************
up_max_bf		= 0.25			!h^-1
up_max_ba		= 0.25			!h^-1
kh				= 200.          !5 mg C m^-3 particle
kdom			= 30.           !20 mg C m^-3 seawater
bge_max			= 0.4			!
basal			= 0.002			!h^-1 at 20°C
max_epsilon		= 0.028         !h^-1
m_bf			= 0.02			!h^-1
m_ba			= 0.04			!h^-1
!******************************************
ehl				= 30.			!30 at 20°C; 72 at 7°C (Steen and Arnosti 2011)
epom            = 118.          !at 20°C, q10=2 (Steen and Arnosti 2011)
rho_Csub        = 1.29E9        !mg m^-3 density of leucine
mr_enz			= 4.9E-9        !Leucine aminopeptidase	(Burley et al. 1990)
mr_h            = 1.5E-10       !m Leucine (Rawn 1989)
Ratio_mmC       = 131./72.      !dimensionless
!******************************************
f_CaCO3			= 0.075			!mg C (mg CaCO3)^-1	(Klaas and Archer 2002)
f_opal			= 0.029			!mg C (mg opal)^-1 (Klaas and Archer 2002)
f_lith			= 0.052			!mg C (mg lithogenic)^-1 (Klaas and Archer 2002)
flux_CaCO3      = 21.			!mg CaCO3 m^-2 d^-1 (Conte et al. 2001)
flux_opal       = 4.34			!mg opal m^-2 d^-1 (Conte et al. 2001)
flux_lith       = 6.77			!mg lith m^-2 d^-1 (Conte et al. 2001)
!******************************************
bc				= 1.38E-23		!m^2 kg s^-2 K^-1
!******************************************
cperb			= 5.0E-12		!mg C (bacteria)^-1
bdiam			= 5.0E-7		!m (Kirchman 2012)
blength			= 1.0E-6		!m (Kirchman 2012)
!******************************************
tt				= 5.			!years (Abell et al. 2000)
!******************************************
temp			= 20.           !degrees C, 19.19 BATS depth 150 m
s				= 36.632		!psu, BATS depth 150 m
rho				= 1026.226		!density, BATS depth 150 m
sldoc			= 174.77		!doc, BATS depth 150 m
!******************************************
!qfrac			= 1.0           !set in namelist
!max_epsilon     = 0.03         !set in namelist
!******************************************


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
    list_enz(1) = enz_initial
	list_denz(1) = denz_initial
	list_h(1) = h_initial
	list_edom(1) = edom_initial
	list_dedom(1) = dedom_initial
	list_time(1) = 0


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
	
	dv_sw = (dv_fw * (1. + E*SQRT(Clv) + F*Clv))*0.1    !kg m^-1 s^-1 (A7, Jumars et al. 1993)


!**Calculate kinematic viscosity**

    kv = (dv_sw/rho)*3600. !m^2 h^-1 (Jumars et al. 1993)


!**Decay of semi-labile to labile DOC per time step**

    decayK = (1./(tt*24.*365.))  !(Abell et al. 2000)


!----------------------------
!  PARTICLE CHARACTERISTICS
!----------------------------

!**Particle radius**

    r = r_top


!**Particle dry weight**

    dw = (17.*(r*2.*1000.)**1.8)  !(Iverson et al. 2010)
        !r*1000 is converting m to mm


!**Initial particulate organic matter (POM)**

    pnum = 84.
        !Number of particles in the shallowest box of the 1D model in steady state.
	
    list_pom(1) = pnum*(dw*0.001*0.12)  !(Iverson et al. 2010)
        !0.12 is the ratio of POC to dry weight (dw)
        !dw*0.001 is converting ug to mg


!**Particle volume**
	
	Vtotal = (4./3.)*3.14*r**3.


!**Reynolds Number of a particle**

	Re = (w*r)/kv
        !includes the effect of sinking speed


!**Quorum-sensing threshold**
	
	quorum = ((4.*3.14*r**2.)/(bdiam*blength))*cperb  ! mg C particle^-1
        !One layer of bacteria covering the surface area of the spherical particle


!**Initial particle-attached bacteria (PB)**

	list_ba(1) = quorum*pnum*qfrac


!------------------
!  START TIME LOOP
!------------------

finalt = 600

DO t = 2, finalt


!------------------------
!    RATE CONSTANTS
!------------------------	

!**Fraction of particle surface area covered by attached bacteria**

    IF(scenario .eq. 3)THEN

        bacover = (((list_ba(t-1)/pnum)/cperb)*(bdiam*blength))/(4.*3.14*r**2.)

        IF (bacover .gt. 1.)THEN

            bacover = 1.

        ENDIF

    ELSE

        bacover = 0.

    ENDIF


!**Particle degradation**
	
	pomp = list_pom(t-1)/(pnum*Vtotal)

	gamma = ((f_CaCO3*flux_CaCO3 + f_opal*flux_opal + f_lith*flux_lith)/24./w) !(Klaas and Archer 2002)

	gammap = gamma/(pnum*Vtotal)

	availC = (pomp-gammap)/pomp

	IF (availC .lt. 0.)THEN
	
        availC = 0.
	
	ENDIF
	
	pom_dgd = tfac_enz*epom*availC


!**Volume of substrate carbon available per particle**

    Vcarbon = (((list_pom(t-1)*availC)/pnum)*(Ratio_mmC))/rho_Csub


!**Bacterial attachment and detachment**

	attach = 0. !No attachment in the 0D model.

	detach = 0. !No detachment in the 0D model.


!**Mass transfer of hydrolysate from particles**

    diff_h = ((bc*(temp+273.16))/(6.*3.14*dv_sw*mr_h))*3600. !h^-1 (Equation 6, Jumars et al. 1993)
	
	sc_h = kv/diff_h !Schmidt Number
	
	sh_h = 1. + 0.619*(Re**0.412)*(sc_h)**(1./3.) !Sherwood Number
	
	h_hdom = (4.*3.14*diff_h*r*sh_h) !h^-1 (Equation 4, Kiorboe et al. 2001)


!**Mass transfer of exoenzymes (active and inactive) from particles**
	
	diff_enz = ((bc*(temp+273.16))/(6.*3.14*dv_sw*mr_enz))*3600. !h^-1 (Equation 6, Jumars et al. 1993)

    sc_enz = kv/diff_enz !Schmidt Number
	
	sh_enz = 1. + 0.619*(Re**0.412)*(sc_enz)**(1./3.) !Sherwood Number
	
	enz_edom = (4.*3.14*diff_enz*r*sh_enz) !h^-1 (Equation 4, Kiorboe et al. 2001)

	denz_dedom = (4.*3.14*diff_enz*r*sh_enz) !h^-1 (Equation 4, Kiorboe et al. 2001)


!**Bacterial uptake and growth**

    up_bf = tfac_b*up_max_bf*(list_hdom(t-1)/(kdom+list_hdom(t-1))) !free-living bacterial uptake

    hp = list_h(t-1)/(pnum*Vtotal)

    up_ba = tfac_b*up_max_ba*(hp/(kh+hp)) !particle-attached bacterial uptake

    epsilon = max_epsilon !exoenzyme production

    grow_bf = (bge_max)*up_bf - tfac_b*basal !free-living bacterial growth

    grow_ba = (bge_max)*up_ba - tfac_b*basal - epsilon !particle-attached bacterial growth


!**Exoenzyme half-life**

	ehl_k = (LOG(0.5)/(tfac_hl*ehl))*(-1.)


!------------------------
!        FLUXES
!------------------------

!**Sinking Rates In**

!No sinking fluxes in 0D
    fpom_in = 0.
	fba_in = 0.
	fenz_in = 0.
	fh_in = 0.
    fdenz_in = 0.


!**Sinking Rates Out* 

!No sinking fluxes in 0D
	fpom_out = 0.
	fenz_out = 0.
	fh_out = 0.
	fba_out = 0.
	fdenz_out = 0.


!------------------------
!        RATES
!------------------------

	jm_ba = m_ba*list_ba(t-1)**2.				!mortality of attached bacteria
	jpom_dgd = pom_dgd*list_enz(t-1)			!degradation of pom by exoenzyme
	jup_ba = up_ba*list_ba(t-1)					!uptake by attached bacteria

    IF(scenario .eq. 2)THEN

        hdom_check = h_hdom*((list_h(t-1)/(Vtotal*pnum)))*pnum - up_ba*list_ba(t-1)!hydrolysate interception

        IF (hdom_check .lt. 0.)THEN

            jh_hdom = 0. !hydrolysate flux rate out

        ELSE

            jh_hdom = hdom_check  !hydrolysate flux rate out

        ENDIF

    ELSE

        jh_hdom = (1.-bacover)*h_hdom*((list_h(t-1)/(Vtotal*pnum)))*pnum  !hydrolysate flux rate out

    ENDIF 

	jgrow_ba = grow_ba*list_ba(t-1)			!growth rate of attached bacteria	
	jattach = attach*list_bf(t-1)*pnum		!attachment rate of bacteria
	jdetach = detach*list_ba(t-1)			!detachment rate of bacteria
	jgrow_bf = grow_bf*list_bf(t-1)			!growth rate of free-living bacteria
	jm_bf = m_bf*list_bf(t-1)**2.			!mortality rate of free-living bacteria
	jepsilon = epsilon*list_ba(t-1)			!exoenzyme production rate by attached bacteria

    IF(scenario .eq. 3)THEN

        jenz_edom = (1.-bacover)*enz_edom*((list_enz(t-1)-list_enz(t-1)*(Vcarbon/Vtotal))/Vtotal) !active exoenzyme flux

    ELSE

        jenz_edom = enz_edom*((list_enz(t-1)-list_enz(t-1)*(Vcarbon/Vtotal))/Vtotal) !active exoenzyme flux

    ENDIF
    
    jm_enz = ehl_k*list_enz(t-1)				!particle exoenzyme deactivation rate
    jup_bf = up_bf*list_bf(t-1)                 !uptake by free-living bacteria
	jdedomdecay = decayK*list_dedom(t-1)		!dissolved inactive exoenzyme decay rate
	jdenzdecay = decayK*list_denz(t-1)          !particle inactive exoenzyme decay rate

    IF(scenario .eq. 3)THEN

        jdenz_dedom = (1.-bacover)*denz_dedom*((list_denz(t-1)/(Vtotal*pnum)))*pnum !inactive exoenzyme flux

    ELSE

        jdenz_dedom = denz_dedom*((list_denz(t-1)/(Vtotal*pnum)))*pnum !inactive exoenzyme flux

    ENDIF

	jm_edom = ehl_k*list_edom(t-1)	!dissolved exoenzyme deactivation rate


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

	list_time(t)    =   (t-1)*60
	list_pom(t)     =	list_pom(t-1) + d_pom
	list_enz(t)     =	list_enz(t-1) + d_enz
	list_h(t)       =	list_h(t-1) + d_h
	list_ba(t)      =	list_ba(t-1) + d_ba
	list_hdom(t)    =	list_hdom(t-1) + d_hdom
	list_bf(t)      =	list_bf(t-1) + d_bf
	list_denz(t)    =   list_denz(t-1) + d_denz
	list_edom(t)    =   list_edom(t-1) + d_edom
	list_dedom(t)   =   list_dedom(t-1) + d_dedom


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

    IF(list_h(t) .lt. 0.)THEN

        PRINT *, "***ERROR Negative Hydrolysate***"
        PRINT *, "CountTime =", t

    STOP
    END IF

!-----------------------------------
!   PRINT FINAL OUTPUT TO FILES
!-----------------------------------

	IF (t == finalt) THEN
			
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


!-----------------
! END TIME LOOP
!-----------------

END DO


END PROGRAM MAIN