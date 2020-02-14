module global
use real_kind_mod
implicit none

! Units:
! From spheroid-abm, the max rates of consumption of oxygen and glucose are:
!   oxygen:  6.25e-17 moles/cell/s
!   glucose: 3.80e-17 moles/cell/s
! We work with mumol/cell/sec, and convert these values by scaling by 1.0e6, to give
!   oxygen:  6.25e-11 mumol/cell/s
!   glucose: 3.80e-11 mumol/cell/s
!
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM


integer, parameter :: ALIVE = 0
integer, parameter :: DEAD = 1
integer, parameter :: CFSE = 0
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: LACTATE = 3
integer, parameter :: GLUTAMINE = 4
integer, parameter :: OTHERNUTRIENT = 5
integer, parameter :: DRUG_A = 6
integer, parameter :: MAX_CELLTYPES = 2
integer, parameter :: MAX_DRUGTYPES = 2
integer, parameter :: MAX_CHEMO = 6	!?????????

integer, parameter :: nfin = 10, nfout = 11, nflog = 12

type metabolism_type
	real(REAL_KIND) :: HIF1
	real(REAL_KIND) :: PDK1
	real(REAL_KIND) :: I_rate_max
	real(REAL_KIND) :: G_rate
	real(REAL_KIND) :: PP_rate
	real(REAL_KIND) :: P_rate 
	real(REAL_KIND) :: L_rate
	real(REAL_KIND) :: A_rate
	real(REAL_KIND) :: I_rate
	real(REAL_KIND) :: O_rate
!	real(REAL_KIND) :: Itotal	! total of intermediates pool			NOT USED
!	real(REAL_KIND) :: I2Divide	! intermediates total needed to divide	NOT USED
	real(REAL_KIND) :: GA_rate
	real(REAL_KIND) :: f_G
	real(REAL_KIND) :: f_P
	real(REAL_KIND) :: C_P
	real(REAL_KIND) :: C_A
	! glutamine
	real(REAL_KIND) :: Gln_rate
	real(REAL_KIND) :: f_Gln
	real(REAL_KIND) :: recalcable   ! VBA doesn't like logical or integer
	! amino acids
	real(REAL_KIND) :: ON_rate
end type


type cell_type
	integer :: ID
	integer :: celltype
	integer :: site(3)
	integer :: ivin
	logical :: active
	integer :: state
	logical :: Iphase
!    integer :: nspheres             ! =1 for Iphase, =2 for Mphase
!	real(REAL_KIND) :: radius(2)	! sphere radii (um) -> cm
!	real(REAL_KIND) :: centre(3,2)  ! sphere centre positions
!	real(REAL_KIND) :: d			! centre separation distance (um) -> cm
	integer :: generation
!	real(REAL_KIND) :: conc(MAX_CHEMO)
	real(REAL_KIND) :: Cin(MAX_CHEMO)
!	real(REAL_KIND) :: Cex(MAX_CHEMO)
	real(REAL_KIND) :: dCdt(MAX_CHEMO)
	real(REAL_KIND) :: dMdt(MAX_CHEMO)      ! mumol/s
	real(REAL_KIND) :: CFSE
	real(REAL_KIND) :: dVdt
	real(REAL_KIND) :: V			! actual volume cm3
	real(REAL_KIND) :: divide_volume	! fractional divide volume (normalised)
	real(REAL_KIND) :: divide_time
	real(REAL_KIND) :: t_divide_last	! these two values are used for colony simulation
	real(REAL_KIND) :: t_divide_next
	real(REAL_KIND) :: birthtime
	real(REAL_KIND) :: t_anoxia
	real(REAL_KIND) :: t_anoxia_die
	real(REAL_KIND) :: t_aglucosia
	real(REAL_KIND) :: t_aglucosia_die
	real(REAL_KIND) :: M
	real(REAL_KIND) :: p_rad_death
	real(REAL_KIND) :: p_drug_death(MAX_DRUGTYPES)
	real(REAL_KIND) :: t_start_mitosis
	real(REAL_KIND) :: mitosis
	logical :: growth_delay
	real(REAL_KIND) :: dt_delay
	real(REAL_KIND) :: t_growth_delay_end			! this is for suppression of growth before first division
	integer :: N_delayed_cycles_left		! decremented by 1 at each cell division
	logical :: radiation_tag, anoxia_tag, aglucosia_tag
	logical :: drug_tag(MAX_DRUGTYPES)
	logical :: G2_M
!	logical :: exists
!	integer :: cnr(3,8)
!	real(REAL_KIND) :: wt(8)

	! Cell cycle 
    integer :: phase
    logical :: G1_flag, G1S_flag, G2_flag, G2M_flag
    real(REAL_KIND) :: G1_time, S_time, G2_time
    real(REAL_KIND) :: G1_V, S_V, G2_V
    real(REAL_KIND) :: G1S_time, G2M_time, M_time
    real(REAL_KIND) :: doubling_time
    integer :: NL1, NL2(2)
    logical :: starved
    
    ! Metabolism (HIF1, ATP_rate not needed here for monolayer)
!    real(REAL_KIND) :: HIF1, ATP_rate, I_rate
	type(metabolism_type) :: metab
    
	integer :: ndt

end type

type(cell_type), target :: cell_list(1), ccell_list(1)
integer :: Ncelltypes, nlist 
real(REAL_KIND) :: DELTA_T, Vcell_cm3, Vdivide0
logical :: colony_simulation, use_metabolism

! Metabolism parameters
real(REAL_KIND) :: f_Gu     ! unconstrained fraction of glycosis (r_G) going to make intermediates
real(REAL_KIND) :: f_Pu     ! unconstrained fraction of pyruvates (r_P) going to make intermediates
real(REAL_KIND) :: f_Glnu   ! unconstrained fraction of glutamine metabolism (r_Gln) going to make intermediates
real(REAL_KIND) :: N_GA		! number of ATP molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GP		! number of pyruvate molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GI		! number of intermediate molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_PA		! number of ATP molecules generated per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_PI		! number of intermediate molecules generated per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_PO		! number of O2 molecules consumed per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_ONI	! number of intermediate molecules generated per ON molecule
real(REAL_KIND) :: f_ATPg	! threshold ATP production rate fractions for cell growth
real(REAL_KIND) :: f_ATPs	! threshold ATP production rate fractions for cell survival
real(REAL_KIND) :: f_ATPramp	! multiplying factor for ramp start for reducing r_G, r_P
real(REAL_KIND) :: r_Ag		! threshold ATP production rates for cell growth
real(REAL_KIND) :: r_As		! threshold ATP production rates for cell survival
real(REAL_KIND) :: CO_H		! threshold O2 for Ofactor
real(REAL_KIND) :: CG_H		! threshold glucose for Gfactor
real(REAL_KIND) :: N_GlnA, N_GlnI, N_GlnO
real(REAL_KIND) :: Km_citrate 
! By cell
type(metabolism_type), target :: metabolic

! From chemokine.f90
type chemokine_type
	character(24) :: name
	logical :: used
	logical :: present
	logical :: use_secretion
	logical :: constant
	logical :: controls_growth
	logical :: controls_death
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: membrane_diff_in
	real(REAL_KIND) :: membrane_diff_out
	logical :: decay
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: decay_rate
	real(REAL_KIND) :: max_cell_rate	! Vmax
	real(REAL_KIND) :: MM_C0			! Km
	real(REAL_KIND) :: Hill_N           ! N
	real(REAL_KIND) :: medium_diff_coef	! diffusion coefficient in the medium
	real(REAL_KIND) :: medium_dlayer	! unstirred layer thickness
	real(REAL_KIND) :: medium_M			! mass of constituent
	real(REAL_KIND) :: medium_U			! total blob uptake rate
	real(REAL_KIND) :: medium_Cext		! far-field concentration
	real(REAL_KIND) :: medium_Cbnd		! boundary concentration
	real(REAL_KIND) :: total_flux_prev
	real(REAL_KIND) :: medium_Cbnd_prev

	real(REAL_KIND) :: diff_reduction_factor
	real(REAL_KIND), allocatable :: conc
end type

type(chemokine_type), target :: chemo(MAX_CHEMO)
real(REAL_KIND) :: Caverage(2*MAX_CHEMO)    ! average cell and average medium concentrations

integer, parameter :: N1D = 0
logical :: noSS = .false.

end module