! Program to test cellular metabolism 

module main_mod
use global
use metabolism

implicit none

type ab_type
	sequence
	real(REAL_KIND) :: a
	real(REAL_KIND) :: b
end type

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine SetupChemo
use global
use metabolism
integer :: ichemo

chemo(OXYGEN)%name = 'Oxygen'
chemo(GLUCOSE)%name = 'Glucose'
chemo(LACTATE)%name = 'Lactate'
chemo(GLUTAMINE)%name = 'Glutamine'
do ichemo = 1,4
	chemo(ichemo)%present = .true.
	chemo(ichemo)%used = .true.
	chemo(ichemo)%decay_rate = 0
	chemo(ichemo)%controls_growth = .false.
	chemo(ichemo)%controls_death = .false.
	chemo(ichemo)%max_cell_rate = chemo(ichemo)%max_cell_rate*1.0e6	! mol/cell/s -> mumol/cell/s
	chemo(ichemo)%MM_C0 = chemo(ichemo)%MM_C0/1000	! uM -> mM
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine setup_from_file(infile)
use global
use metabolism
character*(64) :: infile
integer :: nf

nf = nfin
open(nf,file=infile,status='old')
call ReadMetabolismParameters(nf)
read(nf,*)
read(nf,*) chemo(OXYGEN)%bdry_conc
read(nf,*) chemo(OXYGEN)%MM_C0
read(nf,*) chemo(OXYGEN)%max_cell_rate
read(nf,*) chemo(OXYGEN)%Hill_N
read(nf,*) chemo(GLUCOSE)%bdry_conc
read(nf,*) chemo(GLUCOSE)%MM_C0
read(nf,*) chemo(GLUCOSE)%max_cell_rate
read(nf,*) chemo(GLUCOSE)%Hill_N
read(nf,*) chemo(LACTATE)%bdry_conc
read(nf,*) chemo(LACTATE)%MM_C0
read(nf,*) chemo(LACTATE)%max_cell_rate
read(nf,*) chemo(LACTATE)%Hill_N
close(nfin)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine setup_from_arrays(nparam, nchem, param, chem)
use global
use metabolism
integer :: nparam, nchem
real(REAL_KIND) :: param(:), chem(:)

open(nflog,file='C:\temp\metab.log',status='replace')

f_Gu = param(1)
f_Pu = param(2)
f_Glnu = param(3)
N_GA = param(4)
N_GP = N_GA
N_PA = param(5)
N_GlnA = param(6)
N_GI = param(7)
N_PI = param(8)
N_GlnI = param(9)
N_PO = param(10)
N_GlnO = param(11)
K_H1 = param(12)
K_H2 = param(13)
K_HB = param(14)
K_PDK = param(15)
PDKmin = param(16)
C_O2_norm = param(17)
C_G_norm = param(18)
C_Gln_norm = param(19)
C_L_norm = param(20)
f_N = param(21)
f_ATPs = param(22)
f_ATPg = param(23)
f_ATPramp = param(24)
K_PL = param(25)
K_LP = param(26)
Hill_Km_P = param(27)
Hill_Km_C = param(28)
Gln_baserate = param(29)
Hill_N_P = 1
Hill_Km_P = Hill_Km_P/1000		! uM -> mM 

write(nflog,'(5e12.3)') param(1:nparam)
write(nflog,*)
chemo(OXYGEN)%bdry_conc = chem(1)
chemo(OXYGEN)%MM_C0 = chem(2)
chemo(OXYGEN)%max_cell_rate = chem(3)
chemo(OXYGEN)%Hill_N = int(chem(4))
chemo(GLUCOSE)%bdry_conc = chem(5)
chemo(GLUCOSE)%MM_C0 = chem(6)
chemo(GLUCOSE)%max_cell_rate = chem(7)
chemo(GLUCOSE)%Hill_N = int(chem(8))
chemo(LACTATE)%bdry_conc = chem(9)
chemo(LACTATE)%MM_C0 = chem(10)
chemo(LACTATE)%max_cell_rate = chem(11)
chemo(LACTATE)%Hill_N = int(chem(12))
! glutamine
chemo(GLUTAMINE)%bdry_conc = max(0.0001,chem(13))
chemo(GLUTAMINE)%MM_C0 = chem(14)
chemo(GLUTAMINE)%max_cell_rate = chem(15)
chemo(GLUTAMINE)%Hill_N = int(chem(16))
!C_Gn_test = chem(13)
!Hill_Km_Gln = chem(14)/1000		! uM -> mM
!Gn_maxrate = chem(15)*1.0e6
!Hill_N_Gln = int(chem(16))

write(nflog,'(5e12.3)') chem(1:nchem)
write(nflog,*)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine setup
use global
use metabolism
integer :: kcell
real(REAL_KIND) :: rsite(3)
type(metabolism_type), pointer :: mp
logical :: ok

nlist = 1
Ncelltypes = 1
DELTA_T = 600
Vcell_cm3 = 1.0e-9		! nominal cell volume in cm3
Vdivide0 = 1.6*Vcell_cm3
colony_simulation = .false.
use_metabolism = .true.

call SetupChemo

mp => metabolic
call SetupMetabolism(mp, ok)
!write(nflog,*) 'r_Ou,C_P: ',r_Ou,mp%C_P

kcell = 1
rsite = [0.,0.,0.]
call AddCell(kcell,rsite)
Caverage(OXYGEN) = chemo(OXYGEN)%bdry_conc
Caverage(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
Caverage(LACTATE) = chemo(LACTATE)%bdry_conc
Caverage(GLUTAMINE) = chemo(GLUTAMINE)%bdry_conc
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadMetabolismParameters(nf)
integer :: nf

read(nf,*) f_Gu
read(nf,*) f_Pu
read(nf,*) N_GA
read(nf,*) N_PA
read(nf,*) N_GI
read(nf,*) N_PI
read(nf,*) N_PO
read(nf,*) K_H1
read(nf,*) K_H2
read(nf,*) K_HB
read(nf,*) K_PDK
read(nf,*) PDKmin
read(nf,*) C_O2_norm
read(nf,*) C_G_norm
read(nf,*) C_L_norm
read(nf,*) f_ATPs
read(nf,*) f_ATPg
read(nf,*) f_ATPramp
read(nf,*) K_PL
read(nf,*) K_LP
read(nf,*) Hill_Km_P
Hill_N_P = 1
Hill_Km_P = Hill_Km_P/1000		! uM -> mM
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell, rsite)
use global
use metabolism
integer :: kcell
real(REAL_KIND) :: rsite(3)
integer :: ityp, k, kpar = 0
real(REAL_KIND) :: V0, Tdiv, Vdiv
type(cell_type), pointer :: cp
!type(cycle_parameters_type),pointer :: ccp
	
cp => cell_list(kcell)
!ccp => cc_parameters
cp%ID = kcell
cp%state = ALIVE
cp%generation = 1
!cp%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
cp%celltype = 1
ityp = cp%celltype
!Ncells_type(ityp) = Ncells_type(ityp) + 1
!cp%Iphase = .true.

V0 = Vdivide0/2
!cp%divide_volume = get_divide_volume(ityp, V0, Tdiv)
Vdiv = Vdivide0
Tdiv = 24*60*60
cp%divide_volume = Vdiv
cp%divide_time = Tdiv
cp%V = 0.75*Vdiv
cp%metab = metabolic

cp%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
cp%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
cp%Cin(LACTATE) = chemo(LACTATE)%bdry_conc
cp%Cin(GLUTAMINE) = chemo(GLUTAMINE)%bdry_conc
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine mp_execute(icase, ixl, dC, ng, isolver, iuse_gln, nparam, nchem, param, chem, CC, &
					r_An, mpArr)
implicit none
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE, ALIAS:"mp_execute" :: mp_execute
integer :: icase, ixl, ng, isolver, iuse_gln, nparam, nchem
real(REAL_KIND) :: dC, param(nparam), chem(nchem)
real(REAL_KIND) :: CC(ng), r_An
!real(REAL_KIND) :: Grate(ng), Arate(ng), Irate(ng), Prate(ng), Lrate(ng), Orate(ng)
type(metabolism_type) :: mpArr(ng)
integer :: ityp, i, no, k
real(REAL_KIND) :: C_G, C_O, C_L, C_Gln, C(4), MM_O2, r_G
character*(64) :: infile, outfile
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp
logical :: dbug

from_excel = (ixl == 1 )
!dbug = (abs(chem(1) - 0.01) < 0.0001)
use_glutamine = (iuse_gln == 1)
!outfile = 'PDK.out'
!open(nfout,file=outfile,status='replace')
!write(*,*) 'Enter metabolism parameter file name: '
!read(*,'(a)') infile
if (from_excel) then
	call setup_from_arrays(nparam, nchem, param, chem)
else
	infile = 'metab.inp'
	call setup_from_file(infile)
endif
call setup

cp => cell_list(1)
mp => cp%metab
ityp = cp%celltype
fgp_solver = isolver
write(nflog,*) 'use_glutamine: ',use_glutamine
!write(nflog,*) 'isolver,C_P,O_rate: ',isolver,mp%C_P,mp%O_rate
if (icase == 1) then	!Fix C_O, C_Gln vary C_G
	! set HIF1, PDK1 to steady-state levels
	mp%HIF1 = get_HIF1steadystate(cp%Cin(OXYGEN))
	call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
	MM_O2 = f_MM(cp%Cin(OXYGEN),Hill_Km_O2,int(Hill_N_O2))
	write(nflog,'(13a12)') 'C_G','Arate','Irate','Grate','Glnrate','Prate','Lrate','Orate','f_G','f_P','HIF1','fPDK','C_P'
	mp%O_rate = r_Ou
	do i = ng,1,-1
		C_O = cp%Cin(OXYGEN)
		C_G = i*dC
		C_L = cp%Cin(LACTATE)
		C_Gln = cp%Cin(GLUTAMINE)
		C = [C_O, C_G, C_L, C_Gln]
        mp%recalcable = -1
		if (use_nitrogen) then
    		call f_metab(mp, C_O, C_G, C_L, C_Gln)
		else
		    call set_f_gp(mp,C)
		    if (.not.solved) then
    		    call f_metab(mp, C_O, C_G, C_L, C_Gln)
    	    endif
    	endif
		write(nflog,'(13e12.3)') C_G,mp%A_rate,mp%I_rate,mp%G_rate,mp%Gln_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%f_G,mp%f_P,mp%HIF1,mp%PDK1,mp%C_P
		!write(nflog,'(a)') 'C_G,C_O,C_L,HIF1,PDK1,G_rate,A_rate,I_rate,P_rate,O_rate,L_rate'
		!write(nflog,'(6e12.3)') C_G,C_O,C_L,mp%HIF1,mp%PDK1,mp%G_rate,mp%A_rate,mp%I_rate,mp%P_rate,mp%O_rate,mp%L_rate
		CC(i) = C_G
		mpArr(i) = cp%metab
	enddo
elseif (icase == 2) then	!Fix C_O, C_Gln vary C_G
	! set HIF1, PDK1 to steady-state levels
	mp%HIF1 = get_HIF1steadystate(cp%Cin(OXYGEN))
	call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
!	write(nflog,'(a,2f9.3)') 'HIF1, PDK1: ',mp%HIF1,mp%PDK1
	MM_O2 = f_MM(cp%Cin(OXYGEN),Hill_Km_O2,int(Hill_N_O2))
!	write(nflog,'(a,2e12.3)') 'MM_O2,MM_O2*fPDK: ',MM_O2,MM_O2*mp%PDK1
	write(nflog,'(13a12)') 'C_G','Arate','Irate','Grate','Glnrate','Prate','Lrate','Orate','f_G','f_P','HIF1','fPDK','C_P'
	mp%O_rate = r_Ou
	do i = ng,1,-1
		C_O = cp%Cin(OXYGEN)
		C_G = cp%Cin(GLUCOSE)
		C_L = cp%Cin(LACTATE)
		C_Gln = i*dC
		C = [C_O, C_G, C_L, C_Gln]
        mp%recalcable = -1
		if (use_nitrogen) then
!	        r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,mp%O_rate)  ! Note: this is the previous O_rate
!		    call set_param_set(r_G,C_L)
    		call f_metab(mp, C_O, C_G, C_L, C_Gln)
		else
		    call set_f_gp(mp,C)
		    if (.not.solved) then
    		    call f_metab(mp, C_O, C_G, C_L, C_Gln)
    	    endif
    	endif
		write(nflog,'(13e12.3)') C_Gln,mp%A_rate,mp%I_rate,mp%G_rate,mp%Gln_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%f_G,mp%f_P,mp%HIF1,mp%PDK1,mp%C_P
!        write(nflog,*)
		!write(nflog,'(a)') 'C_G,C_O,C_L,HIF1,PDK1,G_rate,A_rate,I_rate,P_rate,O_rate,L_rate'
		!write(nflog,'(6e12.3)') C_G,C_O,C_L,mp%HIF1,mp%PDK1,mp%G_rate,mp%A_rate,mp%I_rate,mp%P_rate,mp%O_rate,mp%L_rate
		CC(i) = C_Gln
		mpArr(i) = cp%metab
	enddo
elseif (icase == 3) then	!Fix C_G, C_Gln vary C_O
	no = ng
	write(nflog,'(13a12)') 'C_O','Arate','Irate','Grate','Glnrate','Prate','Lrate','Orate','f_G','f_P','HIF1','fPDK','C_P'
	mp%O_rate = r_Ou
	do i = no,1,-1
		C_G = cp%Cin(GLUCOSE)
		C_O = i*dC
		C_L = cp%Cin(LACTATE)
		C_Gln = cp%Cin(GLUTAMINE)
		C = [C_O, C_G, C_L, C_Gln]
		mp%HIF1 = get_HIF1steadystate(C_O)
		call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
        mp%recalcable = -1
		if (use_nitrogen) then
    		call f_metab(mp, C_O, C_G, C_L, C_Gln)
		else
		    do k = 1,2
		    call set_f_gp(mp,C)
		    call f_metab(mp, C_O, C_G, C_L, C_Gln)
		    enddo
		endif
		write(nflog,'(13e12.3)') C_O,mp%A_rate,mp%I_rate,mp%G_rate,mp%Gln_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%f_G,mp%f_P,mp%HIF1,mp%PDK1,mp%C_P
		CC(i) = C_O
		mpArr(i) = cp%metab
	enddo
elseif (icase == 4) then	!Fix C_O, vary C_G and C_Gln
	no = ng
	write(nflog,'(13a12)') 'C_O','Arate','Irate','Grate','Glnrate','Prate','Lrate','Orate','f_G','f_P','HIF1','fPDK','C_P'
	do i = no,1,-1
		C_G = i*dC
		C_O = cp%Cin(OXYGEN)
		C_L = cp%Cin(LACTATE)
		C_Gln = i*dC
		C = [C_O, C_G, C_L, C_Gln]
		mp%HIF1 = get_HIF1steadystate(C_O)
		call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
        mp%recalcable = -1
		if (use_nitrogen) then
    		call f_metab(mp, C_O, C_G, C_L, C_Gln)
		else
		    call set_f_gp(mp,C)
		    call f_metab(mp, C_O, C_G, C_L, C_Gln)
		endif
		write(nflog,'(13e12.3)') C_G,mp%A_rate,mp%I_rate,mp%G_rate,mp%Gln_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%f_G,mp%f_P,mp%HIF1,mp%PDK1,mp%C_P
		CC(i) = C_G
		mpArr(i) = cp%metab
	enddo
endif
!write(nflog,*) 'Done'
close(nflog)
end subroutine



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine test_metab(n,mp)
implicit none
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE, ALIAS:"test_metab" :: test_metab
integer :: n
type(metabolism_type) :: mp(n)

mp(1) = cell_list(1)%metab
end subroutine

end module

#if 0
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
program main
use main_mod
implicit none
integer :: ixl, nparam, nchem
real(REAL_KIND) :: param(1), chem(1)

ixl = 0
nparam = 1
nchem = 1
call execute(ixl, nparam, nchem, param, chem)

end program
#endif

#if 0
dt = DELTA_T
dGdt = -0.00004    ! test rate of change of ambient glucose
!dOdt = -0.0000004    ! test rate of change of ambient O2
dOdt = -(cp%Cin(OXYGEN) - 0.00)/(200*dt)
!write(nfout,'(a)') '    hour        C_G        C_O      HIF-1       PDK1     G_rate     A_rate     I_rate     Itotal    P_rate     O_rate     L_rate'
write(nfout,'(a8,12a11)') 'hour','C_G','C_O','HIF-1','PDK1','G_rate','A_rate','I_rate','Itotal','P_rate','O_rate','L_rate'
do istep = 1,1	!2*6*24
	hr = istep*dt/3600
	if (vary_G) then
	    C_G = C_G + dGdt*dt
		C_G = max(C_G,0.0)
		C_G = min(C_G,5.5)
	endif
	if (vary_O) then
	    C_O = C_O + dOdt*dt
!	    C_O = max(C_O,0.005)
	    C_O = max(C_O,0.00)
		C_O = min(C_O,0.18)
	endif
	cp%Cin(OXYGEN) = C_O
	Caverage(OXYGEN) = C_O
	cp%Cin(GLUCOSE) = C_G
	Caverage(GLUCOSE) = C_G
	cp%Cin(LACTATE) = C_L
	Caverage(LACTATE) = C_L
	write(*,'(a,3f8.4)') 'C_O,C_G,C_L: ',C_O,C_G,C_L
!	call update_metabolism
!	mp%Itotal = mp%Itotal + dt*mp%I_rate
	HIF1 = metabolic%HIF1
!	call analyticSetHIF1(ityp,Caverage(OXYGEN),HIF1,dt)
	metabolic%HIF1 = HIF1
	PDK1 = metabolic%PDK1
!	call analyticSetPDK1(ityp,HIF1,PDK1,dt)
	metabolic%PDK1 = PDK1
	if (mod(istep,1) == 0) then
		write(nfout,'(f8.2,12e11.3)') hr,C_G,C_O,HIF1,PDK1,mp%G_rate,mp%A_rate,mp%I_rate,mp%P_rate,mp%O_rate,mp%L_rate
	endif
enddo
close(nfout)
#endif


