#!/usr/bin/perl 
#Revision date: 01/10/2015

############# Locations for Programs etc ################

# The directory containing the atomistic fragments and itp files
my $cg2atdir = "~/MemProtMD/CG2AT";

# Location for Gromacs
my $gromacs = "/usr/local/gromacs/bin";

# Location for Pymol (including executable) - requires edit for Macs
my $pymol = "/usr/local/pymol/1.4/pymol/pymol.exe";

# Location for Pulchra (including executable) - Only for 'full' conversion
my $pulchra = "/Users/seb/bin/pulchra";

# Location for Modeller (including executable) - Only for 'full' conversion
my $modeller = "/usr/local/bin/mod9.15";

#################### Options ############################

# CG Protein itp file 
# Used to centralise the system around the protein.
# Must include parameters for all proteins present.
my $itp = "protein-cg.itp";

# Number of individual proteins (multiples of protein-cg.itp)
# Do not change unless your system includes >=2 proteins.
my $protnum = 0;

# Lipid deletion method - alchembed, pymol, gmxsolvate or none - (required for "align" method) 
my $deletion = "alchembed";

# Align the original protein with just the TM domains or all residues (tm/all) 
my $align = "tm";

# Delete Waters? - pymol, gmxsolvate, killwater - order of effectiveness at deleting waters
my $delete_waters = "gmxsolvate";

# Optimise Box? Only works if Bilayer has formed along xy plane.
my $optimise = "yes";

# Merge Chains? - useful if inter-subunit disulphide is required (yes/no).
my $merge = "no";

# Have you installed modeller? (yes/no).
my $use_modeller = "no";

# Name of Gromacs (gmx, gmx_sse, gmx_avx, gmx_d etc)
my $gmx = "gmx";

# Solvate system?
my $solvate = "no";

###################################################################
############# Do not edit anything below this line ################
###################################################################

$ENV{GMXLIB} = "$gromacs/../share/gromacs/top";

unless ( $ARGV[0] )
{
	die "\nUsage: cg2at.pl CG-system.pdb [original.pdb] \n\nFurther Options within cg2at.pl\n\n";
}

my $infile = $ARGV[0];

my $method = "full";

if ($ARGV[1]){
$original = $ARGV[1];
$method = "align";
}

$ff = "GROMOS96-53a6";

#`rm -rf CG2AT`;
mkdir('CG2AT');
`cp $cg2atdir/*.itp CG2AT/`;
`cp $itp $infile CG2AT/`;
if ($ARGV[1]){
`grep -v OXT $original > CG2AT/cg2atalignfile.pdb`;
}
chdir "./CG2AT/";
mkdir('itp');


###################################################
#                                                 #
#               PREPARE CG SYSTEM                 #
#                                                 #
###################################################


system("$gromacs/$gmx editconf -resnr 1 -f $infile -o nolabel.gro");
system("$gromacs/$gmx editconf -f nolabel.gro -o nolabel.pdb");

$skip = 0;

$total=0;
$atom = 0;
open (INP, "nolabel.pdb") || die "Cannot open nolabel.pdb\n";
while (<INP>){
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh) = split;
	if ($woa eq ATOM){
	if ((($woc eq CA)||($woc =~ /^B/)||($woc =~ /^0B/)||($woc =~ /^4B/)||($woc =~ /^5B/))&&($wod ne DDM)&&($wod ne BOG)&&($wod ne SQD)&&($wod ne DGD)&&($wod ne LMG)){
		$xhbx[$woe] = $wof;
		$xhby[$woe] = $wog;
		$xhbz[$woe] = $woh;
		$total++;
	}
	if ($woc eq CA){
		$atom++;
		}
	}
}
close (INP);

if ($total > 0){
	if ($protnum==0){
		$protnum=1;
	}
open (ITP, "$itp") || ($skip = 1);
close (ITP);
}


if ($method eq "align"){
open (ORIG, "../$original") || die "\nCannot open $original\nCheck the name and location of your protein structure to be aligned\n\n";
close (ORIG);
}


# bond or martini
if ($atom > 0) {
$forcefield = "bond";
}else{
$forcefield = "martini";
}

$pip2 = 0;
$pip3 = 0;
$popc = 0;
$pope = 0;
$popg = 0;
$pops = 0;
$pvpg = 0;
$pvpe = 0;
$pvpa = 0;
$popa = 0;
$dopc = 0;
$dag = 0;
$mag = 0;
$dspc = 0;
$dspe = 0;
$dspg = 0;
$dppc = 0;
$dppe = 0;
$dppg = 0;
$ppcs = 0;
$dmpc = 0;
$dmpe = 0;
$dmpg = 0;
$dhpc = 0;
$bog = 0;
$ddm = 0;
$lmpg = 0;
$dpc = 0;
$card = 0;
$chol = 0;
$lmg = 0;
$sqd = 0;
$dgd = 0;
$sds = 0;
$sol = 0;
$na = 0;
$cl = 0;

open (INP, "nolabel.pdb") || die "Cannot open nolabel.pdb\n";
$atom = 0;
while (<INP>){
	chomp;
#ATOM  11856  CL- ION  7084      59.310  48.380  17.730  1.00  0.00
#ATOM   5544  C5B POPGZ 772     102.250  78.650  56.630  1.00  0.00
	$woc=substr($_,13,3);
	$wod=substr($_,17,3);
#	($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
	if (($woc eq PO1)&&($wod eq PIP)){
	$pip2++;
	}
	if (($woc eq PO0)&&($wod eq PI3)){
	$pip3++;
	}
	if (($woc eq NC3)&&($wod eq POP)){
	$popc++;
	}
	if (($woc eq NH3)&&($wod eq POP)){
	$pope++;
	}
	if ((($woc eq GLH)||($woc eq GL0))&&($wod eq POP)){
	$popg++;
	}	
	if (($woc eq NH3)&&($wod eq PVP)){
	$pvpe++;
	}
	if (($woc eq GLH)&&($wod eq PVP)){
	$pvpg++;
	}
	if (($woc eq PO3)&&($wod eq PVP)){
	$pvpa++;
	}
	if (($woc eq PO3)&&($wod eq POP)){
	$popa++;
	}
	if ((($woc eq NCO)||($woc eq CNO))&&($wod eq POP)){
	$pops++;
	}
	if (($woc eq NC3)&&($wod eq DOP)){
	$dopc++;
	}
	if (($woc eq NC3)&&($wod eq DPP)){
	$dppc++;
	}
	if (($woc eq NH3)&&($wod eq DPP)){
	$dppe++;
	}	
	if (($woc eq GLH)&&($wod eq DPP)){
	$dppg++;
	}
	if (($woc eq NC3)&&($wod eq PPC)){
	$ppcs++;
	}
	if (($woc eq NC3)&&($wod eq DLP)){
	$dmpc++;
	}
	if (($woc eq NH3)&&($wod eq DLP)){
	$dmpe++;
	}
	if (($woc eq GLH)&&($wod eq DLP)){
	$dmpg++;
	}
	if (($woc eq NC3)&&($wod eq DSP)){
	$dspc++;
	}
	if (($woc eq NH3)&&($wod eq DSP)){
	$dspe++;
	}
	if (($woc eq GLH)&&($wod eq DSP)){
	$dspg++;
	}
	if (($woc eq NC3)&&($wod eq DHP)){
	$dhpc++;
	}
	if (($woc eq B1)&&($wod eq DDM)){
	$ddm++;
	}
	if (($woc eq PO4)&&($wod eq LMP)){
	$lmpg++;
	}
	if (($woc eq B1)&&($wod eq BOG)){
	$bog++;
	}
	if (($woc eq NC3)&&($wod eq DPC)){
	$dpc++;
	}
	if (($woc eq ROH)&&($wod eq CHO)){
	$chol++;
	}		
	if (($woc eq PO1)&&($wod eq CAR)){
	$card++;
	}	
	if (($woc eq B0)&&($wod eq LMG)){
	$lmg++;
	}	
	if (($woc eq S1)&&($wod eq SQD)){
	$sqd++;
	}	
	if (($woc eq B6)&&($wod eq DGD)){
	$dgd++;
	}
	if (($woc eq DOH)&&($wod eq DAG)){
	$dag++;
	}	
	if (($woc eq DOH)&&($wod eq MAG)){
	$mag++;
	}
	if (($woc eq GL1)&&($wod eq SDS)){
	$sds++;
	}
	if ($woc eq W){
	$sol++;
	}
	if ($woc =~ m/^NA/){
	$na++;
	}
	if ($woc =~ m/^CL/){
	$cl++;
	}
}
close (INP);

$lipidnum=($popc+$pope+$popg+$pops+$pvpa+$popa+$dppc+$dppe+$dppg+$ppcs+$card+$pvpe+$pvpg+$dopc+$dspc+$dspe+$dspg+$dmpc+$dmpe+$dmpg+$pip2+$pip3+$dag+$mag+$chol+$dhpc+$dpc+$ddm+$lmpg+$bog+$lmg+$dgd+$sqd+$sds);

open (CGT, "> cg.top") || die "Cannot open cg.top\n";

if ($forcefield eq "bond") {
	print CGT "\n",'#include "ff_v1.4_x.itp"';
	`rm martini_v2.2.itp`;
}else {
	print CGT "\n",'#include "martini_v2.2.itp"';
	`rm ff_v1.4_x.itp`;
}

if ($protnum > 0) {
	print CGT "\n",'#include "',$itp,'"';
}

print CGT "\n",'[ system ]';
print CGT "\n",'[ molecules ]';

if ($protnum > 0) {
print CGT "\nProtein $protnum";
}

unless ($popc == 0){ 
print OUT "\nPOPC $popc";
}


unless ($popc == 0)
{print CGT "\nPOPC $popc";
}unless ($pope == 0){
print CGT "\nPOPE $pope";
}unless ($popg == 0){
print CGT "\nPOPG $popg";
}unless ($pops == 0){
print CGT "\nPOPS $pops";
}unless ($popa == 0){
print CGT "\nPOPA $popa";
}unless ($dppc == 0){
print CGT "\nDPPC $dppc";
}unless ($dppe == 0){
print CGT "\nDPPE $dppe";
}unless ($dppg == 0){
print CGT "\nDPPG $dppg";
}unless ($ppcs == 0){
print CGT "\nPPCS $ppcs";
}unless ($dopc == 0){
print CGT "\nDOPC $dopc";
}unless ($dspc == 0){
print CGT "\nDSPC $dspc";
}unless ($dspe == 0){
print CGT "\nDSPE $dspe";
}unless ($dspg == 0){
print CGT "\nDSPG $dspg";
}unless ($card == 0){
print CGT "\nCARD $card";
}unless ($pvpe == 0){
print CGT "\nPVPE $pvpe";
}unless ($pvpg == 0){
print CGT "\nPVPG $pvpg";
}unless ($pvpa == 0){
print CGT "\nPVPA $pvpa";
}unless ($dmpc == 0){
print CGT "\nDLPC $dmpc";
}unless ($dmpe == 0){
print CGT "\nDLPE $dmpe";
}unless ($dmpg == 0){
print CGT "\nDLPG $dmpg";
}unless ($pip2 == 0){
print CGT "\nPIP2 $pip2";
}unless ($pip3 == 0){
print CGT "\nPI3 $pip3";
}
if ($forcefield eq "martini"){
unless ($dag == 0){
print CGT "\nDAG $dag";
}unless ($mag == 0){
print CGT "\nMAG $mag";
}unless ($chol == 0){
print CGT "\nCHOL $chol";
}unless ($dhpc == 0){
print CGT "\nDHPC $dhpc";
}unless ($dpc == 0){
print CGT "\nDPC $dpc";
}unless ($ddm == 0){
print CGT "\nDDM $ddm";
}unless ($lmpg == 0){
print CGT "\nLMPG $lmpg";
}unless ($bog == 0){
print CGT "\nBOG $bog";
}unless ($lmg == 0){
print CGT "\nLMG $lmg";
}unless ($sqd == 0){
print CGT "\nSQD $sqd";
}unless ($dgd == 0){
print CGT "\nDGD $dgd";
}unless ($sds == 0){
print CGT "\nSDS $sds";
}
}
#print CGT "\nW $sol";
#print CGT "\nNA+ $na";
#print CGT "\nCL- $cl";
close (CGT);

system("mv *lipid*.itp itp/");

open (LIP, "> lipid.top") || die "No Output possible\n";
if ($ff eq "GROMOS96-43a1"){
print LIP "\n",'#include "gromos43a1.ff/forcefield.itp"';
print LIP "\n",'#include "itp/gromos-lipids.itp"';
print LIP "\n",'#include "itp/lipid-berger.itp"';
} elsif ($ff eq "GROMOS96-43a2"){
print LIP "\n",'#include "gromos43a2.ff/forcefield.itp"';
print LIP "\n",'#include "itp/gromos-lipids.itp"';
print LIP "\n",'#include "itp/lipid-berger.itp"';
} elsif ($ff eq "GROMOS96-45a3"){
print LIP "\n",'#include "gromos45a3.ff/forcefield.itp"';
print LIP "\n",'#include "itp/gromos-lipids.itp"';
print LIP "\n",'#include "itp/lipid-berger.itp"';
} elsif ($ff eq "GROMOS96-53a5"){
print LIP "\n",'#include "gromos53a5.ff/forcefield.itp"';
print LIP "\n",'#include "itp/gromos-lipids.itp"';
print LIP "\n",'#include "itp/lipid-berger.itp"';
} elsif ($ff eq "GROMOS96-53a6"){
print LIP "\n",'#include "gromos53a6.ff/forcefield.itp"';
print LIP "\n",'#include "itp/gromos-lipids.itp"';
print LIP "\n",'#include "itp/lipid-gmx53a6.itp"';
} elsif ($ff eq "OPLSUA"){
print LIP "\n",'#include "oplsaa.ff/forcefield.itp"';
print LIP "\n",'#include "itp/oplsua-lipids.itp"';
print LIP "\n",'#include "itp/lipid-oplsua.itp"';
} elsif ($ff eq "OPLSAA"){
print LIP "\n",'#include "oplsaa.ff/forcefield.itp"';
print LIP "\n",'#include "itp/oplsaa-lipids.itp"';
print LIP "\n",'#include "itp/lipid-oplsaa.itp"';
} elsif ($ff eq "CHARMM27"){
print LIP "\n",'#include "charmm27.ff/forcefield.itp"';
print LIP "\n",'#include "itp/charmm27-lipids.itp"';
} elsif ($ff eq "CHARMM36"){
print LIP "\n",'#include "charmm36.ff/forcefield.itp"';
print LIP "\n",'#include "itp/charmm36-lipids.itp"';
`cp -r $cg2atdir/charmm36.ff .`;
} elsif ($ff eq "AMBER"){
print LIP "\n",'#include "amber99.ff/forcefield.itp"';
print LIP "\n",'#include "itp/lipid-amber.itp"';
}
print LIP "\n",'[ system ]';
print LIP "\n",'[ molecules ]';
close (LIP);

if ($protnum > 0) {
	
system("$gromacs/$gmx make_ndx -f nolabel.pdb -o protein.ndx << EOD
del 0 
del 1-20
q
EOD");

system("$gromacs/$gmx editconf -resnr 1 -f nolabel.pdb -o prot.pdb -n protein.ndx");

open (PROT, "prot.pdb") || die "Cannot open file\n";
open (PRO, "> pro.pdb") || die "No Output possible\n";

while (<PROT>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woa eq ATOM){
	print PRO $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
}
close (PRO);
close (PROT);

}
open (INB, "nolabel.pdb") || die "No Output possible\n";
open (POPC, "> popc.pdb") || die "No Output possible\n";
open (POPE, "> pope.pdb") || die "No Output possible\n";
open (POPG, "> popg.pdb") || die "No Output possible\n";
open (PVPE, "> pvpe.pdb") || die "No Output possible\n";
open (PVPG, "> pvpg.pdb") || die "No Output possible\n";
open (POPS, "> pops.pdb") || die "No Output possible\n";
open (POPA, "> popa.pdb") || die "No Output possible\n";
open (PVPA, "> pvpa.pdb") || die "No Output possible\n";
open (PIP2, "> pip2.pdb") || die "No Output possible\n";
open (PIP3, "> pip3.pdb") || die "No Output possible\n";
open (DHPC, "> dhpc.pdb") || die "No Output possible\n";
open (DPPC, "> dppc.pdb") || die "No Output possible\n";
open (DPPE, "> dppe.pdb") || die "No Output possible\n";
open (DPPG, "> dppg.pdb") || die "No Output possible\n";
open (PPCS, "> ppcs.pdb") || die "No Output possible\n";
open (DMPC, "> dmpc.pdb") || die "No Output possible\n";
open (DMPE, "> dmpe.pdb") || die "No Output possible\n";
open (DMPG, "> dmpg.pdb") || die "No Output possible\n";
open (DOPC, "> dopc.pdb") || die "No Output possible\n";
open (DSPC, "> dspc.pdb") || die "No Output possible\n";
open (DSPE, "> dspe.pdb") || die "No Output possible\n";
open (DSPG, "> dspg.pdb") || die "No Output possible\n";
open (DDM, "> ddm.pdb") || die "No Output possible\n";
open (LMPG, "> lmpg.pdb") || die "No Output possible\n";
open (BOG, "> bog.pdb") || die "No Output possible\n";
open (DPC, "> dpc.pdb") || die "No Output possible\n";
open (CHOL, "> chol.pdb") || die "No Output possible\n";
open (CARD, "> card.pdb") || die "No Output possible\n";
open (LMG, "> lmg.pdb") || die "No Output possible\n";
open (SQD, "> sqd.pdb") || die "No Output possible\n";
open (DGD, "> dgd.pdb") || die "No Output possible\n";
open (DAG, "> dag.pdb") || die "No Output possible\n";
open (MAG, "> mag.pdb") || die "No Output possible\n";
open (SDS, "> sds.pdb") || die "No Output possible\n";
open (BOX, "> box.pdb") || die "No Output possible\n";

while (<INB>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woa eq CRYST1){
	print BOX $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	$xbox = $wob/10;
	$ybox = $woc/10;
	$zbox = $wod/10;
	}
	if (($wod eq PIP)||($wod eq PIP2)){
	print PIP2 $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if (($wod eq PI3)||($wod eq PIP3)){
	print PIP3 $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}	
	if (($wod eq POPC)||($wod eq POP)){
		if ($woc eq NC3){
			$popcid = $woe;
			print POPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popcid){
		print POPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq POPE)||($wod eq POP)){
		if ($woc eq NH3){
			$popeid = $woe;
			print POPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popeid){
		print POPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq POPG)||($wod eq POP)){
		if (($woc eq GLH)||($woc eq GL0)){
			$popgid = $woe;
			print POPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popgid){
		print POPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq POPS)||($wod eq POP)){
		if (($woc eq NCO)||($woc eq CNO)){
			$popsid = $woe;
			print POPS $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popsid){
		print POPS $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq POPA)||($wod eq POP)){
		if ($woc eq PO3){
			$popaid = $woe;
			print POPA $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popaid){
		print POPA $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	elsif (($wod eq DOPC)||($wod eq DOP)){
	print DOPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if (($wod eq DSPC)||($wod eq DSP)){
		if ($woc eq NC3){
			$dspcid = $woe;
			print DSPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dspcid){
		print DSPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DSPE)||($wod eq DSP)){
		if ($woc eq NH3){
			$dspeid = $woe;
			print DSPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dspeid){
		print DSPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DSPG)||($wod eq DSP)){
		if ($woc eq GLH){
			$dspgid = $woe;
			print DSPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dspgid){
		print DSPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DPPC)||($wod eq DPP)){
		if ($woc eq NC3){			
			$dppcid = $woe;
			print DPPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;			
			}
		elsif ($woe eq $dppcid){
		print DPPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}	
	if (($wod eq DPPE)||($wod eq DPP)){
		if ($woc eq NH3){			
			$dppeid = $woe;
			print DPPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;			
			}
		elsif ($woe eq $dppeid){
		print DPPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DPPG)||($wod eq DPP)){
		if ($woc eq GLH){			
			$dppgid = $woe;
			print DPPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;			
			}
		elsif ($woe eq $dppgid){
		print DPPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq PPCS)||($wod eq PPC)){
		if ($woc eq NC3){			
			$ppcsid = $woe;
			print PPCS $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;			
			}
		elsif ($woe eq $ppcsid){
		print PPCS $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DLPC)||($wod eq DLP)){
		if ($woc eq NC3){
			$dlpcid = $woe;
			print DMPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dlpcid){
		print DMPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DLPE)||($wod eq DLPE)){
		if ($woc eq NH3){
			$dlpeid = $woe;
			print DMPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dlpeid){
		print DMPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq DLPG)||($wod eq DLP)){
		if ($woc eq GLH){
			$dlpgid = $woe;
			print DMPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dlpgid){
		print DMPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq PVPE)||($wod eq PVP)){
		if ($woc eq NH3){
			$pvpeid = $woe;
			print PVPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $pvpeid){
		print PVPE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq PVPG)||($wod eq PVPG)){
		if ($woc eq GLH){
			$pvpgid = $woe;
			print PVPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $pvpgid){
		print PVPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	if (($wod eq PVPA)||($wod eq PVP)){
		if ($woc eq PO3){
			$popsid = $woe;
			print PVPA $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popsid){
		print PVPA $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		}
	}
	elsif (($wod eq DHPC)||($wod eq DHP)){
	print DHPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq DDM){
	print DDM $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif (($wod eq LMP)||($wod eq LMPG)){
	print LMPG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq BOG){
	print BOG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq DPC){
	print DPC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif (($wod eq CHO)||($wod eq CHOL)){
	print CHOL $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif (($wod eq CAR)||($wod eq CARD)){
	print CARD $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq LMG){
	print LMG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq SQD){
	print SQD $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq DGD){
	print DGD $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq DAG){
	print DAG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq MAG){
	print MAG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	elsif ($wod eq SDS){
	print SDS $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
}

close (INB);
close (POPC);
close (POPE);
close (POPG);
close (POPS);
close (PVPG);
close (PVPE);
close (POPA);
close (PVPA);
close (PIP2);
close (PIP3);
close (DHPC);
close (DPPC);
close (DPPE);
close (DPPG);
close (PPCS);
close (DMPC);
close (DMPE);
close (DMPG);
close (DOPC);
close (DSPC);
close (DSPE);
close (DSPG);
close (DDM);
close (LMPG);
close (BOG);
close (DPC);
close (CHOL);
close (CARD);
close (LMG);
close (SQD);
close (DGD);
close (DAG);
close (MAG);
close (SDS);
close (BOX);

`cat box.pdb pro.pdb popc.pdb pope.pdb popg.pdb pops.pdb popa.pdb dppc.pdb dppe.pdb dppg.pdb ppcs.pdb card.pdb pvpe.pdb pvpg.pdb pvpa.pdb dopc.pdb dspc.pdb dspe.pdb dspg.pdb dmpc.pdb dmpe.pdb dmpg.pdb pip2.pdb pip3.pdb dag.pdb mag.pdb chol.pdb dhpc.pdb dpc.pdb ddm.pdb lmpg.pdb bog.pdb lmg.pdb sqd.pdb dgd.pdb sds.pdb > order.pdb`;

open (EM, "> em.mdp") || die "Cannot open em.mdp\n";
print EM "cpp                      = /lib/cpp";
print EM "\nintegrator               = steep" ;
print EM "\nemstep                   = 0.001";
print EM "\nemtol                    = 1000";
print EM "\nepsilon_r                = 65" ;
print EM "\nnsteps                   = 5000";
close (EM);

open (EM, "> em-bfgs.mdp") || die "Cannot open em-bfgs.mdp\n";
print EM "cpp                      = /lib/cpp";
print EM "\nintegrator               = l-bfgs" ;
print EM "\nemstep                   = 0.001";
print EM "\nemtol                    = 1000";
print EM "\nepsilon_r                = 65" ;
print EM "\nnsteps                   = 5000";
close (EM);

open (EM, "> em-cg.mdp") || die "Cannot open em-cg.mdp\n";
print EM "cpp                      = /lib/cpp";
print EM "\nintegrator               = cg" ;
print EM "\nemstep                   = 0.001";
print EM "\nemtol                    = 1000";
print EM "\nepsilon_r                = 65" ;
print EM "\nnsteps                   = 5000";
close (EM);

open (EM, "> em-lipid1.mdp") || die "Cannot open em-lipid1.mdp\n";
print EM "cpp                      = /lib/cpp";
print EM "\ndefine               = -DPOSRES_LIPID" ;
print EM "\nintegrator               = cg" ;
print EM "\nemstep                   = 0.001";
print EM "\nemtol                    = 100";
print EM "\nepsilon_r                = 65" ;
print EM "\nnsteps                   = 5000";
close (EM);

open (EM, "> em-lipid2.mdp") || die "Cannot open em-lipid2.mdp\n";
print EM "cpp                      = /lib/cpp";
print EM "\nintegrator               = steep" ;
print EM "\nemstep                   = 0.001";
print EM "\nemtol                    = 10";
print EM "\nepsilon_r                = 65" ;
print EM "\nnsteps                   = 5000";
close (EM);

system("$gromacs/$gmx make_ndx -f order.pdb -o order.ndx  << EOD
del 2-30
1|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
!1|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
q
EOD");

if ($skip == 0){

system("$gromacs/$gmx editconf -resnr 1 -f order.pdb -c -o order-center1.pdb -n order.ndx -box $xbox $ybox 30 << EOD
1
0
EOD");

system("$gromacs/$gmx grompp -maxwarn 10 -f em.mdp -c order-center1.pdb -o em_inbox1.tpr -p cg.top");

open (ORD, "em_inbox1.tpr") || die "\n########################################################################\n\n Cannot set-up system for centralisation - Check your $itp file\n\n";

system("$gromacs/$gmx trjconv -f order-center1.pdb -s em_inbox1.tpr -center -pbc mol -o em_inbox1.pdb<< EOD
1
0
EOD");

system("$gromacs/$gmx editconf -resnr 1 -f em_inbox1.pdb -o em_labelZ.gro");
system("$gromacs/$gmx editconf -resnr 1 -f em_labelZ.gro -o em_nolabel.pdb");

system("$gromacs/$gmx traj -ox -f em_nolabel.pdb -s em_nolabel.pdb -n order.ndx -com <<EOD
2
EOD");

if ($optimise eq "yes"){

$com=0;

open(COM, "coord.xvg");
	while (<COM>) {
	chomp;
	($xoa, $xob, $xoc, $xod) = split;
	if (($xoa ne '#')||($xoa ne '@')){
	$com=$xod;
		}
	}
close (COM);

$top = ($com*10);
$bot = ($com*10);

open (TMP, "em_nolabel.pdb");
	while (<TMP>) {
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($woa eq "ATOM")&&(($woc eq "BB")||($woc eq "CA")||($woc =~ "^B")||($woc =~ /^4B/)||($woc =~ /^0B/)||($woc eq "C1A")||($woc eq "PO4")||($woc =~ /^5B/)||($woc eq "NC3"))){
			if ($woh > $top){
			$top = $woh;
			}
			elsif ($woh < $bot){
			$bot = $woh;
			}
		}
	}
	
$topdist = ((sqrt($top*$top))/10); 
$botdist = ((sqrt($bot*$bot))/10); 
$zaxis = (($topdist-$botdist)+4);

sub round {
	my($zaxis) = shift;
	return int($zaxis + 0.05 * ($zaxis <=> 0));
}

$bigz = round($zaxis);

system("$gromacs/$gmx editconf -resnr 1 -f em_inbox1.pdb -o order-center.pdb -box $xbox $ybox $bigz -c -n order.ndx << EOD
0
0
EOD");

system("$gromacs/$gmx grompp -maxwarn 10 -f em.mdp -c order-center.pdb -o em_inbox -p cg.top");
system("$gromacs/$gmx mdrun -deffnm em_inbox");
system("$gromacs/$gmx trjconv -f em_inbox.gro -s em_inbox.tpr -o em_inbox.pdb -center -pbc res <<EOD
1
0
EOD");

}
elsif ($optimise ne "yes"){
	system("$gromacs/$gmx editconf -resnr 1 -f em_nolabel.pdb -o em_inbox.pdb -box $xbox $ybox $zbox");
}
}
elsif ($skip == 1) {
system("$gromacs/$gmx editconf -resnr 1 -f order.pdb -o em_inbox.pdb -box $xbox $ybox $zbox");
}

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (BOX, "> box.pdb") || die "No Output possible\n";
while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woa eq CRYST1){
	print BOX $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
close (EMI);
close (BOX);

$total=0;
$atom = 0;
open (INP, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
while (<INP>){
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh) = split;
	if ($woa eq ATOM){
	if ((($woc eq CA)||($woc =~ /^B/)||($woc =~ /^0B/)||($woc =~ /^4B/)||($woc =~ /^5B/))&&($wod ne DDM)&&($wod ne BOG)){
		$xhbx[$woe] = $wof;
		$xhby[$woe] = $wog;
		$xhbz[$woe] = $woh;
		$total++;
	}
	if ($woc eq CA){
		$atom++;
		}
	}
}
close (INP);

@chain_num=();
$xi=0;
while ($xi<$total)
	{
	$xii=($xi+1);
	$distx=($xhbx[$xi]-$xhbx[$xii])**2;
	$disty=($xhby[$xi]-$xhby[$xii])**2;
	$distz=($xhbz[$xi]-$xhbz[$xii])**2;
	$dist=sqrt($distx+$disty+$distz);
	$dist=$dist/10;
	$dist = sprintf("%.3f", $dist);	
	if ($dist > 0.8){	
	push (@chain_num, $xii);
	}
	$xi++;
}

$chain = @chain_num;

if ($align eq "tm"){

system ("$gromacs/$gmx make_ndx -f em_inbox.pdb -o po4.ndx << EOD
del 0
del 1-40
aPO4
q
EOD");

system("$gromacs/$gmx traj -ox -f em_inbox.pdb  -s em_inbox.pdb  -n po4.ndx -com <<EOD
1
EOD");

$com=0;

open(COM, "coord.xvg");
	while (<COM>) {
	chomp;
	($xoa, $xob, $xoc, $xod) = split;
	if (($xoa ne '#')||($xoa ne '@')){
	$com=$xod;
		}
	}
close (COM);

open (TMP, "em_inbox.pdb");
open (OUT, "> tm.list");
print OUT "del 0\ndel 1-100\n";
my $resnum = -1;
my $resend = "no";
my $lastres = 0;
my $upper = ((($com)+2.0)*10);
my $lower = ((($com)-2.0)*10);

my $tm=0;
	while (<TMP>) {
	chomp;
	$woc=substr($_,12,4);
	$woe=substr($_,21,5);
	$woh=substr($_,46,8);
	if (($woc=~ m/^B/)||($woc=~ m/^0B/)||($woc=~ m/^4B/)||($woc=~ m/^5B/)||($woc=~ m/CA/)||($woc=~ m/BB/)){
		if (($woh < $upper)&&($woh > $lower)) {
		unless($resend eq "yes"){		
		print OUT "r",$woe,"-";
		$resend = "yes";
		$tm++;
			}
		if (($woe==$atom)&&($resend eq "yes")){
			print OUT $woe,"|";
			$resend = "no";		
			}
		}
		elsif((($woh >= $upper)||($woh <= $lower))&&($resend eq "yes")){
		$lastres=$woe-1;		
		print OUT $lastres,"|";
		$resend = "no";
		}
		$resnum++;
	}
}
if ($tm>0){
print OUT "\n0&1&aB*\n0&1&aCA\nq\nEOD\n";
}else{
print OUT "aB*&0\naCA&0\n1\n1\nq\nEOD\n";
} 
close(OUT);
close(TMP);
}



`cp $infile CG-initial-system.pdb`;
`cp em_inbox.pdb CG-centered.pdb`;

`cp box.pdb bilayer-box.pdb`;
`cp lipid.top popg.top`;
`cp lipid.top dppe.top`;
`cp lipid.top dppg.top`;
`cp lipid.top dspg.top`;
`cp lipid.top dmpg.top`;
`cp lipid.top pops.top`;
`cp lipid.top pvpg.top`;
`cp lipid.top ddm.top`;
`cp lipid.top dopc.top`;
`cp lipid.top lmg.top`;
`cp lipid.top dgd.top`;
`cp lipid.top sqd.top`;

unless (-e "atomistic-lipid.pdb"){

###################################################
#                                                 #
#              ADD LIPID FRAGMENTS                #
#                                                 #
###################################################

################## POPC GROMOS ####################

if (($popc > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPC, "> popc-cg.pdb") || die "No output possible\n";

$popcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPC)){
		if ($woc eq NC3){
			$popcid = $woe;	
			print POPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popcid){
		if ($woc eq PO4){	
			print POPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popc-cg.pdb -o popc-cg.pdb");

$num = 1;

while ($num <= $popc){

$loop = sprintf("%04d", $num);

open (ATM, "popc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPC $popc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popc-bilayer.pdb");
system("cat popc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPC/`;
}

############## POPC GROMOS END ######################
################ POPC GROMOS ########################

if (($popc > 0)&&($ff =~ m/^GROMOS/)&&($ff ne "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPC, "> popc-cg.pdb") || die "No output possible\n";

$popcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPC)){
		if ($woc eq NC3){
			$popcid = $woe;	
			print POPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popcid){
		if ($woc eq PO4){	
			print POPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popc-cg.pdb -o popc-cg.pdb");

$num = 1;

while ($num <= $popc){

$loop = sprintf("%04d", $num);

open (ATM, "popc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPC $popc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popc-bilayer.pdb");
system("cat popc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPC/`;
}

################## POPC GROMOS END ######################
#################### POPC OPLSUA ########################

if (($popc > 0)&&($ff eq "OPLSUA")) {

`cp -R $cg2atdir/POPC-oplsua .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPC, "> popc-cg.pdb") || die "No output possible\n";

$popcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPC)){
		if ($woc eq NC3){
			$popcid = $woe;	
			print POPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popcid){
		if ($woc eq PO4){	
			print POPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O35",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C32",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C40",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C44",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C48",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C51",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popc-cg.pdb -o popc-cg.pdb");

$num = 1;

while ($num <= $popc){

$loop = sprintf("%04d", $num);

open (ATM, "popc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPC $popc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-oplsua/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPC-oplsua/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'ATOM' | grep -v C26 | grep -v C29 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C34' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O35' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O35 | grep -v C40 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C43' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C44' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C44 | grep -v C48 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C51' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C52' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popc-bilayer.pdb");
system("cat popc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPC-oplsua/`;
}

################## POPC OPLSUA END ######################
#################### POPC CHARMM ########################

if (($popc > 0)&&($ff =~ m/^CHARMM/)) {

`cp -R $cg2atdir/POPC-CHARMM .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPC, "> popc-cg.pdb") || die "No output possible\n";

$popcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPC)){
		if ($woc eq NC3){
			$popcid = $woe;	
			print POPC $spa,$woa,$spb,$wob,$spc,"N  ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popcid){
		if ($woc eq PO4){	
			print POPC $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPC $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C28",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C2A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C2B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C2C",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C34",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C3A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPC $spa,$woa,$spb,$wob,$spc,"C3B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popc-cg.pdb -o popc-cg.pdb");

$num = 1;

while ($num <= $popc){

$loop = sprintf("%04d", $num);

open (ATM, "popc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPC $popc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPC-CHARMM/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'P1' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'O1' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O2' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPC-CHARMM/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep ATOM | grep -v 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O4 ' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep ATOM | grep 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'C3 ' | grep -v 'HX '| grep -v 'HY '| grep -v 'O31' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'O21'|grep -v 'C23' | grep -v 'H3R'| grep -v 'H3S'| grep -v 'C24'| grep -v 'H4R'| grep -v 'H4S' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'O2 ' | grep -v 'C1 '| grep -v 'HA '| grep -v 'HB ' | grep -v 'C2 '| grep -v 'HS '| grep -v 'O21'>> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C33' | grep -v 'O31'| grep -v 'H3X'| grep -v 'H3Y'| grep -v 'C34'| grep -v 'H4X'| grep -v 'H4Y' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'C21' | grep -v 'O21'| grep -v 'O22'| grep -v 'C22'| grep -v 'H2R'| grep -v 'H2S' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C2A' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' | grep -v 'C2B' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C31' | grep -v 'O31'| grep -v 'O32'| grep -v 'C32'| grep -v 'H2X'| grep -v 'H2Y' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C34' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'C38' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C3A' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popc-bilayer.pdb");
system("cat popc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPC-CHARMM/`;
}

############## POPC CHARMM END ###################
################ POPE GROMOS #####################

if (($pope > 0)&&($ff =~ m/^GROMOS/)&&($ff ne "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPE .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPE, "> pope-cg.pdb") || die "No output possible\n";

$popeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPE)){
		if ($woc eq NH3){
			$popeid = $woe;	
			print POPE $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popeid){
		if ($woc eq PO4){	
			print POPE $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPE $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPE $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pope-cg.pdb -o pope-cg.pdb");

$num = 1;

while ($num <= $pope){

$loop = sprintf("%04d", $num);

open (ATM, "pope-cg.pdb") || die "Cannot open pope-cg.pdb\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPE $pope";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pope){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPE/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NH3.pdb | grep 'H1' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'H2' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'H3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pope-bilayer.pdb");
system("cat pope-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPE/`;
}

################## POPE GROMOS END ######################
#################### POPE GROMOS ########################

if (($pope > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPE .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPE, "> pope-cg.pdb") || die "No output possible\n";

$popeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPE)){
		if ($woc eq NH3){
			$popeid = $woe;	
			print POPE $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popeid){
		if ($woc eq PO4){	
			print POPE $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPE $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPE $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pope-cg.pdb -o pope-cg.pdb");

$num = 1;

while ($num <= $pope){

$loop = sprintf("%04d", $num);

open (ATM, "pope-cg.pdb") || die "Cannot open pope-cg.pdb\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPE $pope";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pope){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPE/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NH3.pdb | grep 'H1' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'H2' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'H3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NH3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pope-bilayer.pdb");
system("cat pope-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPE/`;
}

################## POPE GROMOS END ######################
#################### POPE CHARMM ########################

if (($pope > 0)&&($ff =~ m/^CHARMM/)) {

`cp -R $cg2atdir/POPE-CHARMM .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPE, "> pope-cg.pdb") || die "No output possible\n";

$popeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPE)){
		if ($woc eq NH3){
			$popeid = $woe;	
			print POPE $spa,$woa,$spb,$wob,$spc,"N  ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popeid){
		if ($woc eq PO4){	
			print POPE $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPE $spa,$woa,$spb,$wob,$spc,"O21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPE $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C28",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C2A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C2B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C2C",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C34",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C3A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPE $spa,$woa,$spb,$wob,$spc,"C3B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pope-cg.pdb -o pope-cg.pdb");

$num = 1;

while ($num <= $pope){

$loop = sprintf("%04d", $num);

open (ATM, "pope-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPE $pope";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pope){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPE-CHARMM/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'P1' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'O1' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O2' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPE-CHARMM/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NH3.pdb | grep ATOM | grep -v 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O4 ' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep ATOM | grep 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'C3 ' | grep -v 'HX '| grep -v 'HY '| grep -v 'O31' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'O21'|grep -v 'C23' | grep -v 'H3R'| grep -v 'H3S'| grep -v 'C24'| grep -v 'H4R'| grep -v 'H4S' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'O2 ' | grep -v 'C1 '| grep -v 'HA '| grep -v 'HB ' | grep -v 'C2 '| grep -v 'HS '| grep -v 'O21'>> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C33' | grep -v 'O31'| grep -v 'H3X'| grep -v 'H3Y'| grep -v 'C34'| grep -v 'H4X'| grep -v 'H4Y' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'C21' | grep -v 'O21'| grep -v 'O22'| grep -v 'C22'| grep -v 'H2R'| grep -v 'H2S' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C2A' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' | grep -v 'C2B' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C31' | grep -v 'O31'| grep -v 'O32'| grep -v 'C32'| grep -v 'H2X'| grep -v 'H2Y' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C34' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'C38' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C3A' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > POPE-bilayer.pdb");
system("cat POPE-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPE-CHARMM/`;
}

############## POPE CHARMM END ###################
################ POPG GROMOS #####################

if (($popg > 0)&&($ff =~ m/^GROMOS/)&&($ff ne "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPG, "> popg-cg.pdb") || die "No output possible\n";

$popgid = 0;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPG)){
		if (($woc eq GLH)||($woc eq GL0)){
			$popgid = $woe;	
			print POPG $spa,$woa,$spb,$wob,$spc,"C2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popgid){
		if ($woc eq PO4){	
			print POPG $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPG $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popg-cg.pdb -o popg-cg.pdb");

$num = 1;

while ($num <= $popg){

$loop = sprintf("%04d", $num);

open (ATM, "popg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPG $popg";
close (LIP);

open (POPG, ">> popg.top") || die "No output possible\n";
print POPG "\nPOPG $popg";
close (POPG);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

`less fit_GLH.pdb | grep 'O7' > headgroup_chemistry.pdb`;
`less fit_GLH.pdb | grep 'P8' >> headgroup_chemistry.pdb`;
`less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb`;

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_GLH.pdb | grep 'H0' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O1' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C2 ' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O4' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'H5' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popg-bilayer.pdb");
system("cat box.pdb popg-bilayer.pdb > popg-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c popg-box.pdb -p popg.top -o popg_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm popg_em");
system("$gromacs/$gmx editconf -resnr 1 -f popg_em.gro -o popg_em.pdb");
system("cat popg_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPG/`;
}

############### POPG GROMOS END ##################
################# POPG GROMOS ####################

if (($popg > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPG, "> popg-cg.pdb") || die "No output possible\n";

$popgid = 0;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPG)){
		if (($woc eq GLH)||($woc eq GL0)){
			$popgid = $woe;	
			print POPG $spa,$woa,$spb,$wob,$spc,"C2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popgid){
		if ($woc eq PO4){	
			print POPG $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPG $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popg-cg.pdb -o popg-cg.pdb");

$num = 1;

while ($num <= $popg){

$loop = sprintf("%04d", $num);

open (ATM, "popg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPG $popg";
close (LIP);

open (POPG, ">> popg.top") || die "No output possible\n";
print POPG "\nPOPG $popg";
close (POPG);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

`less fit_GLH.pdb | grep 'O7' > headgroup_chemistry.pdb`;
`less fit_GLH.pdb | grep 'P8' >> headgroup_chemistry.pdb`;
`less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb`;

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_GLH.pdb | grep 'H0' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O1' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C2 ' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O4' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'H5' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popg-bilayer.pdb");
system("cat box.pdb popg-bilayer.pdb > popg-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c popg-box.pdb -p popg.top -o popg_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm popg_em");
system("$gromacs/$gmx editconf -resnr 1 -f popg_em.gro -o popg_em.pdb");
system("cat popg_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPG/`;
}

################## POPG GROMOS END ######################
#################### POPG CHARMM ########################

if (($popg > 0)&&($ff =~ m/^CHARMM/)) {

`cp -R $cg2atdir/POPG-CHARMM .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPG, "> popg-cg.pdb") || die "No output possible\n";

$popgid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPG)){
		if ($woc eq GLH){
			$popgid = $woe;	
			print POPG $spa,$woa,$spb,$wob,$spc,"C12",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popgid){
		if ($woc eq PO4){	
			print POPG $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPG $spa,$woa,$spb,$wob,$spc,"O21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPG $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C28",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C2A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C2B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C2C",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C34",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C3A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPG $spa,$woa,$spb,$wob,$spc,"C3B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popg-cg.pdb -o popg-cg.pdb");

$num = 1;

while ($num <= $popg){

$loop = sprintf("%04d", $num);

open (ATM, "popg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPG $popg";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPG-CHARMM/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_GLH.pdb | grep 'P1' > headgroup_chemistry.pdb");
system("less fit_GLH.pdb | grep 'O1' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O2' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPG-CHARMM/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_GLH.pdb | grep ATOM | grep -v 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O4 ' >> complete_${loop_main}.pdb");
system("less fit_GLH.pdb | grep ATOM | grep 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'C3 ' | grep -v 'HX '| grep -v 'HY '| grep -v 'O31' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'O21'|grep -v 'C23' | grep -v 'H3R'| grep -v 'H3S'| grep -v 'C24'| grep -v 'H4R'| grep -v 'H4S' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'O2 ' | grep -v 'C1 '| grep -v 'HA '| grep -v 'HB ' | grep -v 'C2 '| grep -v 'HS '| grep -v 'O21'>> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C33' | grep -v 'O31'| grep -v 'H3X'| grep -v 'H3Y'| grep -v 'C34'| grep -v 'H4X'| grep -v 'H4Y' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'C21' | grep -v 'O21'| grep -v 'O22'| grep -v 'C22'| grep -v 'H2R'| grep -v 'H2S' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C2A' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' | grep -v 'C2B' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C31' | grep -v 'O31'| grep -v 'O32'| grep -v 'C32'| grep -v 'H2X'| grep -v 'H2Y' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C34' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'C38' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C3A' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popg-bilayer.pdb");
system("cat popg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPG-CHARMM/`;
}

############## POPG CHARMM END ###################
############## POPS GROMOS OLD ###################

if (($pops > 0)&&($ff =~ m/^GROMOS/)&&($ff ne "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPS .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPS, "> pops-cg.pdb") || die "No output possible\n";

$popsid = 0;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPS)){
		if (($woc eq NCO)||($woc eq CNO)){
			$popsid = $woe;	
			print POPS $spa,$woa,$spb,$wob,$spc,"C4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popsid){
		if ($woc eq PO4){	
			print POPS $spa,$woa,$spb,$wob,$spc,"P11",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPS $spa,$woa,$spb,$wob,$spc,"O17",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPS $spa,$woa,$spb,$wob,$spc,"O38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C22",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C32",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C36",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C43",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C47",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C51",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C54",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPS);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pops-cg.pdb -o pops-cg.pdb");

$num = 1;

while ($num <= $pops){

$loop = sprintf("%04d", $num);

open (ATM, "pops-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPS $pops";
close (LIP);

open (POPS, ">> pops.top") || die "No output possible\n";
print POPS "\nPOPS $pops";
close (POPS);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pops){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/NCO.pdb -name -o fit_${loop_main}/fit_NCO.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

`less fit_NCO.pdb | grep 'O10' > headgroup_chemistry.pdb`;
`less fit_NCO.pdb | grep 'P11' >> headgroup_chemistry.pdb`;
`less fit_CCP.pdb | grep 'O14' >> headgroup_chemistry.pdb`;

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPS/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NCO.pdb | grep -v 'P11' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P11' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O12' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O13' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C15' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C16' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O17' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O17'| grep -v 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C24' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C25' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C26 | grep -v C29 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C34' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C35' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C36' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C37' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O38' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O38 | grep -v C43 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C43' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C44' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C45' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C47 | grep -v C51 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C51' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C52' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C53' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C54' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C55' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pops-bilayer.pdb");
system("cat box.pdb pops-bilayer.pdb > pops-box.pdb");
system("$gromacs/$gmx grompp -f em.mdp -c pops-box.pdb -p pops.top -o pops_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm pops_em");
system("$gromacs/$gmx editconf -resnr 1 -f pops_em.gro -o pops_em.pdb");
system("cat pops_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPS/`;
}

######### POPS GROMOS OLD END ############
############# POPS GROMOS ################

if (($pops > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPS-GRO .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPS, "> pops-cg.pdb") || die "No output possible\n";

$popsid = 0;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPS)){
		if (($woc eq NCO)||($woc eq CNO)){
			$popsid = $woe;	
			print POPS $spa,$woa,$spb,$wob,$spc,"C4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popsid){
		if ($woc eq PO4){	
			print POPS $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print POPS $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POPS $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POPS $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POPS $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(POPS);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pops-cg.pdb -o pops-cg.pdb");

$num = 1;

while ($num <= $pops){

$loop = sprintf("%04d", $num);

open (ATM, "pops-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPS $pops";
close (LIP);

open (POPS, ">> pops.top") || die "No output possible\n";
print POPS "\nPOPS $pops";
close (POPS);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pops){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/NCO.pdb -name -o fit_${loop_main}/fit_NCO.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-GRO/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

`less fit_NCO.pdb | grep 'O7' > headgroup_chemistry.pdb`;
`less fit_NCO.pdb | grep 'P8' >> headgroup_chemistry.pdb`;
`less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb`;

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPS-GRO/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NCO.pdb | grep -v 'P8' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pops-bilayer.pdb");
system("cat box.pdb pops-bilayer.pdb > pops-box.pdb");
system("$gromacs/$gmx grompp -f em.mdp -c pops-box.pdb -p pops.top -o pops_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm pops_em");
system("$gromacs/$gmx editconf -resnr 1 -f pops_em.gro -o pops_em.pdb");
system("cat pops_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPS/`;
}

################## POPS GROMOS END ######################
#################### POPS CHARMM ########################

if (($pops > 0)&&($ff =~ m/^CHARMM/)) {
`cp -R $cg2atdir/POPS-CHARMM .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POPS, "> pops-cg.pdb") || die "No output possible\n";

$popsid = 0 ;

while (<EMI>){
        chomp;
        ($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
        ($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
        if (($wod eq POP)||($wod eq POPS)){
                if (($woc eq NCO)||($woc eq CNO)){
                        $popsid = $woe;
                        print POPS $spa,$woa,$spb,$wob,$spc,"C12",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                elsif ($woe eq $popsid){
                if ($woc eq PO4){
                        print POPS $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq GL1){
                        print POPS $spa,$woa,$spb,$wob,$spc,"O21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq GL2){
                        print POPS $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C1B){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C2B){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C28",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq D3B){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C2A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C4B){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C2B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C5B){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C2C",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C1A){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C34",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C2A){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C3A){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C3A",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                if ($woc eq C4A){
                        print POPS $spa,$woa,$spb,$wob,$spc,"C3B",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                }
        }
}
close(EMI);
close(POPS);

#######`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pops-cg.pdb -o pops-cg.pdb");

$num = 1;

while ($num <= $pops){

$loop = sprintf("%04d", $num);

open (ATM, "pops-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
        chomp;
        ($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
        ($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
        if ($woe eq $num){
                if ((length($woe))==4){
                        print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
                        }
                elsif ((length($woe))==3){
                        print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
                        }
                elsif ((length($woe))==2){
                        print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
                        }
                elsif ((length($woe))==1){
                        print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
                        }
                }
        }
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPS $pops";
close (LIP);


open (POPS, ">> pops.top") || die "No output possible\n";
print POPS "\nPOPS $pops";
close (POPS);


# For every CG lipid
$num_main = 1;
while ($num_main <= $pops){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/NCO.pdb -name -o fit_${loop_main}/fit_NCO.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPS-CHARMM/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NCO.pdb | grep 'P1' > headgroup_chemistry.pdb");
system("less fit_NCO.pdb | grep 'O1' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O2' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../POPS-CHARMM/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NCO.pdb | grep ATOM | grep -v 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3 ' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O4 ' >> complete_${loop_main}.pdb");
system("less fit_NCO.pdb | grep ATOM | grep 'O1 ' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'C3 ' | grep -v 'HX '| grep -v 'HY '| grep -v 'O31' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'O21'|grep -v 'C23' | grep -v 'H3R'| grep -v 'H3S'| grep -v 'C24'| grep -v 'H4R'| grep -v 'H4S' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1' | grep 'ATOM'| grep -v 'O2 ' | grep -v 'C1 '| grep -v 'HA '| grep -v 'HB ' | grep -v 'C2 '| grep -v 'HS '| grep -v 'O21'>> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C33' | grep -v 'O31'| grep -v 'H3X'| grep -v 'H3Y'| grep -v 'C34'| grep -v 'H4X'| grep -v 'H4Y' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'C21' | grep -v 'O21'| grep -v 'O22'| grep -v 'C22'| grep -v 'H2R'| grep -v 'H2S' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C2A' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' | grep -v 'C2B' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C31' | grep -v 'O31'| grep -v 'O32'| grep -v 'C32'| grep -v 'H2X'| grep -v 'H2Y' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C34' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'C38' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C3A' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pops-bilayer.pdb");
system("cat box.pdb pops-bilayer.pdb > pops-box.pdb");
system("$gromacs/$gmx grompp -f em.mdp -c pops-box.pdb -p pops.top -o pops_em -maxwarn 5");
system("$gromacs/$gmx mdrun -deffnm pops_em");
system("$gromacs/$gmx editconf -resnr 1 -f pops_em.gro -o pops_em.pdb");
system("cat pops_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPS-CHARMM/`;
}

################### POPS CHARMM END #####################

#################### POPA GROMOS ########################

if (($popa > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/POPA .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (POP, "> popa-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq POP)||($wod eq POPA)){
		if ($woc eq PO3){
			$popaid = $woe;	
			print POP $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $popaid){
		if ($woc eq GL1){	
			print POP $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print POP $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print POP $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print POP $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print POP $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print POP $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print POP $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print POP $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print POP $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print POP $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print POP $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}
			}
		}
	}

close(EMI);
close(POP);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f popa-cg.pdb -o popa-cg.pdb");

$num = 1;

while ($num <= $popa){

$loop = sprintf("%04d", $num);

open (ATM, "popa-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPOPA $popa";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $popa){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 POPA/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

`less fit_CCP.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P1' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > popa-bilayer.pdb");
system("cat popa-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ POPA/`;
}

######################POPA_END###########################

##################### DPPC GROMOS #######################

if (($dppc > 0)&&(($ff =~ m/^GROMOS/)||($ff eq "OPLSUA"))) {

`cp -R $cg2atdir/DPPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPPC, "> dppc-cg.pdb") || die "No output possible\n";

$dppcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DPP)||($wod eq DPPC)){
		if ($woc eq NC3){
			$dppcid = $woe;	
			print DPPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dppcid){
		if ($woc eq PO4){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C27",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C30",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}
			}			
		}
	}

close(EMI);
close(DPPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dppc-cg.pdb -o dppc-cg.pdb");

$num = 1;

while ($num <= $dppc){

$loop = sprintf("%04d", $num);

open (ATM, "dppc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPPC $dppc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dppc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");


system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'N4' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'P8' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dppc-bilayer.pdb");
system("cat dppc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPPC/`;
}

################## DPPC GROMOS END ######################
#################### DPPC OPLSAA ########################

if (($dppc > 0)&&($ff =~ m/^OPLSAA/)) {

`cp -R $cg2atdir/DPPC-oplsaa .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPPC, "> dppc-cg.pdb") || die "No output possible\n";

$dppcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DPP)||($wod eq DPPC)){
		if ($woc eq NC3){
			$dppcid = $woe;	
			print DPPC $spa,$woa,$spb,$wob,$spc,"N1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}						
		elsif ($woe eq $dppcid){
		if ($woc eq PO4){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"O5 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"O7 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C28",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C32",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C36",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C40",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C11",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C15",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DPPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}
			}			
		}
	}

close(EMI);
close(DPPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dppc-cg.pdb -o dppc-cg.pdb");
Cannot Open
$num = 1;

while ($num <= $dppc){

$loop = sprintf("%04d", $num);

open (ATM, "dppc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPPC $dppc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dppc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");


system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPC-oplsaa/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O1' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P1' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O4' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPPC-oplsaa/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_C3A.pdb | grep 'C23' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep -v 'C23' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep -v 'C19' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v 'C15' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep -v 'C11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C7' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O4' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'P1' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep -v 'P1' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O2' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'P1'| grep -v 'O5'| grep -v 'C7' | grep -v 'C6'| grep -v 'O4'>> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep -v 'O7' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep -v 'C28' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep -v 'C32' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep -v 'C36' >> complete_${loop_main}.pdb");


chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dppc-bilayer.pdb");
system("cat dppc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPPC-oplsaa/`;
}

################## DPPC OPLSAA END ######################

##################### DPPE GROMOS #######################

if (($dppe > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DPPE .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPPE, "> dppe-cg.pdb") || die "No output possible\n";

$dppeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DPP)||($wod eq DPPE)){
		if ($woc eq NH3){
			$dppeid = $woe;	
			print DPPE $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dppeid){
		if ($woc eq PO4){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C27",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DPPE $spa,$woa,$spb,$wob,$spc,"C30",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}			
			}			
		}
	}

close(EMI);
close(DPPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dppe-cg.pdb -o dppe-cg.pdb");

$num = 1;

while ($num <= $dppe){

$loop = sprintf("%04d", $num);

open (ATM, "dppe-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPPE $dppe";
close (LIP);

open (DPPE, ">> dppe.top") || die "No output possible\n";
print DPPE "\nDPPE $dppe";
close (DPPE);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dppe){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");


system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPE/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPPE/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NH3.pdb | grep 'H1' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'H2' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'H3' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'N4' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'C5' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'P8' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dppe-bilayer.pdb");
system("cat box.pdb dppe-bilayer.pdb > dppe-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c dppe-box.pdb -p dppe.top -o dppe_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dppe_em");
system("$gromacs/$gmx editconf -resnr 1 -f dppe_em.gro -o dppe_em.pdb");
system("cat dppe_em.pdb | grep ATOM >> bilayer-box.pdb");

#system("cat fit_*/complete_*.pdb | grep 'ATOM' > dppe-bilayer.pdb");
#system("cat dppe-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPPE/`;
}

################## DPPE GROMOS END ######################

#################### DPPG GROMOS ########################

if (($dppg > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DPPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPPG, "> dppg-cg.pdb") || die "No output possible\n";

$dppgid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DPP)||($wod eq DPPG)){
		if ($woc eq GLH){
			$dppgid = $woe;	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dppgid){
		if ($woc eq PO4){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C27",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C30",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}
			}			
		}
	}

close(EMI);
close(DPPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dppg-cg.pdb -o dppg-cg.pdb");

$num = 1;

while ($num <= $dppg){

$loop = sprintf("%04d", $num);

open (ATM, "dppg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPPG $dppg";
close (LIP);

open (DPPG, ">> dppg.top") || die "No output possible\n";
print DPPG "\nDPPG $dppg";
close (DPPG);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dppg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");


system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_GLH.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_GLH.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_GLH.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dppg-bilayer.pdb");
system("cat box.pdb dppg-bilayer.pdb > dppg-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c dppg-box.pdb -p dppg.top -o dppg_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dppg_em");
system("$gromacs/$gmx editconf -resnr 1 -f dppg_em.gro -o dppg_em.pdb");
system("cat dppg_em.pdb | grep ATOM >> bilayer-box.pdb");

#system("cat fit_*/complete_*.pdb | grep 'ATOM' > dppg-bilayer.pdb");
#system("cat dppg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPPG/`;
}

################## DPPG GROMOS END #####################
#################### PPCS GROMOS #######################

if (($ppcs > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/PPCS .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (PPCS, "> ppcs-cg.pdb") || die "No output possible\n";

$ppcsid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq PPC)||($wod eq PPCS)){
		if ($woc eq NC3){
			$ppcsid = $woe;	
			print PPCS $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $ppcsid){
		if ($woc eq PO4){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq AM1){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C16",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq AM2){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C25",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C32",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D1B){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C39",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C43",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C47",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print PPCS $spa,$woa,$spb,$wob,$spc,"C50",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}			
	}
}

close(EMI);
close(PPCS);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f ppcs-cg.pdb -o ppcs-cg.pdb");

$num = 1;

while ($num <= $ppcs){

$loop = sprintf("%04d", $num);

open (ATM, "ppcs-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPPCS $ppcs";
close (LIP);

open (PPCS, ">> ppcs.top") || die "No output possible\n";
print PPCS "\nPPCS $ppcs";
close (PPCS);

# For every CG lipid
$num_main = 1;
while ($num_main <= $ppcs){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");


system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PPCS/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../PPCS/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'N14' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'H15' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C16' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O17' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C16' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C21'>> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C25' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' | grep -v 'C29' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C33' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O34' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'H35' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C33' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'C39' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C43' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'ATOM' | grep -v 'C47' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

#system("cat fit_*/complete_*.pdb | grep 'ATOM' > ppcs-bilayer.pdb");
#system("cat box.pdb ppcs-bilayer.pdb > ppcs-box.pdb");
#system("$gromacs/$gmx grompp -f em-cg.mdp -c ppcs-box.pdb -p ppcs.top -o ppcs_em -maxwarn 5"); 
#system("$gromacs/$gmx mdrun -nt 1 -deffnm ppcs_em");
#system("$gromacs/$gmx editconf -resnr 1 -f ppcs_em.gro -o ppcs_em.pdb");
#system("cat ppcs_em.pdb | grep ATOM >> bilayer-box.pdb");

system("cat fit_*/complete_*.pdb | grep 'ATOM' > ppcs-bilayer.pdb");
system("cat ppcs-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ PPCS/`;
}

################## PPCS GROMOS END #####################
#################### DPPG AMBER ########################

if (($dppg > 0)&&($ff eq "AMBER")) {

`cp -R $cg2atdir/DPPG-AMBER .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPPG, "> dppg-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DPP)||($wod eq DPPG)){
		if ($woc eq GLH){
			print DPPG $spa,$woa,$spb,$wob,$spc,"D2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"O7 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"O8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C13",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C17",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C25",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C37",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DPPG $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}

close(EMI);
close(DPPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dppg-cg.pdb -o dppg-cg.pdb");

$num = 1;

while ($num <= $dppg){

$loop = sprintf("%04d", $num);

open (ATM, "dppg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPPG $dppg";
close (LIP);

open (LIP, ">> dppg.top") || die "No output possible\n";
print LIP "\nDPPG 1";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dppg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPPG-AMBER/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_CCP.pdb | grep 'O1' > headgroup_chemistry.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P1' >> headgroup_chemistry.pdb");
system("less fit_GLH.pdb | grep 'O3' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPPG-AMBER/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_GL1.pdb | grep 'C12' >> C1A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C13' >> C1A-link.pdb`;
system("less fit_C1A.pdb | grep 'C14' >> C1A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1A-link.pdb -f2 ../DPPG-AMBER/C1AH.pdb -name -o fit_C1AH.pdb << EOD
1
1
EOD");

system("less fit_C1A.pdb | grep 'C16' >> C2A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C17' >> C2A-link.pdb`;
system("less fit_C2A.pdb | grep 'C18' >> C2A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2A-link.pdb -f2 ../DPPG-AMBER/C2AH.pdb -name -o fit_C2AH.pdb << EOD
1
1
EOD");

system("less fit_C2A.pdb | grep 'C20' >> C3A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C21' >> C3A-link.pdb`;
system("less fit_C3A.pdb | grep 'C22' >> C3A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3A-link.pdb -f2 ../DPPG-AMBER/C3AH.pdb -name -o fit_C3AH.pdb << EOD
1
1
EOD");

system("less fit_GL2.pdb | grep 'C28' >> C1B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> C1B-link.pdb`;
system("less fit_C1B.pdb | grep 'C30' >> C1B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1B-link.pdb -f2 ../DPPG-AMBER/C1BH.pdb -name -o fit_C1BH.pdb << EOD
1
1
EOD");

system("less fit_C1B.pdb | grep 'C32' >> C2B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> C2B-link.pdb`;
system("less fit_C2B.pdb | grep 'C34' >> C2B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2B-link.pdb -f2 ../DPPG-AMBER/C2BH.pdb -name -o fit_C2BH.pdb << EOD
1
1
EOD");

system("less fit_C2B.pdb | grep 'C36' >> C3B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> C3B-link.pdb`;
system("less fit_C3B.pdb | grep 'C38' >> C3B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3B-link.pdb -f2 ../DPPG-AMBER/C3BH.pdb -name -o fit_C3BH.pdb << EOD
1
1
EOD");


`cat ../box.pdb > complete_${loop_main}.pdb`;
system("less fit_CCP.pdb | grep ' C7' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' O1' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P1'| grep 'ATOM' >> complete_${loop_main}.pdb`;
system("less fit_GLH.pdb | grep -v 'D2'| grep -v 'P1'| grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep ' O4' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep ' O5' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H46' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C8' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
system("less fit_GL1.pdb | grep 'ATOM'| grep -v ' O7'| grep -v ' D7' | grep -v ' O9' | grep -v 'C13' |grep -v ' H' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
system("less fit_C1A.pdb | grep -v ' H' | grep 'ATOM'| grep -v 'C13' | grep -v 'C17' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C17' >> complete_${loop_main}.pdb`;
system("less fit_C2A.pdb | grep -v ' H' | grep 'ATOM'| grep -v 'C21' | grep -v 'C17' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
system("less fit_C3A.pdb | grep -v ' H' | grep 'ATOM' | grep -v 'C21' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3AH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C2AH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1AH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep ' O9' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H13' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C9' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H47' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H14' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O8' >> complete_${loop_main}.pdb`;
system("less fit_GL2.pdb | grep 'ATOM'| grep -v 'O8' | grep -v 'C29' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
system("less fit_C1BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep -v ' C33' |  grep -v ' C29'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
system("less fit_C2BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep -v ' C37'| grep -v ' C33' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> complete_${loop_main}.pdb`;
system("less fit_C3BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep -v ' C37' |grep 'ATOM' >> complete_${loop_main}.pdb");

system("$gromacs/$gmx grompp -f ../em-lipid1.mdp -c complete_${loop_main}.pdb -p ../dppg.top -o dppg_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dppg_${loop_main}");
system("$gromacs/$gmx grompp -f ../em-lipid2.mdp -c dppg_${loop_main}.gro -p ../dppg.top -o dppg_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dppg_${loop_main}");
system("$gromacs/$gmx editconf -resnr 1 -f dppg_${loop_main}.gro -o dppg_${loop_main}.pdb");

chdir "..";

$num_main++;
}


system("cat fit_*/dppg_*.pdb | grep 'ATOM' > dppg-bilayer.pdb");
system("cat dppg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPPG-AMBER/`;
}

#################### DPPG AMBER END ########################
################## CARDIOLIPIN GROMOS ######################

if (($card > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/CARD .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (CARD, "> card-cg.pdb") || die "No output possible\n";

$cardid = 0;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq CAR)||($wod eq CARD)){
		if ($woc eq GL0){
			$cardid = $woe;	
			print CARD $spa,$woa,$spb,$wob,$spc,"O4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $cardid){
		if ($woc eq PO1){	
			print CARD $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print CARD $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print CARD $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print CARD $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print CARD $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO2){	
			print CARD $spa,$woa,$spb,$wob," P8'  ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL3){	
			print CARD $spa,$woa,$spb,$wob," O14' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL4){	
			print CARD $spa,$woa,$spb,$wob," O33' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1D){	
			print CARD $spa,$woa,$spb,$wob," C19' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2D){	
			print CARD $spa,$woa,$spb,$wob," C23' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3D){	
			print CARD $spa,$woa,$spb,$wob," C26' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4D){	
			print CARD $spa,$woa,$spb,$wob," C29' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5D){	
			print CARD $spa,$woa,$spb,$wob," CA1' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1C){	
			print CARD $spa,$woa,$spb,$wob," C38' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2C){	
			print CARD $spa,$woa,$spb,$wob," C42' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3C){	
			print CARD $spa,$woa,$spb,$wob," C46' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4C){	
			print CARD $spa,$woa,$spb,$wob," C49' ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(CARD);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f card-cg.pdb -o card-cg.pdb");

$num = 1;

while ($num <= $card){

$loop = sprintf("%04d", $num);

open (ATM, "card-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nCARD $card";
close (LIP);

open (CARD, ">> card.top") || die "No output possible\n";
print CARD "\nCARD $card";
close (CARD);

# For every CG lipid
$num_main = 1;
while ($num_main <= $card){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/GLH1.pdb -name -o fit_${loop_main}/fit_GLH1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/CCP1.pdb -name -o fit_${loop_main}/fit_CCP1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/GLH2.pdb -name -o fit_${loop_main}/fit_GLH2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/CCP2.pdb -name -o fit_${loop_main}/fit_CCP2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/GL3.pdb -name -o fit_${loop_main}/fit_GL3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/GL4.pdb -name -o fit_${loop_main}/fit_GL4.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C1C.pdb -name -o fit_${loop_main}/fit_C1C.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C2C.pdb -name -o fit_${loop_main}/fit_C2C.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C3C.pdb -name -o fit_${loop_main}/fit_C3C.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C4C.pdb -name -o fit_${loop_main}/fit_C4C.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C1D.pdb -name -o fit_${loop_main}/fit_C1D.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C2D.pdb -name -o fit_${loop_main}/fit_C2D.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CARD/C3D.pdb -name -o fit_${loop_main}/fit_C3D.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

`less fit_GLH1.pdb | grep 'O7' > headgroup_chemistry1.pdb`;
`less fit_GLH1.pdb | grep 'P8' >> headgroup_chemistry1.pdb`;
`less fit_CCP1.pdb | grep 'O11' >> headgroup_chemistry1.pdb`;

`less fit_CCP2.pdb | grep 'O11' > headgroup_chemistry2.pdb`;
`less fit_GLH2.pdb | grep 'P8' >> headgroup_chemistry2.pdb`;
`less fit_GLH2.pdb | grep 'O7' >> headgroup_chemistry2.pdb`;

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry1.pdb -f2 ../CARD/PO1.pdb -name -o fit_PO1.pdb << EOD
1
1
EOD");
system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry2.pdb -f2 ../CARD/PO2.pdb -name -o fit_PO2.pdb << EOD
1
1
EOD");

`less fit_C3D.pdb | grep -v 'C46' >> complete_${loop_main}.pdb`;
`less fit_C2D.pdb | grep -v 'C42' >> complete_${loop_main}.pdb`;
`less fit_C1D.pdb | grep -v 'C38' >> complete_${loop_main}.pdb`;
`less fit_GL4.pdb | grep -v 'O33' >> complete_${loop_main}.pdb`;
`less fit_CCP2.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_CCP2.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less fit_C4C.pdb | grep -v 'C29' >> complete_${loop_main}.pdb`;
`less fit_C3C.pdb | grep -v 'C26' >> complete_${loop_main}.pdb`;
`less fit_C2C.pdb | grep -v 'C23' >> complete_${loop_main}.pdb`;
`less fit_C1C.pdb | grep -v 'C19' >> complete_${loop_main}.pdb`;
`less fit_GL3.pdb | grep -v 'O14' >> complete_${loop_main}.pdb`;
`less fit_CCP2.pdb | grep -v 'P8' | grep -v 'C32'| grep -v 'O33' >> complete_${loop_main}.pdb`;
`less fit_PO2.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_PO2.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_GLH2.pdb | grep -v 'O4' >> complete_${loop_main}.pdb`;
`less fit_GLH1.pdb  >> complete_${loop_main}.pdb`;
`less fit_PO1.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO1.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP1.pdb | grep -v 'P8' | grep -v 'C32'| grep -v 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep -v 'O14' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep -v 'C19' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep -v 'C23' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep -v 'C26' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep -v 'C29' >> complete_${loop_main}.pdb`;
`less fit_CCP1.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less fit_CCP1.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep -v 'O33' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep -v 'C38' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep -v 'C42' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep -v 'C46' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > card-bilayer.pdb");
system("cat box.pdb card-bilayer.pdb > card-box.pdb");
system("cat card-box.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ CARD/`;
}

################## CARD GROMOS END ######################
#################### PVPE GROMOS ########################

if (($pvpe > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/PVPE .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (PVPE, "> pvpe-cg.pdb") || die "No output possible\n";

$pvpeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq PVP)||($wod eq PVPE)){
		if ($woc eq NH3){
			$pvpeid = $woe;	
			print PVPE $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $pvpeid){
		if ($woc eq PO4){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print PVPE $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(PVPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pvpe-cg.pdb -o pvpe-cg.pdb");

$num = 1;

while ($num <= $pvpe){

$loop = sprintf("%04d", $num);

open (ATM, "pvpe-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPVPE $pvpe";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pvpe){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPE/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../PVPE/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NH3.pdb | grep 'H1' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'H2' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'H3' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'N4' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'C5' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'P8' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA1' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pvpe-bilayer.pdb");
system("cat pvpe-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ PVPE/`;
}

############### PVPE GROMOS END ######################
################# PVPG GROMOS ########################

if (($pvpg > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/PVPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (PVPG, "> pvpg-cg.pdb") || die "No output possible\n";

$pvpgid = 0;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq PVP)||($wod eq PVPG)){
		if ($woc eq GLH){
			$pvpgid = $woe;	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $pvpgid){
		if ($woc eq PO4){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print PVPG $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}
}
close(EMI);
close(PVPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pvpg-cg.pdb -o pvpg-cg.pdb");

$num = 1;

while ($num <= $pvpg){

$loop = sprintf("%04d", $num);

open (ATM, "pvpg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPVPG $pvpg";
close (LIP);

open (PVPG, ">> pvpg.top") || die "No output possible\n";
print PVPG "\nPVPG $pvpg";
close (PVPG);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pvpg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PVPG/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

`less fit_GLH.pdb | grep 'O7' > headgroup_chemistry.pdb`;
`less fit_GLH.pdb | grep 'P8' >> headgroup_chemistry.pdb`;
`less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb`;

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../PVPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_GLH.pdb | grep 'H0' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O1' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C2 ' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O4' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'H5' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pvpg-bilayer.pdb");
system("cat box.pdb pvpg-bilayer.pdb > pvpg-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c pvpg-box.pdb -p pvpg.top -o pvpg_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm pvpg_em");
system("$gromacs/$gmx editconf -resnr 1 -f pvpg_em.gro -o pvpg_em.pdb");
system("cat pvpg_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ PVPG/`;
}

################## PVPG GROMOS END ######################
#################### DOPC GROMOS ########################

if (($dopc > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DOPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DOPC, "> dopc-cg.pdb") || die "No output possible\n";

$dopcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DOP)||($wod eq POPC)){
		if ($woc eq NC3){
			$dopcid = $woe;	
			print DOPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dopcid){
		if ($woc eq PO4){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C45",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C48",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		if ($woc eq C5A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"CA3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}
			}
		}
	}

close(EMI);
close(DOPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dopc-cg.pdb -o dopc-cg.pdb");

$num = 1;

while ($num <= $dopc){

$loop = sprintf("%04d", $num);

open (ATM, "dopc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDOPC $dopc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dopc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DOPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C45 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C45' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4B.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA3' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA4' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dopc-bilayer.pdb");
system("cat dopc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DOPC/`;
}

################## DOPC GROMOS END ######################
#################### DOPC OPLSUA ########################

if (($dopc > 0)&&($ff eq "OPLSUA")) {

`cp -R $cg2atdir/DOPC-oplsua .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DOPC, "> dopc-cg.pdb") || die "No output possible\n";

$dopcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DOP)||($wod eq DOPC)){
		if ($woc eq NC3){
			print DOPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"O35",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C40",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C44",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C47",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C50",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C53",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		if ($woc eq C5A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C32",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}

close(EMI);
close(DOPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dopc-cg.pdb -o dopc-cg.pdb");

$num = 1;

while ($num <= $dopc){

$loop = sprintf("%04d", $num);

open (ATM, "dopc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDOPC $dopc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dopc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-oplsua/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DOPC-oplsua/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'N4 ' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb`;
`less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C34' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O35' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O35 | grep -v C40 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C43' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C44' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C44 | grep -v C47 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C50' >> complete_${loop_main}.pdb`;
`less fit_C4B.pdb | grep 'C51' >> complete_${loop_main}.pdb`;
`less fit_C4B.pdb | grep 'C52' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C53' >> complete_${loop_main}.pdb`;
`less fit_C4B.pdb | grep 'C54' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dopc-bilayer.pdb");
system("cat dopc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DOPC-oplsua/`;
}

################## DOPC OPLSUA END #####################
#################### DOPC AMBER ########################

if (($dopc > 0)&&($ff eq "AMBER")) {

`cp -R $cg2atdir/DOPC-AMBER .`;
`cp $cg2atdir/DOPC-AMBER/posre-dopc.itp ./itp/`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DOPC, "> dopc-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DOP)||($wod eq DOPC)){
		if ($woc eq NC3){
			print DOPC $spa,$woa,$spb,$wob,$spc,"DN1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"DP1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"O6 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"D7 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C30",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C34",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C37",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C40",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5A){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"C44",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"D4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"D8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"D11",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"D14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		if ($woc eq C5B){	
			print DOPC $spa,$woa,$spb,$wob,$spc,"D18",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}

close(EMI);
close(DOPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dopc-cg.pdb -o dopc-cg.pdb");

$num = 1;

while ($num <= $dopc){

$loop = sprintf("%04d", $num);

open (ATM, "dopc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDOPC $dopc";
close (LIP);

open (LIP, ">> dopc.top") || die "No output possible\n";
print LIP "\nDOPC 1";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dopc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DOPC-AMBER/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_CCP.pdb | grep 'O4' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P1' >> headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'O5' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DOPC-AMBER/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_C4A.pdb | grep 'C15' >> C4A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D14' >> C4A-link.pdb`;
system("less fit_C3A.pdb | grep 'C13' >> C4A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C4A-link.pdb -f2 ../DOPC-AMBER/C4AH.pdb -name -o fit_C4AH.pdb << EOD
1
1
EOD");

system("less fit_C3A.pdb | grep 'C12' >> C3A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D11' >> C3A-link.pdb`;
system("less fit_C2A.pdb | grep 'C10' >> C3A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3A-link.pdb -f2 ../DOPC-AMBER/C3AH.pdb -name -o fit_C3AH.pdb << EOD
1
1
EOD");

system("less fit_C2A.pdb | grep 'C9' >> C2A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D8' >> C2A-link.pdb`;
system("less fit_C1A.pdb | grep 'C7' >> C2A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2A-link.pdb -f2 ../DOPC-AMBER/C2AH.pdb -name -o fit_C2AH.pdb << EOD
1
1
EOD");

system("less fit_C1A.pdb | grep 'C5' >> C1A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D4' >> C1A-link.pdb`;
system("less fit_GL1.pdb | grep 'C3' >> C1A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1A-link.pdb -f2 ../DOPC-AMBER/C1AH.pdb -name -o fit_C1AH.pdb << EOD
1
1
EOD");

system("less fit_GL2.pdb | grep 'C29' >> C1B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C30' >> C1B-link.pdb`;
system("less fit_C1B.pdb | grep 'C31' >> C1B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1B-link.pdb -f2 ../DOPC-AMBER/C1BH.pdb -name -o fit_C1BH.pdb << EOD
1
1
EOD");

system("less fit_C1B.pdb | grep 'C33' >> C2B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C34' >> C2B-link.pdb`;
system("less fit_C2B.pdb | grep 'C35' >> C2B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2B-link.pdb -f2 ../DOPC-AMBER/C2BH.pdb -name -o fit_C2BH.pdb << EOD
1
1
EOD");

system("less fit_C2B.pdb | grep 'C36' >> C3B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> C3B-link.pdb`;
system("less fit_C3B.pdb | grep 'C38' >> C3B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3B-link.pdb -f2 ../DOPC-AMBER/C3BH.pdb -name -o fit_C3BH.pdb << EOD
1
1
EOD");

system("less fit_C3B.pdb | grep 'C39' >> C4B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C40' >> C4B-link.pdb`;
system("less fit_C4B.pdb | grep 'C41' >> C4B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C4B-link.pdb -f2 ../DOPC-AMBER/C4BH.pdb -name -o fit_C4BH.pdb << EOD
1
1
EOD");

`cat ../box.pdb > complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D18' >> complete_${loop_main}.pdb`;
system("less fit_C4A.pdb | grep -v ' D' | grep -v 'C14' | grep -v 'C18'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D14' >> complete_${loop_main}.pdb`;
system("less fit_C4AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep -v ' D' | grep -v 'C11' | grep -v 'C14'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D11' >> complete_${loop_main}.pdb`;
system("less fit_C3AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep -v ' D' | grep -v 'C8' | grep -v 'C11' |grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D8' >> complete_${loop_main}.pdb`;
system("less fit_C2AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v ' D' | grep -v 'C4' | grep -v 'C8' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D4' >> complete_${loop_main}.pdb`;
system("less fit_C1AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep -v ' D' | grep -v 'O7' | grep -v 'C4' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D7' >> complete_${loop_main}.pdb`;
system("less fit_CCP.pdb | grep -v ' D' | grep -v ' H4' | grep -v 'O6' | grep -v 'O7'| grep -v 'P1' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'DP1' >> complete_${loop_main}.pdb`;
system("less fit_PO4.pdb | grep 'O2' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep -v ' D' | grep -v 'P1' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H4' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O6' >> complete_${loop_main}.pdb`;
system("less fit_GL2.pdb | grep 'ATOM'| grep -v 'O6' | grep -v 'C30' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
system("less fit_C1BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep -v ' C30' |  grep -v ' C34'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C34' >> complete_${loop_main}.pdb`;
system("less fit_C2BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep -v ' C34'| grep -v ' C37' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> complete_${loop_main}.pdb`;
system("less fit_C3BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep -v ' C37' | grep -v ' C40'|grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
system("less fit_C4BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep -v ' C40' | grep 'ATOM' >> complete_${loop_main}.pdb");

system("$gromacs/$gmx grompp -f ../em-lipid1.mdp -c complete_${loop_main}.pdb -p ../dopc.top -o dopc_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dopc_${loop_main}");
system("$gromacs/$gmx grompp -f ../em-lipid2.mdp -c dopc_${loop_main}.gro -p ../dopc.top -o dopc_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dopc_${loop_main}");
system("$gromacs/$gmx editconf -resnr 1 -f dopc_${loop_main}.gro -o dopc_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/dopc_*.pdb | grep 'ATOM' > dopc-bilayer.pdb");
system("cat box.pdb dopc-bilayer.pdb > bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DOPC-AMBER/`;
}

################## DOPC AMBER END #######################
#################### DSPC GROMOS ########################

if (($dspc > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DSPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DSPC, "> dspc-cg.pdb") || die "No output possible\n";

$dspcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DSP)||($wod eq DSPC)){
		if ($woc eq NC3){
			$dspcid = $woe;	
			print DSPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dspcid){
		if ($woc eq PO4){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C45",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C48",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"CA3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		if ($woc eq C5A){	
			print DSPC $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}			
			}
		}
	}

close(EMI);
close(DSPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dspc-cg.pdb -o dspc-cg.pdb");

$num = 1;

while ($num <= $dspc){

$loop = sprintf("%04d", $num);

open (ATM, "dspc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDSPC $dspc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dspc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPC/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DSPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA1' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C46' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C50' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'CA3' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'CA4' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dspc-bilayer.pdb");
system("cat dspc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DSPC/`;
}

################### DSPC GROMOS END #####################
#################### DSPE GROMOS ########################

if (($dspe > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DSPE .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DSPE, "> dspe-cg.pdb") || die "No output possible\n";

$dspeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DSP)||($wod eq DSPE)){
		if ($woc eq NH3){
			$dspeid = $woe;	
			print DSPE $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dspeid){
		if ($woc eq PO4){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C45",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C48",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"CA3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		if ($woc eq C5A){	
			print DSPE $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}			
			}
		}
	}

close(EMI);
close(DSPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dspe-cg.pdb -o dspe-cg.pdb");

$num = 1;

while ($num <= $dspe){

$loop = sprintf("%04d", $num);

open (ATM, "dspe-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDSPE $dspe";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dspe){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPE/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DSPE/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NH3.pdb | grep 'H1' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'H2' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'H3' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'N4' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'C5' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NH3.pdb | grep 'P8' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA1' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C46' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C50' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'CA3' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'CA4' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dspe-bilayer.pdb");
system("cat dspe-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DSPE/`;
}

################### DSPE GROMOS END #####################
#################### DSPG GROMOS ########################

if (($dspg > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DSPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DSPG, "> dspg-cg.pdb") || die "No output possible\n";

$dspgid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DSP)||($wod eq DSPG)){
		if ($woc eq GLH){
			$dspgid = $woe;	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dspgid){
		if ($woc eq PO4){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C45",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C48",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"CA3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		if ($woc eq C5A){	
			print DSPG $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
				}			
			}
		}
	}

close(EMI);
close(DSPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dspg-cg.pdb -o dspg-cg.pdb");

$num = 1;

while ($num <= $dspg){

$loop = sprintf("%04d", $num);

open (ATM, "dspg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDSPG $dspg";
close (LIP);

open (DSPG, ">> dspg.top") || die "No output possible\n";
print DSPG "\nDSPG $dspg";
close (DSPG);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dspg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DSPG/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_GLH.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_GLH.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DSPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_GLH.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA1' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C46' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C48' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C49' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'C50' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'CA3' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'CA4' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dspg-bilayer.pdb");
system("cat box.pdb dspg-bilayer.pdb > dspg-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c dspg-box.pdb -p dspg.top -o dspg_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dspg_em");
system("$gromacs/$gmx editconf -resnr 1 -f dspg_em.gro -o dspg_em.pdb");
system("cat dspg_em.pdb | grep ATOM >> bilayer-box.pdb");

#system("cat fit_*/complete_*.pdb | grep 'ATOM' > dspg-bilayer.pdb");
#system("cat dspg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DSPG/`;
}

################### DSPG GROMOS END #####################
##################### DMPC GROMOS #######################

if ($dmpc > 0) {

`cp -R $cg2atdir/DMPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DMPC, "> dmpc-cg.pdb") || die "No output possible\n";

$dmpcid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DLP)||($wod eq DLPC)){
		if ($woc eq NC3){
			$dmpcid = $woe;	
			print DMPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dmpcid){
		if ($woc eq PO4){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"C36",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DMPC $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
}
close(EMI);
close(DMPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dmpc-cg.pdb -o dmpc-cg.pdb");

$num = 1;

while ($num <= $dmpc){

$loop = sprintf("%04d", $num);

open (ATM, "dmpc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDMPC $dmpc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dmpc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPC/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DMPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'N4' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'P8' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C23' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C37' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C38' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dmpc-bilayer.pdb");
system("cat dmpc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DMPC/`;
}

################## DMPC GROMOS END ####################
##################### DMPE GROMOS #######################

if ($dmpe > 0) {

`cp -R $cg2atdir/DMPE .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DMPE, "> dmpe-cg.pdb") || die "No output possible\n";

$dmpeid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DLP)||($wod eq DLPE)){
		if ($woc eq NH3){
			$dmpeid = $woe;	
			print DMPE $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dmpeid){
		if ($woc eq PO4){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"C36",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DMPE $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
}
close(EMI);
close(DMPE);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dmpe-cg.pdb -o dmpe-cg.pdb");

$num = 1;

while ($num <= $dmpe){

$loop = sprintf("%04d", $num);

open (ATM, "dmpe-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDMPE $dmpe";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dmpe){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/NH3.pdb -name -o fit_${loop_main}/fit_NH3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPE/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_NH3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NH3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DMPE/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NH3.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C23' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C37' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C38' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dmpe-bilayer.pdb");
system("cat dmpe-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DMPE/`;
}

################## DMPE GROMOS END ####################
#################### DMPG GROMOS ######################

if ($dmpg > 0) {

`cp -R $cg2atdir/DMPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DMPG, "> dmpg-cg.pdb") || die "No output possible\n";

$dmpgid = 0 ;

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DLP)||($wod eq DLPG)){
		if ($woc eq GLH){
			$dmpgid = $woe;	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ($woe eq $dmpgid){			
		if ($woc eq PO4){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"O31",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C36",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DMPG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
}
close(EMI);
close(DMPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dmpg-cg.pdb -o dmpg-cg.pdb");

$num = 1;

while ($num <= $dmpg){

$loop = sprintf("%04d", $num);

open (ATM, "dmpg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDMPG $dmpg";
close (LIP);

open (DMPG, ">> dmpg.top") || die "No output possible\n";
print DMPG "\nDMPG $dmpg";
close (DMPG);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dmpg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/GLH.pdb -name -o fit_${loop_main}/fit_GLH.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DMPG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_GLH.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_GLH.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DMPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

`less fit_GLH.pdb | grep 'H0' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O1' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C2' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C3' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O4' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'H5' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'C6' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'O7' >> complete_${loop_main}.pdb`;
`less fit_GLH.pdb | grep 'P8' >> complete_${loop_main}.pdb`;
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C23' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C30' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C37' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C38' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dmpg-bilayer.pdb");
system("cat box.pdb dmpg-bilayer.pdb > dmpg-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c dmpg-box.pdb -p dmpg.top -o dmpg_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dmpg_em");
system("$gromacs/$gmx editconf -resnr 1 -f dmpg_em.gro -o dmpg_em.pdb");
system("cat dmpg_em.pdb | grep ATOM >> bilayer-box.pdb");

#system("cat fit_*/complete_*.pdb | grep 'ATOM' > dmpg-bilayer.pdb");
#system("cat dmpg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DMPG/`;
}

################## DMPG GROMOS END #######################
#################### PIP2 GROMOS #########################

if (($pip2 > 0)&&($ff ne "GROMOS96-53a6")){

`cp -R $cg2atdir/PIP2 .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (PIP2, "> pip2-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod =~ /PIP/)||($wod =~ /PIP2/)){
		if ($woc eq PO1){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"PBN",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO2){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"PBH",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO3){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"PAW",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP1){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C5 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP2){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C3 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP3){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"OBT",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"OAR",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CBY",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CCC",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CCG",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CCK",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CCO",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CAM",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CAI",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CAE",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CAA",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
	}
}
close(EMI);
close(PIP2);

# Directory to store all the best-fitting AT lipids, for concatenation into AT bilayer structure
`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pip2-cg.pdb -o pip2-cg.pdb");

$num = 1;

while ($num <= $pip2){

$loop = sprintf("%04d", $num);

open (ATM, "pip2-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPIP2 $pip2";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pip2){

$loop_main = sprintf("%04d", $num_main);

mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

open (CCP, "fit_CCP.pdb") || die "Cannot open file\n";
open (RP, "fit_RP.pdb") || die "Cannot open file\n";
open (RP2, "fit_RP.pdb") || die "Cannot open file\n";
open (RP3, "fit_RP.pdb") || die "Cannot open file\n";
open (LOOP1, "../Res_CG/res_cg_${loop_main}.pdb") || die "Cannot open file\n";
open (LOOP2, "../Res_CG/res_cg_${loop_main}.pdb") || die "Cannot open file\n";
open (O1, "> O1.pdb") || die "No Output possible\n";
open (O2, "> O2.pdb") || die "No Output possible\n";
open (O3, "> O3.pdb") || die "No Output possible\n";

while (<CCP>){
	if ($_ =~ m/OAV/){ 
		print O1 $_ ;
		}
	if ($_ =~ m/PAW/){ 
		print O1 $_ ;
		}
	}

while (<RP>){
	if ($_ =~ m/OAZ/){ 
		print O1 $_ ;
		}
	if ($_ =~ m/OBM/){ 
		print O1 $_ ;
		}
	}


while (<LOOP1>){
	if ($_ =~ m/PBN/){ 
		print O2 $_ ;
		}
	}

while (<RP2>){
	if ($_ =~ m/OBM/){ 
		print O2 $_ ;
		}
	if ($_ =~ m/OAZ/){ 
		print O2 $_ ;
		}
	if ($_ =~ m/OBG/){ 
		print O3 $_ ;
		}
	}

while (<LOOP2>){
	if ($_ =~ m/PBH/){ 
		print O3 $_ ;
		}
	}

while (<RP3>){
	if ($_ =~ m/OAZ/){ 
		print O3 $_ ;
		}
	}

close (CCP);
close(RP);
close(RP2);
close(RP3);
close (LOOP1);
close (LOOP2);	
close (O1);
close (O2); 
close (O3);

system("$gromacs/$gmx confrms -nice 0 -one -f1 O1.pdb -f2 ../PIP2/PO1.pdb -name -o fit_PO1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O2.pdb -f2 ../PIP2/PO2.pdb -name -o fit_PO2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O3.pdb -f2 ../PIP2/PO3.pdb -name -o fit_PO3.pdb << EOD
1
1
EOD");

system("less fit_C4B.pdb | grep 'ATOM' | grep -v 'CCK' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'CCG'  >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'CCC'  >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'CBY'  >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'CBY'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'CBX'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'CBW'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'CBU'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'OBV'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'OBT'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'CAT'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'CAU'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'OAV'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'PAW'  >> complete_${loop_main}.pdb");
system("less fit_PO1.pdb | grep 'OAX'  >> complete_${loop_main}.pdb");
system("less fit_PO1.pdb | grep 'OAY'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OAZ'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C1 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C6 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBS'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'HBS'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C5 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBM'  >> complete_${loop_main}.pdb");
system("less fit_PO2.pdb | grep 'ATOM' | grep -v 'OBM'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C4 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBG'  >> complete_${loop_main}.pdb");
system("less fit_PO3.pdb | grep 'ATOM' | grep -v 'OBG'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C3 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBE'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'HBE'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C2 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBC'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'HBC'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'CAS'  >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'CAL' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'CAK' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'CAJ' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'CAE' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

`cat fit_*/complete_*.pdb | grep 'ATOM' > pip2-bilayer.pdb`;
`cat pip2-bilayer.pdb >> bilayer-box.pdb`;

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ PIP2/`;
}

################## PIP2 GROMOS END ######################
################## PIP2 GROMOS-CKP #######################

if (($pip2 > 0)&&($ff eq "GROMOS96-53a6"))  {

`cp -R $cg2atdir/PIP2-GRO .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (PIP2, "> pip2-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod =~ /PIP/)||($wod =~ /PIP2/)){
		if ($woc eq PO1){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"P4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO2){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"P5 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO3){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP1){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C5 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP2){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C3 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP3){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"O9 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"O28",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C18",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C24",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C37",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print PIP2 $spa,$woa,$spb,$wob,$spc,"C44",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
	}
}
close(EMI);
close(PIP2);

# Directory to store all the best-fitting AT lipids, for concatenation into AT bilayer structure
`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pip2-cg.pdb -o pip2-cg.pdb");

$num = 1;

while ($num <= $pip2){

$loop = sprintf("%04d", $num);

open (ATM, "pip2-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPIP2 $pip2";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pip2){

$loop_main = sprintf("%04d", $num_main);

mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP2-GRO/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");



chdir "fit_${loop_main}";

open (CCP, "fit_CCP.pdb") || die "Cannot open file\n";
open (RP, "fit_RP.pdb") || die "Cannot open file\n";
open (RP2, "fit_RP.pdb") || die "Cannot open file\n";
open (RP3, "fit_RP.pdb") || die "Cannot open file\n";
open (LOOP1, "../Res_CG/res_cg_${loop_main}.pdb") || die "Cannot open file\n";
open (LOOP2, "../Res_CG/res_cg_${loop_main}.pdb") || die "Cannot open file\n";
open (O1, "> O1.pdb") || die "No Output possible\n";
open (O2, "> O2.pdb") || die "No Output possible\n";
open (O3, "> O3.pdb") || die "No Output possible\n";

while (<RP>){
	if ($_ =~ m/1OP1/){ 
		print O1 $_ ;
		}
	if ($_ =~ m/1OP4/){ 
		print O2 $_ ;
		}
	}

while (<CCP>){
	if ($_ =~ m/4OP1/){ 
		print O1 $_ ;
		}
	if ($_ =~ m/P1/){ 
		print O1 $_ ;
		}
	}

while (<LOOP1>){
	if ($_ =~ m/P4/){ 
		print O2 $_ ;
		}
	}

while (<RP2>){
	if ($_ =~ m/1OP1/){ 
		print O2 $_ ;
		}
#	if ($_ =~ m/1OP1/){ 
#		print O2 $_ ;
#		}
	if ($_ =~ m/1OP5/){ 
		print O3 $_ ;
		}
	}

while (<LOOP2>){
	if ($_ =~ m/P5/){ 
		print O3 $_ ;
		}
	}

while (<RP3>){
	if ($_ =~ m/1OP1/){ 
		print O3 $_ ;
		}
	}

close (CCP);
close (RP);
close (RP2);
close (RP3);
close (LOOP1);
close (LOOP2);	
close (O1);
close (O2); 
close (O3);

system("$gromacs/$gmx confrms -nice 0 -one -f1 O3.pdb -f2 ../PIP2-GRO/PO1.pdb -name -o fit_PO1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O2.pdb -f2 ../PIP2-GRO/PO2.pdb -name -o fit_PO2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O1.pdb -f2 ../PIP2-GRO/PO3.pdb -name -o fit_PO3.pdb << EOD
1
1
EOD");

system("less fit_RP.pdb | grep 'C1'  >> complete_${loop_main}.pdb");
system("less fit_PO3.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C7'  >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C8'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v 'C14'| grep -v 'C18' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C21'| grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C27'  >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C33'| grep -v 'C37' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C41'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C2'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'O2'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C3'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'O3'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C4'  >> complete_${loop_main}.pdb");
system("less fit_PO1.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C5'  >> complete_${loop_main}.pdb");
system("less fit_PO2.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C6'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'O6'  >> complete_${loop_main}.pdb");
chdir "..";

$num_main++;
}

`cat fit_*/complete_*.pdb | grep 'ATOM' > pip2-bilayer.pdb`;
`cat pip2-bilayer.pdb >> bilayer-box.pdb`;

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ PIP2-GRO/`;
}

################ PIP2 GROMOS-CKP END ####################

#################### PIP3 GROMOS ########################

if ($pip3 > 0) {

`cp -R $cg2atdir/PIP3 .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (PIP3, "> pip3-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq PIP3)||($wod eq PI3)){
		if ($woc eq PO1){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"PBN",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO2){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"PBH",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO3){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"PAW",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO0){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"PBE",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP1){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C5 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP2){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C3 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq RP3){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"CA2",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print PIP3 $spa,$woa,$spb,$wob,$spc,"C50",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
	}
}
close(EMI);
close(PIP3);

# Directory to store all the best-fitting AT lipids, for concatenation into AT bilayer structure
`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f pip3-cg.pdb -o pip3-cg.pdb");

$num = 1;

while ($num <= $pip3){

$loop = sprintf("%04d", $num);

open (ATM, "pip3-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nPIP3 $pip3";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $pip3){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid

mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 PIP3/C4B.pdb -name -o fit_${loop_main}/fit_C4B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_CCP.pdb | grep 'OAV' > O1.pdb");
system("less fit_CCP.pdb | grep 'PAW' >> O1.pdb");
system("less fit_RP.pdb | grep 'OAZ' >> O1.pdb");
system("less fit_RP.pdb | grep 'OBM' >> O1.pdb");

system("less fit_RP.pdb | grep 'OBM' > O2.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'PBN' >> O2.pdb");
system("less fit_RP.pdb | grep 'OAZ' >> O2.pdb");

system("less fit_RP.pdb | grep 'OBG' > O3.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'PBH' >> O3.pdb");
system("less fit_RP.pdb | grep 'OAZ' >> O3.pdb");

system("less fit_RP.pdb | grep 'OBE' > O4.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'PBE' >> O4.pdb");
system("less fit_RP.pdb | grep 'OAZ' >> O4.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O1.pdb -f2 ../PIP3/PO1.pdb -name -o fit_PO1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O2.pdb -f2 ../PIP3/PO2.pdb -name -o fit_PO2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O3.pdb -f2 ../PIP3/PO3.pdb -name -o fit_PO3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 O4.pdb -f2 ../PIP3/PO0.pdb -name -o fit_PO0.pdb << EOD
1
1
EOD");

system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C46'  >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C42'  >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C38'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep -v 'O33'  >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' | grep 'O33'     >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'O33'   >> complete_${loop_main}.pdb");
system("less fit_PO1.pdb | grep 'OAX'  >> complete_${loop_main}.pdb");
system("less fit_PO1.pdb | grep 'OAY'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OAZ'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C1 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C6 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBS'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'HBS'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C5 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBM'  >> complete_${loop_main}.pdb");
system("less fit_PO2.pdb | grep 'ATOM' | grep -v 'OBM'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C4 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBG'  >> complete_${loop_main}.pdb");
system("less fit_PO3.pdb | grep 'ATOM' | grep -v 'OBG'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C3 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBE'  >> complete_${loop_main}.pdb");
system("less fit_PO0.pdb | grep 'PBE'  >> complete_${loop_main}.pdb");
system("less fit_PO0.pdb | grep 'OCA'  >> complete_${loop_main}.pdb");
system("less fit_PO0.pdb | grep 'OCB'  >> complete_${loop_main}.pdb");
system("less fit_PO0.pdb | grep 'OCC'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'C2 '  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'OBC'  >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'HBC'  >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C19' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep -v 'C23' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C26' >> complete_${loop_main}.pdb");
system("less fit_C4B.pdb | grep 'ATOM' | grep -v 'C29' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > pip3-bilayer.pdb");
system("cat pip3-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ PIP3/`;
}

################## PIP3 GROMOS END #####################
#################### DAG GROMOS ########################

if (($dag > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/DAG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DAG, "> dag-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq DAG){
		if ($woc eq DOH){	
			print DAG $spa,$woa,$spb,$wob,$spc,"O11",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DAG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DAG $spa,$woa,$spb,$wob,$spc,"O33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq D3B){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5B){	
			print DAG $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C38",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C42",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C46",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DAG $spa,$woa,$spb,$wob,$spc,"C49",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}

close(EMI);
close(DAG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dag-cg.pdb -o dag-cg.pdb");

$num = 1;

while ($num <= $dag){

$loop = sprintf("%04d", $num);

open (ATM, "dag-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDAG $dag";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dag){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DAG/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

`less fit_CCP.pdb | grep 'H1' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_GL2.pdb | grep 'ATOM' | grep -v O33 | grep -v C38 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C38' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C39' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C40' >> complete_${loop_main}.pdb`;
`less fit_C1B.pdb | grep 'C41' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C42' >> complete_${loop_main}.pdb`;
`less fit_C2B.pdb | grep 'ATOM' | grep -v C42 | grep -v C46 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C46' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C47' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C48' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C49' >> complete_${loop_main}.pdb`;
`less fit_C3B.pdb | grep 'C50' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dag-bilayer.pdb");
system("cat dag-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DAG/`;
}

######################DAG_END###########################

#################### MAG GROMOS ########################

if (($mag > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/MAG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (MAG, "> mag-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq MAG){
		if ($woc eq DOH){	
			print MAG $spa,$woa,$spb,$wob,$spc,"O11",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print MAG $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print MAG $spa,$woa,$spb,$wob,$spc,"C19",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print MAG $spa,$woa,$spb,$wob,$spc,"C23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print MAG $spa,$woa,$spb,$wob,$spc,"C26",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print MAG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C5A){	
			print MAG $spa,$woa,$spb,$wob,$spc,"CA1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}

close(EMI);
close(MAG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f mag-cg.pdb -o mag-cg.pdb");

$num = 1;

while ($num <= $mag){

$loop = sprintf("%04d", $num);

open (ATM, "mag-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nMAG $mag";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $mag){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 MAG/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 MAG/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 MAG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 MAG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 MAG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 MAG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

`less fit_CCP.pdb | grep 'H1A' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O14' >> complete_${loop_main}.pdb`;
`less fit_GL1.pdb | grep 'ATOM' | grep -v 'O14'| grep -v 'C19' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C19' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
`less fit_C1A.pdb | grep 'C22' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> complete_${loop_main}.pdb`;
`less fit_C2A.pdb | grep 'ATOM' | grep -v C23 | grep -v C26 >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C26' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C27' >> complete_${loop_main}.pdb`;
`less fit_C3A.pdb | grep 'C28' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C30' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'C31' >> complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'CA1' >> complete_${loop_main}.pdb`;
`less fit_C4A.pdb | grep 'CA2' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'C32' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'O33' >> complete_${loop_main}.pdb`;
`less fit_CCP.pdb | grep 'H1B' >> complete_${loop_main}.pdb`;

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > mag-bilayer.pdb");
system("cat mag-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ MAG/`;
}

######################MAG_END###########################

########################################################
#                                                      #
#                     STEROLS                          #
#                                                      #
########################################################

############## CHOLESTEROL GROMOS OLD ##################

if (($chol > 0)&&($ff ne "GROMOS96-53a6")) {

`cp -R $cg2atdir/CHOL .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (CHO, "> chol-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq CHO)||($wod eq CHOL)){
		if ($woc eq ROH){
			print CHO $spa,$woa,$spb,$wob,$spc,"C5  CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R1){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C5  CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R2){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C11 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R3){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C14 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R4){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C19 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R5){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C16 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C24 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C27 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(CHO);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f chol-cg.pdb -o chol-cg.pdb");

$num = 1;

while ($num <= $chol){

$loop = sprintf("%04d", $num);

open (ATM, "chol-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nCHOL $chol";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $chol){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CHOL/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_RP.pdb | grep 'ATOM' > chol-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C24' >> chol-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C27' >> chol-rings+tail.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 chol-rings+tail.pdb -f2 ../CHOL/C1A.pdb -name -o fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 chol-rings+tail.pdb -f2 ../CHOL/C2A.pdb -name -o fit_C2A.pdb << EOD
1
1
EOD");

system("less fit_RP.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C16'| grep -v 'C20'| grep -v 'C21'| grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > chol-bilayer.pdb");
system("cat chol-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ CHOL/`;
}

############# CHOLESTEROL GROMOS OLD END ###############
################# CHOLESTEROL GROMOS ###################

if (($chol > 0)&&($ff eq "GROMOS96-53a6")) {

`cp -R $cg2atdir/CHOL-GRO .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (CHO, "> chol-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq CHO)||($wod eq CHOL)){
		if ($woc eq ROH){
			print CHO $spa,$woa,$spb,$wob,$spc,"O1  CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R1){	
			print CHO $spa,$woa,$spb,$wob,$spc,"D10 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R2){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C7  CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R3){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C11 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R4){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C15 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq R5){	
			print CHO $spa,$woa,$spb,$wob,$spc,"D13 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C22 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print CHO $spa,$woa,$spb,$wob,$spc,"C25 CHO",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(CHO);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f chol-cg.pdb -o chol-cg.pdb");

$num = 1;

while ($num <= $chol){

$loop = sprintf("%04d", $num);

open (ATM, "chol-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nCHOL $chol";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $chol){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 CHOL-GRO/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_RP.pdb | grep 'ATOM' | grep -v 'D13'| grep -v 'D10' > chol-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C22' >> chol-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C25' >> chol-rings+tail.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 chol-rings+tail.pdb -f2 ../CHOL-GRO/C1A.pdb -name -o fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 chol-rings+tail.pdb -f2 ../CHOL-GRO/C2A.pdb -name -o fit_C2A.pdb << EOD
1
1
EOD");

system("less fit_RP.pdb | grep -v 'D13'| grep -v 'D10' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' |  grep -v 'C17'| grep -v 'C22' | grep -v 'C16'>> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > chol-bilayer.pdb");
system("cat chol-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ CHOL-GRO/`;
}

############## CHOLESTEROL GROMOS END ##################

########################################################
#                                                      #
#                    Detergents                        #
#                                                      #
########################################################

#################### DHPC GROMOS #######################

if (($dhpc > 0)&&($ff ne "OPLSAA")) {

`cp -R $cg2atdir/DHPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DHPC, "> dhpc-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DHP)||($wod eq DHPC)){
		if ($woc eq NC3){
			print DHPC $spa,$woa,$spb,$wob,$spc,"N4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"P8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"O14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"O23",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C27",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C30",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C18",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(DHPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dhpc-cg.pdb -o dhpc-cg.pdb");

$num = 1;

while ($num <= $dhpc){

$loop = sprintf("%04d", $num);

open (ATM, "dhpc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDHPC $dhpc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dhpc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DHPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep 'C1' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C2' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C3' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'N4' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C5' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'C6' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep 'P8' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'O11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C19' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'C21' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep 'C22' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C28' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C29' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'C30' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dhpc-bilayer.pdb");
system("cat dhpc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DHPC/`;
}

#################### DHPC GROMOS END #######################
###################### DHPC OPLSAA #########################

if (($dhpc > 0)&&($ff eq "OPLSAA")) {

`cp -R $cg2atdir/DHPC-oplsaa .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DHPC, "> dhpc-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq DHP)||($wod eq DHPC)){
		if ($woc eq NC3){
			print DHPC $spa,$woa,$spb,$wob,$spc,"N1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"P1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"O3 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"O2 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C4 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C1 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C11",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DHPC $spa,$woa,$spb,$wob,$spc,"C14",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(DHPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dhpc-cg.pdb -o dhpc-cg.pdb");

$num = 1;

while ($num <= $dhpc){

$loop = sprintf("%04d", $num);

open (ATM, "dhpc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDHPC $dhpc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dhpc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC-oplsaa/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC-oplsaa/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC-oplsaa/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC-oplsaa/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC-oplsaa/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DHPC-oplsaa/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_CCP.pdb | grep 'O5' > headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'P1' >> headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'O8' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DHPC-oplsaa/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_C1B.pdb | grep 'ATOM' | grep -v 'C4' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM' | grep 'C4' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'C4' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'C15' | grep -v 'H26'| grep -v 'H27'| grep -v 'O5'| grep -v 'P1'|grep -v 'O2'>> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep -v 'O3' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v 'C11' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'C7' | grep -v 'H12'| grep -v 'H13'| grep -v 'C8'| grep -v 'H14'|grep -v 'O3'|grep -v 'O2'>> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O6' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O7' >> complete_${loop_main}.pdb");
system("less fit_NC3.pdb | grep -v 'P1' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dhpc-bilayer.pdb");
system("cat dhpc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DHPC-oplsaa/`;
}

################## DHPC OPLS END ######################
################### BOG GROMOS ########################

if (($bog > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/BOG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (BOG, "> bog-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq BOG){
		if ($woc eq B1){
			print BOG $spa,$woa,$spb,$wob,$spc,"C5",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B3){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C4B ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C8B ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(BOG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f bog-cg.pdb -o bog-cg.pdb");

$num = 1;

while ($num <= $bog){

$loop = sprintf("%04d", $num);

open (ATM, "bog-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nBOG $bog";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $bog){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 BOG/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 BOG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 BOG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C4B' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM'| grep -v 'C3 ' >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'ATOM' | grep 'O1' >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'ATOM' | grep 'C1' >> complete_${loop_main}.pdb");
system("less fit_RP.pdb | grep 'ATOM' | grep -v 'O1' | grep -v 'C1' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > bog-bilayer.pdb");
system("cat bog-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ BOG/`;
}

################## BOG GROMOS END ######################
#################### BOG OPLSAA ########################

if (($bog > 0)&&($ff eq "OPLSAA")) {

`cp -R $cg2atdir/BOG-oplsaa .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (BOG, "> bog-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq BOG){
		if ($woc eq B1){
			print BOG $spa,$woa,$spb,$wob,$spc,"C5",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B3){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C10 ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print BOG $spa,$woa,$spb,$wob,$spc,"C14 ",$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(BOG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f bog-cg.pdb -o bog-cg.pdb");

$num = 1;

while ($num <= $bog){

$loop = sprintf("%04d", $num);

open (ATM, "bog-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nBOG $bog";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $bog){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 BOG-oplsaa/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 BOG-oplsaa/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 BOG-oplsaa/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_RP.pdb | grep 'ATOM' | grep -v 'C1'  >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' |  grep -v 'C10' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > bog-bilayer.pdb");
system("cat bog-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
#`rm '#'*`;
#`rm -rf Res_CG/ fit_*/ BOG-oplsaa/`;
}

###################### BOG OPLSAA END ########################
######################## DDM OPLSAA ##########################

if (($ddm > 0)&&($ff eq "OPLSAA")) {

`cp -R $cg2atdir/DDM .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DDM, "> ddm-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq DDM){
		if ($woc eq B1){
			print DDM $spa,$woa,$spb,$wob,$spc,"C1  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C6  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B3){	
			print DDM $spa,$woa,$spb,$wob,$spc,"O4  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B4){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C12 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B5){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C9  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B6){	
			print DDM $spa,$woa,$spb,$wob,$spc,"O9  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C16 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C20 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C23 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(DDM);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f ddm-cg.pdb -o ddm-cg.pdb");

$num = 1;

while ($num <= $ddm){

$loop = sprintf("%04d", $num);

open (ATM, "ddm-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDDM $ddm";
close (LIP);

open (DDM, ">> ddm.top") || die "No output possible\n";
print DDM "\nDDM $ddm";
close (DDM);

# For every CG lipid
$num_main = 1;
while ($num_main <= $ddm){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DDM/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";
system("cat fit_RP.pdb | grep 'ATOM' > ddm-rings+tail.pdb");

system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C16' >> ddm-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C20' >> ddm-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C23' >> ddm-rings+tail.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 ddm-rings+tail.pdb -f2 ../DDM/C1A.pdb -name -o fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 ddm-rings+tail.pdb -f2 ../DDM/C2A.pdb -name -o fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 ddm-rings+tail.pdb -f2 ../DDM/C3A.pdb -name -o fit_C3A.pdb << EOD
1
1
EOD");

system("less fit_RP.pdb | grep 'ATOM'>> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'O11' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C16' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C20' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > ddm-bilayer.pdb");
system("cat box.pdb ddm-bilayer.pdb > ddm-box.pdb");
system("$gromacs/$gmx grompp -f em-cg.mdp -c ddm-box.pdb -p ddm.top -o ddm_em -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm ddm_em");
system("$gromacs/$gmx editconf -resnr 1 -f ddm_em.gro -o ddm_em.pdb");
system("cat ddm_em.pdb | grep ATOM >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DDM/`;
}

###################### DDM OPLSAA END #########################
######################## DDM GROMOS ###########################

if (($ddm > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DDM-GRO .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DDM, "> ddm-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq DDM){
		if ($woc eq B1){
			print DDM $spa,$woa,$spb,$wob,$spc,"DU1 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C6  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B3){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C1  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B4){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C12 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B5){	
			print DDM $spa,$woa,$spb,$wob,$spc,"O10 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B6){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C9  DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C16 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C20 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3){	
			print DDM $spa,$woa,$spb,$wob,$spc,"C24 DDM",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(DDM);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f ddm-cg.pdb -o ddm-cg.pdb");

$num = 1;

while ($num <= $ddm){

$loop = sprintf("%04d", $num);

open (ATM, "ddm-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDDM $ddm";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $ddm){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DDM-GRO/RP1.pdb -name -o fit_${loop_main}/fit_RP1.pdb << EOD
1
1
EOD");

system("less Res_CG/res_cg_${loop_main}.pdb | grep 'C12' > fit_${loop_main}/ddm_rp1.pdb");
system("less Res_CG/res_cg_${loop_main}.pdb | grep 'O10' >> fit_${loop_main}/ddm_rp1.pdb");
system("less Res_CG/res_cg_${loop_main}.pdb | grep 'C9' >> fit_${loop_main}/ddm_rp1.pdb");
system("cat fit_${loop_main}/fit_RP1.pdb >> fit_${loop_main}/ddm_rp1.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 fit_${loop_main}/ddm_rp1.pdb -f2 DDM-GRO/RP2.pdb -name -o fit_${loop_main}/fit_RP2.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";


system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C24' > ddm-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C20' >> ddm-rings+tail.pdb");
system("less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C16' >> ddm-rings+tail.pdb");
system("cat fit_RP1.pdb fit_RP2.pdb | grep -v DUM | grep 'ATOM' >> ddm-rings+tail.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 ddm-rings+tail.pdb -f2 ../DDM-GRO/C1A.pdb -name -o fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 ddm-rings+tail.pdb -f2 ../DDM-GRO/C2A.pdb -name -o fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 ddm-rings+tail.pdb -f2 ../DDM-GRO/C3A.pdb -name -o fit_C3A.pdb << EOD
1
1
EOD");



system("less fit_C3A.pdb | grep 'ATOM' | grep 'C24' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C24' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C16'| grep -v 'C20' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM'  >> complete_${loop_main}.pdb");

system("less fit_RP2.pdb | grep 'ATOM'| grep -v 'O11' | grep -v ' DU2' >> complete_${loop_main}.pdb");
system("less fit_RP1.pdb | grep 'ATOM'| grep -v ' DU1'  | grep -v ' O6 '>> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > ddm-bilayer.pdb");
system("cat ddm-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DDM-GRO/`;
}

################# DDM GROMOS END ######################

######################## LMPG GROMOS ###########################

if (($lmpg > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/LMPG .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (LMPG, "> lmpg-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($wod eq LMP)||($wod eq LMPG)){
		if ($woc eq GLH){
			print LMPG $spa,$woa,$spb,$wob,$spc,"O2  LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print LMPG $spa,$woa,$spb,$wob,$spc,"P9  LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print LMPG $spa,$woa,$spb,$wob,$spc,"C17 LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print LMPG $spa,$woa,$spb,$wob,$spc,"C21 LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print LMPG $spa,$woa,$spb,$wob,$spc,"C25 LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print LMPG $spa,$woa,$spb,$wob,$spc,"C29 LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print LMPG $spa,$woa,$spb,$wob,$spc,"C33 LMP",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(LMPG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f lmpg-cg.pdb -o lmpg-cg.pdb");

$num = 1;

while ($num <= $lmpg){

$loop = sprintf("%04d", $num);

open (ATM, "lmpg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nLMPG $lmpg";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $lmpg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMPG/P1A.pdb -name -o fit_${loop_main}/fit_P1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMPG/P2A.pdb -name -o fit_${loop_main}/fit_P2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMPG/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMPG/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMPG/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMPG/C4A.pdb -name -o fit_${loop_main}/fit_C4A.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_P1A.pdb | grep 'O8' > headgroup_chemistry.pdb");
system("less fit_P1A.pdb | grep 'P9' >> headgroup_chemistry.pdb");
system("less fit_P2A.pdb | grep 'O12' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../LMPG/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_P1A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O11'  >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O12'  >> complete_${loop_main}.pdb");
system("less fit_P2A.pdb | grep 'ATOM' | grep -v 'P9' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep 'ATOM' | grep -v 'C17' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' | grep -v 'C21' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep 'ATOM' | grep -v 'C25' >> complete_${loop_main}.pdb");
system("less fit_C4A.pdb | grep 'ATOM' | grep -v 'C29' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > lmpg-bilayer.pdb");
system("cat lmpg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ LMPG/`;
}

################# LMPG GROMOS END ######################

#################### DPC OPLS #########################

if (($dpc > 0)&&($ff =~ m/^OPLS/)) {

`cp -R $cg2atdir/DPC .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPC, "> dpc-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq DPC){
		if ($woc eq NC3){
			print DPC $spa,$woa,$spb,$wob,$spc,"N1  DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DPC $spa,$woa,$spb,$wob,$spc,"P1  DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print DPC $spa,$woa,$spb,$wob,$spc,"C9  DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print DPC $spa,$woa,$spb,$wob,$spc,"C13 DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3){	
			print DPC $spa,$woa,$spb,$wob,$spc,"C17 DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(DPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dpc-cg.pdb -o dpc-cg.pdb");

$num = 1;

while ($num <= $dpc){

$loop = sprintf("%04d", $num);

open (ATM, "dpc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPC $dpc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dpc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O1' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P1' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O4' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPC/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb | grep -v 'P1' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O2' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O3' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'C9' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v 'C13' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dpc-bilayer.pdb");
system("cat dpc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPC/`;
}

################### DPC OPLS END ##########################
#################### DPC GROMOS ###########################

if (($dpc > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/DPC-GRO .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DPC, "> dpc-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq DPC){
		if ($woc eq NC3){
			print DPC $spa,$woa,$spb,$wob,$spc,"N4  DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq PO4){	
			print DPC $spa,$woa,$spb,$wob,$spc,"P8  DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1){	
			print DPC $spa,$woa,$spb,$wob,$spc,"C14 DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2){	
			print DPC $spa,$woa,$spb,$wob,$spc,"C18 DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3){	
			print DPC $spa,$woa,$spb,$wob,$spc,"C22 DPC",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(DPC);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dpc-cg.pdb -o dpc-cg.pdb");

$num = 1;

while ($num <= $dpc){

$loop = sprintf("%04d", $num);

open (ATM, "dpc-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDPC $dpc";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dpc){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC-GRO/NC3.pdb -name -o fit_${loop_main}/fit_NC3.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC-GRO/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC-GRO/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DPC-GRO/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");


chdir "fit_${loop_main}";

system("less fit_NC3.pdb | grep 'O7' > headgroup_chemistry.pdb");
system("less fit_NC3.pdb | grep 'P8' >> headgroup_chemistry.pdb");
system("less fit_CCP.pdb | grep 'O11' >> headgroup_chemistry.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 headgroup_chemistry.pdb -f2 ../DPC-GRO/PO4.pdb -name -o fit_PO4.pdb << EOD
1
1
EOD");

system("less fit_NC3.pdb  >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O9' >> complete_${loop_main}.pdb");
system("less fit_PO4.pdb | grep 'O10' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep -v 'C14' | grep -v 'P8' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v 'C18' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep 'ATOM' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > dpc-bilayer.pdb");
system("cat dpc-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DPC-GRO/`;
}

################## DPC GROMOS END #####################

#################### SDS GROMOS ###########################

if (($sds > 0)&&($ff =~ m/^GROMOS/)) {

`cp -R $cg2atdir/SDS .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (SDS, "> sds-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq SDS){
		if ($woc eq GL1){
			print SDS $spa,$woa,$spb,$wob,$spc,"DAO SDS",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print SDS $spa,$woa,$spb,$wob,$spc,"DAK SDS",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print SDS $spa,$woa,$spb,$wob,$spc,"DAG SDS",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print SDS $spa,$woa,$spb,$wob,$spc,"CAC SDS",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		}
	}
	
close(EMI);
close(SDS);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f sds-cg.pdb -o sds-cg.pdb");

$num = 1;

while ($num <= $sds){

$loop = sprintf("%04d", $num);

open (ATM, "sds-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nSDS $sds";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $sds){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SDS/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");


system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SDS/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SDS/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";

system("less fit_C2A.pdb | grep -v 'DAG'  >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v 'DAK'| grep -v 'DAG'>> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep -v 'DAO'  | grep -v 'DAK' >> complete_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/complete_*.pdb | grep 'ATOM' > sds-bilayer.pdb");
system("cat sds-bilayer.pdb >> bilayer-box.pdb");


# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ SDS/`;
}

################## SDS GROMOS END #####################

#################### LMG AMBER ########################

if (($lmg > 0)&&($ff eq "AMBER")) {

`cp -R $cg2atdir/LMG-AMBER .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (LMG, "> lmg-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq LMG){
		if ($woc eq B0){
			print LMG $spa,$woa,$spb,$wob,$spc,"D3",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print LMG $spa,$woa,$spb,$wob,$spc,"O1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}		
		if ($woc eq B3){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C6",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print LMG $spa,$woa,$spb,$wob,$spc,"D7 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print LMG $spa,$woa,$spb,$wob,$spc,"O8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C13",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C17",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C25",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C37",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print LMG $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}

close(EMI);
close(LMG);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f lmg-cg.pdb -o lmg-cg.pdb");

$num = 1;

while ($num <= $lmg){

$loop = sprintf("%04d", $num);

open (ATM, "lmg-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nLMG $lmg";
close (LIP);

open (LIP, ">> lmg.top") || die "No output possible\n";
print LIP "\nLMG 1";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $lmg){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 LMG-AMBER/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";


system("less fit_GL1.pdb | grep 'C12' >> C1A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C13' >> C1A-link.pdb`;
system("less fit_C1A.pdb | grep 'C14' >> C1A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1A-link.pdb -f2 ../LMG-AMBER/C1AH.pdb -name -o fit_C1AH.pdb << EOD
1
1
EOD");

system("less fit_C1A.pdb | grep 'C16' >> C2A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C17' >> C2A-link.pdb`;
system("less fit_C2A.pdb | grep 'C18' >> C2A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2A-link.pdb -f2 ../LMG-AMBER/C2AH.pdb -name -o fit_C2AH.pdb << EOD
1
1
EOD");

system("less fit_C2A.pdb | grep 'C20' >> C3A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C21' >> C3A-link.pdb`;
system("less fit_C3A.pdb | grep 'C22' >> C3A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3A-link.pdb -f2 ../LMG-AMBER/C3AH.pdb -name -o fit_C3AH.pdb << EOD
1
1
EOD");

system("less fit_GL2.pdb | grep 'C28' >> C1B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> C1B-link.pdb`;
system("less fit_C1B.pdb | grep 'C30' >> C1B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1B-link.pdb -f2 ../LMG-AMBER/C1BH.pdb -name -o fit_C1BH.pdb << EOD
1
1
EOD");

system("less fit_C1B.pdb | grep 'C32' >> C2B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> C2B-link.pdb`;
system("less fit_C2B.pdb | grep 'C34' >> C2B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2B-link.pdb -f2 ../LMG-AMBER/C2BH.pdb -name -o fit_C2BH.pdb << EOD
1
1
EOD");

system("less fit_C2B.pdb | grep 'C36' >> C3B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> C3B-link.pdb`;
system("less fit_C3B.pdb | grep 'C38' >> C3B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3B-link.pdb -f2 ../LMG-AMBER/C3BH.pdb -name -o fit_C3BH.pdb << EOD
1
1
EOD");

`cat ../box.pdb > complete_${loop_main}.pdb`;
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D7' >> complete_${loop_main}.pdb`;
system("less fit_CCP.pdb | grep ' C8' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C7' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O1'| grep 'ATOM' >> complete_${loop_main}.pdb`;
system("less fit_RP.pdb | grep -v 'D3'| grep -v 'O1'| grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H46' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C9' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O8'| grep 'ATOM' >> complete_${loop_main}.pdb`;
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'O8'|grep -v ' H' | grep -v 'C29'>> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
system("less fit_C1B.pdb | grep -v ' H' | grep 'ATOM'| grep -v 'C29' | grep -v 'C33' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
system("less fit_C2B.pdb | grep -v ' H' | grep 'ATOM'| grep -v 'C33' | grep -v 'C37' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> complete_${loop_main}.pdb`;
system("less fit_C3B.pdb | grep -v ' H' | grep 'ATOM' | grep -v 'C37' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3BH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C2BH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1BH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep ' H'  |grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H14' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H47' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H13' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb |grep -v ' D7' |grep -v ' O7'| grep -v ' C13' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
system("less fit_C1AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v ' C13' |  grep -v ' C17'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C17' >> complete_${loop_main}.pdb`;
system("less fit_C2AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep -v ' C17'| grep -v ' C21' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
system("less fit_C3AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep -v ' C21' |grep 'ATOM' >> complete_${loop_main}.pdb");

system("$gromacs/$gmx grompp -f ../em-lipid1.mdp -c complete_${loop_main}.pdb -p ../lmg.top -o lmg_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm lmg_${loop_main}");
system("$gromacs/$gmx grompp -f ../em-lipid2.mdp -c lmg_${loop_main}.gro -p ../lmg.top -o lmg_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm lmg_${loop_main}");
system("$gromacs/$gmx editconf -resnr 1 -f lmg_${loop_main}.gro -o lmg_${loop_main}.pdb");

chdir "..";

$num_main++;
}


system("cat fit_*/lmg_*.pdb | grep 'ATOM' > lmg-bilayer.pdb");
system("cat lmg-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ LMG-AMBER/`;
}

#################### LMG AMBER END ####################
###################### SQD AMBER ######################

if (($sqd > 0)&&($ff eq "AMBER")) {

`cp -R $cg2atdir/SQD-AMBER .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (SQD, "> sqd-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq SQD){
		if ($woc eq B1){
			print SQD $spa,$woa,$spb,$wob,$spc,"D50 SQD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print SQD $spa,$woa,$spb,$wob,$spc,"O1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B3){
			print SQD $spa,$woa,$spb,$wob,$spc,"D43 SQD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq S1){	
			print SQD $spa,$woa,$spb,$wob,$spc,"S1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print SQD $spa,$woa,$spb,$wob,$spc,"D7 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print SQD $spa,$woa,$spb,$wob,$spc,"O8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print SQD $spa,$woa,$spb,$wob,$spc,"D13",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print SQD $spa,$woa,$spb,$wob,$spc,"D17",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print SQD $spa,$woa,$spb,$wob,$spc,"D21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print SQD $spa,$woa,$spb,$wob,$spc,"D25",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print SQD $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print SQD $spa,$woa,$spb,$wob,$spc,"C33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print SQD $spa,$woa,$spb,$wob,$spc,"C37",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print SQD $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}

close(EMI);
close(SQD);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f sqd-cg.pdb -o sqd-cg.pdb");

$num = 1;

while ($num <= $sqd){

$loop = sprintf("%04d", $num);

open (ATM, "sqd-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nSQD $sqd";
close (LIP);

open (LIP, ">> sqd.top") || die "No output possible\n";
print LIP "\nSQD 1";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $sqd){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 SQD-AMBER/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";
system("less fit_C1A.pdb | grep 'C14' >> C1A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D13' >> C1A-link.pdb`;
system("less fit_GL1.pdb | grep 'C12' >> C1A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1A-link.pdb -f2 ../SQD-AMBER/C1AH.pdb -name -o fit_C1AH.pdb << EOD
1
1
EOD");

system("less fit_C2A.pdb | grep 'C18' >> C2A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D17' >> C2A-link.pdb`;
system("less fit_C1A.pdb | grep 'C16' >> C2A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2A-link.pdb -f2 ../SQD-AMBER/C2AH.pdb -name -o fit_C2AH.pdb << EOD
1
1
EOD");

system("less fit_C3A.pdb | grep 'C22' >> C3A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D21' >> C3A-link.pdb`;
system("less fit_C2A.pdb | grep 'C20' >> C3A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3A-link.pdb -f2 ../SQD-AMBER/C3AH.pdb -name -o fit_C3AH.pdb << EOD
1
1
EOD");

system("less fit_GL2.pdb | grep 'C28' >> C1B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> C1B-link.pdb`;
system("less fit_C1B.pdb | grep 'C30' >> C1B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1B-link.pdb -f2 ../SQD-AMBER/C1BH.pdb -name -o fit_C1BH.pdb << EOD
1
1
EOD");

system("less fit_C1B.pdb | grep 'C32' >> C2B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> C2B-link.pdb`;
system("less fit_C2B.pdb | grep 'C34' >> C2B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2B-link.pdb -f2 ../SQD-AMBER/C2BH.pdb -name -o fit_C2BH.pdb << EOD
1
1
EOD");

system("less fit_C2B.pdb | grep 'C36' >> C3B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> C3B-link.pdb`;
system("less fit_C3B.pdb | grep 'C38' >> C3B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3B-link.pdb -f2 ../SQD-AMBER/C3BH.pdb -name -o fit_C3BH.pdb << EOD
1
1
EOD");

`cat ../box.pdb > complete_${loop_main}.pdb`;

system("less fit_C3A.pdb | grep -v ' C21' | grep -v ' D'|grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D21' >> complete_${loop_main}.pdb`;
system("less fit_C3AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep -v ' C17'| grep -v ' D'| grep -v ' C21' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D17' >> complete_${loop_main}.pdb`;
system("less fit_C2AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep -v ' C13' |  grep -v ' C17'| grep -v ' D'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D13' >> complete_${loop_main}.pdb`;
system("less fit_C1AH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb |grep -v ' D' | grep -v ' O7'| grep -v ' C13' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'D7' >> complete_${loop_main}.pdb`;
system("less fit_CCP.pdb | grep ' C8' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H13' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C7' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H46' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H12' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O1'| grep 'ATOM' >> complete_${loop_main}.pdb`;
system("less fit_RP.pdb | grep -v 'D43'| grep -v 'O1' | grep -v 'D50' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C9' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H14' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H47' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O8'| grep 'ATOM' >> complete_${loop_main}.pdb`;
system("less fit_GL2.pdb | grep 'ATOM' | grep -v 'O8' | grep -v 'C29'>> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
system("less fit_C1BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep 'ATOM'| grep -v 'C29' | grep -v 'C33' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
system("less fit_C2BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM'| grep -v 'C33' | grep -v 'C37' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> complete_${loop_main}.pdb`;
system("less fit_C3BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep -v 'C37' >> complete_${loop_main}.pdb");

system("$gromacs/$gmx grompp -f ../em-lipid1.mdp -c complete_${loop_main}.pdb -p ../sqd.top -o sqd_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm sqd_${loop_main}");
system("$gromacs/$gmx grompp -f ../em-lipid2.mdp -c sqd_${loop_main}.gro -p ../sqd.top -o sqd_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm sqd_${loop_main}");
system("$gromacs/$gmx editconf -resnr 1 -f sqd_${loop_main}.gro -o sqd_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/sqd_*.pdb | grep 'ATOM' > sqd-bilayer.pdb");
system("cat sqd-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ SQD-AMBER/`;
}

#################### SQD AMBER END ####################
###################### DGD AMBER ######################

if (($dgd > 0)&&($ff eq "AMBER")) {

`cp -R $cg2atdir/DGD-AMBER .`;

open (EMI, "em_inbox.pdb") || die "Cannot open em_inbox.pdb\n";
open (DGD, "> dgd-cg.pdb") || die "No output possible\n";

while (<EMI>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($wod eq DGD){
		if ($woc eq B1){
			print DGD $spa,$woa,$spb,$wob,$spc,"D72 DGD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B2){	
			print DGD $spa,$woa,$spb,$wob,$spc,"O1",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B3){
			print DGD $spa,$woa,$spb,$wob,$spc,"C63 DGD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B4){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C54 DGD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B5){
			print DGD $spa,$woa,$spb,$wob,$spc,"D73 DGD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq B6){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C64 DGD",$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL1){	
			print DGD $spa,$woa,$spb,$wob,$spc,"D7 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq GL2){	
			print DGD $spa,$woa,$spb,$wob,$spc,"O8 ",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1A){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C13",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2A){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C17",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3A){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C21",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4A){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C25",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C1B){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C29",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C2B){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C33",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C3B){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C37",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		if ($woc eq C4B){	
			print DGD $spa,$woa,$spb,$wob,$spc,"C41",$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}			
		}
	}

close(EMI);
close(DGD);

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');
system("$gromacs/$gmx editconf -resnr 1 -f dgd-cg.pdb -o dgd-cg.pdb");

$num = 1;

while ($num <= $dgd){

$loop = sprintf("%04d", $num);

open (ATM, "dgd-cg.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $num){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);
close(RESCG);

$num++;
}

open (LIP, ">> lipid.top") || die "No output possible\n";
print LIP "\nDGD $dgd";
close (LIP);

open (LIP, ">> dgd.top") || die "No output possible\n";
print LIP "\nDGD 1";
close (LIP);

# For every CG lipid
$num_main = 1;
while ($num_main <= $dgd){

$loop_main = sprintf("%04d", $num_main);

# Temporary directory to store all fitted AT lipids in library to the current CG lipid
mkdir("fit_${loop_main}");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/RP.pdb -name -o fit_${loop_main}/fit_RP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/CCP.pdb -name -o fit_${loop_main}/fit_CCP.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/GL1.pdb -name -o fit_${loop_main}/fit_GL1.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/GL2.pdb -name -o fit_${loop_main}/fit_GL2.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/C1A.pdb -name -o fit_${loop_main}/fit_C1A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/C2A.pdb -name -o fit_${loop_main}/fit_C2A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/C3A.pdb -name -o fit_${loop_main}/fit_C3A.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/C1B.pdb -name -o fit_${loop_main}/fit_C1B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/C2B.pdb -name -o fit_${loop_main}/fit_C2B.pdb << EOD
1
1
EOD");

system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 DGD-AMBER/C3B.pdb -name -o fit_${loop_main}/fit_C3B.pdb << EOD
1
1
EOD");

chdir "fit_${loop_main}";


system("less fit_GL1.pdb | grep 'C12' >> C1A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C13' >> C1A-link.pdb`;
system("less fit_C1A.pdb | grep 'C14' >> C1A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1A-link.pdb -f2 ../DGD-AMBER/C1AH.pdb -name -o fit_C1AH.pdb << EOD
1
1
EOD");

system("less fit_C1A.pdb | grep 'C16' >> C2A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C17' >> C2A-link.pdb`;
system("less fit_C2A.pdb | grep 'C18' >> C2A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2A-link.pdb -f2 ../DGD-AMBER/C2AH.pdb -name -o fit_C2AH.pdb << EOD
1
1
EOD");

system("less fit_C2A.pdb | grep 'C20' >> C3A-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C21' >> C3A-link.pdb`;
system("less fit_C3A.pdb | grep 'C22' >> C3A-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3A-link.pdb -f2 ../DGD-AMBER/C3AH.pdb -name -o fit_C3AH.pdb << EOD
1
1
EOD");

system("less fit_GL2.pdb | grep 'C28' >> C1B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> C1B-link.pdb`;
system("less fit_C1B.pdb | grep 'C30' >> C1B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C1B-link.pdb -f2 ../DGD-AMBER/C1BH.pdb -name -o fit_C1BH.pdb << EOD
1
1
EOD");

system("less fit_C1B.pdb | grep 'C32' >> C2B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> C2B-link.pdb`;
system("less fit_C2B.pdb | grep 'C34' >> C2B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C2B-link.pdb -f2 ../DGD-AMBER/C2BH.pdb -name -o fit_C2BH.pdb << EOD
1
1
EOD");

system("less fit_C2B.pdb | grep 'C36' >> C3B-link.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> C3B-link.pdb`;
system("less fit_C3B.pdb | grep 'C38' >> C3B-link.pdb");

system("$gromacs/$gmx confrms -nice 0 -one -f1 C3B-link.pdb -f2 ../DGD-AMBER/C3BH.pdb -name -o fit_C3BH.pdb << EOD
1
1
EOD");

`cat ../box.pdb > complete_${loop_main}.pdb`;
system("less fit_GL2.pdb | grep ' C26' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep ' O10' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O8'| grep 'ATOM' >> complete_${loop_main}.pdb`;
system("less fit_CCP.pdb | grep ' C9' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C8' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep ' D7 '| grep 'ATOM' >> complete_${loop_main}.pdb`;

system("less fit_GL1.pdb | grep 'ATOM'| grep -v ' O7'| grep -v ' D7' | grep -v 'C13' |grep -v ' H' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C13' >> complete_${loop_main}.pdb`;
system("less fit_C1A.pdb | grep -v ' H' | grep 'ATOM'| grep -v 'C13' | grep -v 'C17' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C17' >> complete_${loop_main}.pdb`;
system("less fit_C2A.pdb | grep -v ' H' | grep 'ATOM'| grep -v 'C21' | grep -v 'C17' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C21' >> complete_${loop_main}.pdb`;
system("less fit_C3A.pdb | grep -v ' H' | grep 'ATOM' | grep -v 'C21' >> complete_${loop_main}.pdb");
system("less fit_C3A.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C3AH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C2A.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C2AH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1A.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_C1AH.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_GL1.pdb | grep ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' C7' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'O1' >> complete_${loop_main}.pdb`;
system("less fit_RP.pdb | grep -v ' O1' | grep -v 'D64' | grep -v 'D63' |grep 'ATOM' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H46' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H12' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H13' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H14' >> complete_${loop_main}.pdb");
system("less fit_CCP.pdb | grep ' H47' >> complete_${loop_main}.pdb");
system("less fit_GL2.pdb | grep 'ATOM'| grep -v 'O8' | grep -v 'C29' | grep -v 'O10' | grep -v 'C26' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C29' >> complete_${loop_main}.pdb`;
system("less fit_C1BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C1B.pdb | grep -v ' C33' |  grep -v ' C29'| grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C33' >> complete_${loop_main}.pdb`;
system("less fit_C2BH.pdb | grep 'ATOM'| grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep -v ' C37'| grep -v ' C33' | grep -v ' H' | grep 'ATOM' >> complete_${loop_main}.pdb");
`less ../Res_CG/res_cg_${loop_main}.pdb | grep 'C37' >> complete_${loop_main}.pdb`;
system("less fit_C3B.pdb | grep -v ' C37' |grep 'ATOM'| grep -v ' H' >> complete_${loop_main}.pdb");
system("less fit_C3B.pdb | grep 'ATOM' | grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C3BH.pdb | grep 'ATOM' | grep ' H' >> complete_${loop_main}.pdb");
system("less fit_C2B.pdb | grep 'ATOM' | grep ' H' >> complete_${loop_main}.pdb");

system("$gromacs/$gmx grompp -f ../em-lipid1.mdp -c complete_${loop_main}.pdb -p ../dgd.top -o dgd_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dgd_${loop_main}");
system("$gromacs/$gmx grompp -f ../em-lipid2.mdp -c dgd_${loop_main}.gro -p ../dgd.top -o dgd_${loop_main} -maxwarn 5"); 
system("$gromacs/$gmx mdrun -deffnm dgd_${loop_main}");
system("$gromacs/$gmx editconf -resnr 1 -f dgd_${loop_main}.gro -o dgd_${loop_main}.pdb");

chdir "..";

$num_main++;
}

system("cat fit_*/dgd_*.pdb | grep 'ATOM' > dgd-bilayer.pdb");
system("cat dgd-bilayer.pdb >> bilayer-box.pdb");

# Cleaning up...
`rm '#'*`;
`rm -rf Res_CG/ fit_*/ DGD-AMBER/`;
}

################### DGD AMBER END ##################

####################################################
#                                                  #
#                LIPID FRAGMENTS END               #
#                                                  #
####################################################

if ($lipidnum > 0) {
system("$gromacs/$gmx grompp -maxwarn 10 -f em.mdp -o lipid_em -c bilayer-box.pdb -p lipid.top");
system("$gromacs/$gmx mdrun -deffnm lipid_em");
system("$gromacs/$gmx trjconv -f lipid_em.gro -s lipid_em.tpr -center -pbc res -o atomistic-lipid.pdb << EOD
1
0
EOD");
open ATM, "<", 'atomistic-lipid.pdb' or die "Energy minimisation of AT lipid failed\n";
close(ATM);
}
}
###################################################
#                                                 #
#               PROTEIN CONVERSION                #
#                                                 #
################################################### 

if ($protnum > 0) {

`cp -R $cg2atdir/Sidechains .`;

system("$gromacs/$gmx make_ndx -f em_inbox.pdb -o full.ndx << EOD
del 0
del 1-30
q
EOD");

system("$gromacs/$gmx editconf -resnr 1 -f em_inbox.pdb -o protein.pdb -n full.ndx");

open PRO,"<", 'protein.pdb' or die;
open ATM, ">", 'atomistic-particles.pdb' or die;

if ($forcefield eq "martini") { 

$count=1;

while (<PRO>){
		chomp;
		($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
		($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
		if ($woa eq "ATOM"){
#### Newer Martini ####
		 if ($woc=~ m/^BB/){$woc="CA";
		 	} elsif (($woc=~ m/^SC1/)&&($wod eq "CYS")){$woc="SG ";
		 	} elsif (($woc=~ m/^SC1/)&&(($wod eq "GLU")||($wod eq "GLN"))) {$woc="CD ";
			} elsif (($woc=~ m/^SC1/)&&($wod eq "ILE")){$woc="CG1";
			} elsif (($woc=~ m/^SC1/)&&($wod eq "MET")){$woc="SD ";
			} elsif (($woc=~ m/^SC1/)&&($wod eq "VAL")){$woc="CB ";
			} elsif (($woc=~ m/^SC1/)&&($wod eq "SER")){$woc="OG ";
			} elsif (($woc=~ m/^SC1/)&&($wod eq "THR")){$woc="OG1";
		 	} elsif ($woc=~ m/^SC1/){$woc="CG ";
			} elsif (($woc=~ m/^SC2/)&&($wod eq "ARG")){$woc="NE ";
			} elsif (($woc=~ m/^SC2/)&&($wod eq "LYS")){$woc="CE ";
			} elsif (($wod eq "HIS")&&($woc=~ m/^SC2/)){$woc="ND1";
			} elsif (($wod eq "HIS")&&($woc=~ m/^SC3/)){$woc="NE2";
			} elsif (($wod eq "TYR")&&($woc=~ m/^SC2/)){$woc="CE1";
			} elsif (($wod eq "TYR")&&($woc=~ m/^SC3/)){$woc="CE2";
			} elsif (($wod eq "PHE")&&($woc=~ m/^SC2/)){$woc="CE1";
			} elsif (($wod eq "PHE")&&($woc=~ m/^SC3/)){$woc="CE2";
			} elsif (($wod eq "TRP")&&($woc=~ m/^SC2/)){$woc="CE2";
			} elsif (($wod eq "TRP")&&($woc=~ m/^SC3/)){$woc="CE3";
			} elsif (($wod eq "TRP")&&($woc=~ m/^SC4/)){$woc="CH2";

#### Older Martini ####
			} elsif (($woc=~ m/^B/)||($woc=~ m/^0B/)||($woc=~ m/^4B/)||($woc=~ m/^5B/)){$woc=" CA ";
		 	} elsif ($wod eq "CYS"){$woc=" SG ";
		 	} elsif (($wod eq "GLU")||($wod eq "GLN")) {$woc=" CD ";
			} elsif ($wod eq "ILE"){$woc=" CG1";
			} elsif ($wod eq "MET"){$woc=" SD ";
			} elsif ($wod eq "VAL"){$woc=" CB ";
			} elsif ($wod eq "SER"){$woc=" OG ";
			} elsif ($wod eq "THR"){$woc=" OG1";
			} elsif (($wod eq "ARG")&&($woc=~ m/Qd/)){$woc=" NE ";
			} elsif (($wod eq "LYS")&&($woc=~ m/Qd/)){$woc=" CE ";
			} elsif (($wod eq "HIS")&&($woc=~ m/^S/)&&($count==1)){$woc=" CG ";
				$count++;
			} elsif (($wod eq "HIS")&&($woc=~ m/^S/)&&($count==2)){$woc=" ND1";
				$count++;
			} elsif (($wod eq "HIS")&&($woc=~ m/^S/)&&($count==3)){$woc=" NE2";
				$count=1;
			} elsif (($wod eq "TYR")&&($woc=~ m/^S/)&&($count==1)){$woc=" CG ";
				$count++;
			} elsif (($wod eq "TYR")&&($woc=~ m/^S/)&&($count==2)){$woc=" CE1";
				$count++;
			} elsif (($wod eq "TYR")&&($woc=~ m/^S/)&&($count==3)){$woc=" CE2";
				$count=1;
			} elsif (($wod eq "PHE")&&($woc=~ m/^S/)&&($count==1)){$woc=" CG ";
				$count++;
			} elsif (($wod eq "PHE")&&($woc=~ m/^S/)&&($count==2)){$woc=" CE1";
				$count++;
			} elsif (($wod eq "PHE")&&($woc=~ m/^S/)&&($count==3)){$woc=" CE2";
				$count=1;
			} elsif (($wod eq "TRP")&&($woc=~ m/^S/)&&($count==1)){$woc=" CG ";
				$count++;
			} elsif (($wod eq "TRP")&&($woc=~ m/^S/)&&($count==2)){$woc=" CE2";
				$count++;
			} elsif (($wod eq "TRP")&&($woc=~ m/^S/)&&($count==3)){$woc=" CE3";
				$count++;
			} elsif (($wod eq "TRP")&&($woc=~ m/^S/)&&($count==4)){$woc=" CH2";
				$count=1;
				} else {$woc=" CG ";
			}
		print ATM $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
		} else {
		print ATM $_,"\n";
		}
	}
close (PRO);
close (ATM);
}
elsif ($forcefield eq "bond") { 

my @prot = <PRO>;
s/CB  ASP/CG  ASP/g for @prot;
s/CB  GLU/CD  GLU/g for @prot;
s/CB  PHE/CE1 PHE/g for @prot;
s/CG  PHE/CE2 PHE/g for @prot;
s/CB  HIS/CD2 HIS/g for @prot;
s/CG  HIS/ND1 HIS/g for @prot;
s/CB  ILE/CG1 ILE/g for @prot;
s/CB  LEU/CG  LEU/g for @prot;
s/CG  ARG/CZ  ARG/g for @prot;
s/CB  ARG/CG  ARG/g for @prot;
s/CB  SER/OG  SER/g for @prot;
s/CB  LYS/CG  LYS/g for @prot;
s/CG  LYS/NZ  LYS/g for @prot;
s/CB  TRP/CE3 TRP/g for @prot;
s/CG  TRP/NE1 TRP/g for @prot;
s/CB  CYS/SG  CYS/g for @prot;
s/CB  MET/SD  MET/g for @prot;
s/CB  GLN/CD  GLN/g for @prot;
s/CB  ASN/CG  ASN/g for @prot;
s/CB  THR/OG1 THR/g for @prot;
s/CB  VAL/CB  VAL/g for @prot;
s/CB  TYR/CE1 TYR/g for @prot;
s/CG  TYR/CE2 TYR/g for @prot;
s/CB  PRO/CG  PRO/g for @prot;
print ATM @prot;
close (ATM);
close (PRO);
}

open ATM, "<", 'atomistic-particles.pdb' or die;

$resnum = 0;

while (<ATM>){
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
	if (($woc eq CA)&&($woa eq ATOM)){
	$resnum++;
		}
	}
close (ATM);	

###################################################
#                                                 #
#                  ALIGN PROTEIN                  #
#                                                 #
################################################### 

if ($method eq "align"){

system("$gromacs/$gmx make_ndx -f cg2atalignfile.pdb -o nohydrogens.ndx << EOD
del 0-1
del 1-20
q
EOD
");

system("$gromacs/$gmx editconf -resnr 1 -f cg2atalignfile.pdb -o nohydrogens.pdb -n nohydrogens.ndx << EOD
0
EOD
");

system("$gromacs/$gmx editconf -resnr 1 -f nohydrogens.pdb -o labelZ.gro");
system("$gromacs/$gmx editconf -f labelZ.gro -o nolabel.pdb");

if ($protnum > 1){

system("rm chain_*.pdb");

@chainid = 1 .. $protnum;

@chain_end = map($_-1,@chain_num);
push (@chain_end,$resnum);

system("$gromacs/$gmx make_ndx -f atomistic-particles.pdb -o split-chain.ndx << EOD
del 1-30
r${chain_num[0]}-${chain_end[1]}&aCA
q
EOD
");

for($i = 1; $i <= $protnum ; $i++) {
system("$gromacs/$gmx make_ndx -f atomistic-particles.pdb -n split-chain.ndx  -o split-chain.ndx << EOD
r${chain_num[$i]}-${chain_end[$i+1]}&aCA
q
EOD");
}

for($i = 1; $i <= $protnum ; $i++) {
system("$gromacs/$gmx confrms -nice 0 -one -f1 atomistic-particles.pdb -f2 nolabel.pdb -name -o chain_${chainid[$i]}.pdb -n1 split-chain.ndx << EOD
$i
3
EOD");

unless (-e "chain_${chainid[$i]}.pdb"){

open (PYMOL, "> align.pml") || die "Cannot open align.pml\n";
print PYMOL "set retain_order,1\n";
print PYMOL "set pdb_retain_ids,1\n";
print PYMOL "load atomistic-particles.pdb\n";
print PYMOL "load nolabel.pdb\n";
print PYMOL "align nolabel and name CA and ${chainid[$i]}, atomistic-particles and ${chainid[$i]} and name CA\n";
print PYMOL "save chain_${chainid[$i]}.pdb, nolabel\n";
close (PYMOL);

system("$pymol -cq align.pml");

	}
}


system("cat chain_*.pdb | grep ATOM > protein_aligned.pdb");

system("$gromacs/$gmx editconf -resnr 1 -f protein_aligned.pdb -o protein-sorted.pdb");

open SID, "<", 'protein-sorted.pdb' or die;
open CHA, ">", 'aligned_chain.pdb' or die;

$oldres=-999;
$nchain=0;
$chain=();
@refchain = 'A' .. 'Z';

while (<SID>) {
           $len=length($_);
           $str1=(substr($_,0,21));
	   $str2=(substr($_,27,($len-26)));
	   chomp;
	   ($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
           $res=$woe;
	   if (($res==$chain_num[$nchain]) && ($res ne $oldres) && ($woa eq "ATOM")) {
             $oldres=$res;
             $chain=$refchain[$nchain] ;         
             $nchain++;
             }
           if (length($res)==4) {
             print CHA $str1,"",$chain,"",$res,"",$str2;
             } elsif (length($res)==3) {
             	print CHA $str1,"",$chain," ",$res,"",$str2;
             	} elsif (length($res)==2) {
                    print CHA $str1,"",$chain,"  ",$res,"",$str2;
            		} elsif (length($res)==1) {
             		print CHA $str1,"",$chain,"   ",$res,"",$str2;
           		 } else {print CHA $_,"\n"};
	   $oldres=$res;
           	}
close(CHA);
close(SID);

open OPT, ">", 'optimize.py' or die;
print OPT "\nfrom modeller import *";
print OPT "\nfrom modeller.scripts import complete_pdb";
print OPT "\nfrom modeller.optimizers import conjugate_gradients, molecular_dynamics, actions";
print OPT "\nenv = environ()";
print OPT "\nenv.io.atom_files_directory = ['../atom_files']";
print OPT "\nenv.edat.dynamic_sphere = True";
print OPT "\nenv.libs.topology.read(file='\$(LIB)/top_heav.lib')";
print OPT "\nenv.libs.parameters.read(file='\$(LIB)/par.lib')";
print OPT "\ncode = 'aligned_chain.pdb'";
print OPT "\nmdl = complete_pdb(env, code)";
print OPT "\nmdl.write(file=code+'.ini')";
print OPT "\natmsel = selection(mdl)";
print OPT "\nmdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)";
print OPT "\nmdl.restraints.write(file=code+'.rsr')";
print OPT "\nmpdf = atmsel.energy()";
print OPT "\ncg = conjugate_gradients(output='REPORT')";
print OPT "\nmd = molecular_dynamics(output='REPORT')";
print OPT "\ntrcfil = file(code+'.D00000001', 'w')";
print OPT "\ncg.optimize(atmsel, max_iterations=5)";
print OPT "\nmd.optimize(atmsel, temperature=300, max_iterations=50,";
print OPT "\nactions=[actions.write_structure(10, code+'.D9999%04d.pdb'),";
print OPT "\nactions.trace(10, trcfil)])";
print OPT "\ncg.optimize(atmsel, max_iterations=20,";
print OPT "\nactions=[actions.trace(5, trcfil)])";


print OPT "\nmpdf = atmsel.energy()";
print OPT "\nmdl.write(file='model.pdb')";
close(OPT);

if ($use_modeller = "yes"){
system("$modeller optimize.py");
}
if ($use_modeller = "no"){
system("$gromacs/$gmx editconf -f aligned_chain.pdb -o model.pdb");
}

}

if ($protnum == 1){

open SID, "<", 'nolabel.pdb' or die;
open CHA, ">", 'chainlabel.pdb' or die;

$oldres=-999;
$nchain=0;
$chain=();
@refchain='A' .. 'Z'; 

while (<SID>) {
           $len=length($_);
           $str1=(substr($_,0,21));
	   $str2=(substr($_,26,($len-26)));
	   chomp;
	   ($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
           $res=$woe;
	   if (($res==$chain_num[$nchain]) && ($res ne $oldres) && ($woa eq "ATOM")) {
             $oldres=$res;
             $chain=$refchain[$nchain];         
             $nchain++;
             }
           if (length($res)==4) {
             print CHA $str1,"",$chain,"",$res,"",$str2;
             } elsif (length($res)==3) {
             	print CHA $str1,"",$chain," ",$res,"",$str2;
             	} elsif (length($res)==2) {
                    print CHA $str1,"",$chain,"  ",$res,"",$str2;
            		} elsif (length($res)==1) {
             		print CHA $str1,"",$chain,"   ",$res,"",$str2;
           		 } ;
	   $oldres=$res;
           }
close(CHA);
close(SID);

if ($align eq "all"){
#system("$gromacs/$gmx make_ndx -f chainlabel.pdb -o cgtm.ndx");
#system("$gromacs/$gmx make_ndx -f atomistic-particles.pdb -o attm.ndx");
system("$gromacs/$gmx confrms -nice 0 -one -f1 atomistic-particles.pdb -f2 chainlabel.pdb -name -o model.pdb << EOD
3
3
EOD");

unless (-e "model.pdb"){

open (PYMOL, "> align.pml") || die "Cannot open align.pml\n";
print PYMOL "set retain_order,1\n";
print PYMOL "set pdb_retain_ids,1\n";
print PYMOL "load atomistic-particles.pdb\n";
print PYMOL "load chainlabel.pdb\n";
print PYMOL "align chainlabel and name CA, atomistic-particles and name CA\n";
print PYMOL "select OXT, name OXT or name O2\n";
print PYMOL "remove OXT\n";
print PYMOL "save model.pdb, chainlabel\n";
close (PYMOL);

system("$pymol -cq align.pml");

        }

}
elsif ($align eq "tm"){
system("$gromacs/$gmx make_ndx -f chainlabel.pdb -o cgtm.ndx < tm.list");
system("$gromacs/$gmx make_ndx -f atomistic-particles.pdb -o attm.ndx < tm.list");
system("$gromacs/$gmx confrms -nice 0 -one -f1 atomistic-particles.pdb -f2 chainlabel.pdb -name -o model.pdb -n1 attm.ndx -n2 cgtm.ndx<< EOD
3
3
EOD");

unless (-e "model.pdb"){

open (PYMOL, "> align.pml") || die "Cannot open align.pml\n";
print PYMOL "set retain_order,1\n";
print PYMOL "set pdb_retain_ids,1\n";
print PYMOL "load atomistic-particles.pdb\n";
print PYMOL "load chainlabel.pdb\n";
print PYMOL "align chainlabel and name CA, atomistic-particles and name CA\n";
print PYMOL "select OXT, name OXT or name O2\n";
print PYMOL "remove OXT\n";
print PYMOL "save model.pdb, chainlabel\n";
close (PYMOL);

system("$pymol -cq align.pml");

        }

}

	}
}

###################################################
#                                                 #
#            FULL PROTEIN CONVERSION              #
#                                                 #
################################################### 

elsif ($method eq "full"){

`rm -rf Res_CG fit_*`;
mkdir('Res_CG');

$loop = 1;

while ($loop <= $resnum){

open (ATM, "atomistic-particles.pdb") || die "Cannot Open\n";
open (RESCG, "> Res_CG/res_cg_${loop}.pdb") || die "No output possible in Res_CG\n";

while (<ATM>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe eq $loop){
		if ((length($woe))==4){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"   1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
			}
		elsif ((length($woe))==3){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"  1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==2){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe," 1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		elsif ((length($woe))==1){
			print RESCG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,"1",$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n";
			}
		}
	}			
close(ATM);

print RESCG "ATOM      3  DA  PRO     1       1.000   1.000   1.000  1.00  0.00","\n" ;

close(RESCG);

$loop++;
}

$loop_main = 1;

while ($loop_main <= $resnum){

mkdir("fit_${loop_main}");

open (RES, "< Res_CG/res_cg_${loop_main}.pdb") || die;
open (CAF, "> fit_${loop_main}/residue_${loop_main}.pdb") || die "No output possible\n";

while (<RES>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($woc eq CA)&&(($wod eq "ALA")||($wod eq "CYS")||($wod eq "GLU")||($wod eq "GLN")||($wod eq "ASN")||($wod eq "ARG")||($wod eq "VAL")||($wod eq "ASP")||($wod eq "ILE")||($wod eq "LEU")||($wod eq "MET")||($wod eq "GLY")||($wod eq "LYS")||($wod eq "SER")||($wod eq "THR"))){
	`cp Res_CG/res_cg_${loop_main}.pdb fit_${loop_main}/fit_${loop_main}.pdb`;
	}	
	elsif (($woc eq CA)&&($wod eq "PRO")) {
	system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 Sidechains/PRO.pdb -name -o fit_${loop_main}/fit_${loop_main}.pdb << EOD
	1
	1
	EOD");
	}
	elsif (($woc eq CA)&&($wod eq "PHE")) {
	system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 Sidechains/PHE.pdb -name -o fit_${loop_main}/fit_${loop_main}.pdb << EOD
	1
	1
	EOD");
	}
	elsif (($woc eq CA)&&($wod eq "HIS")) {
	system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 Sidechains/HIS.pdb -name -o fit_${loop_main}/fit_${loop_main}.pdb << EOD
	1
	1
	EOD");
	}
	elsif (($woc eq CA)&&($wod eq "TRP")) {
	system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 Sidechains/TRP.pdb -name -o fit_${loop_main}/fit_${loop_main}.pdb << EOD
	1	
	1
	EOD");
	}
	elsif (($woc eq CA)&&($wod eq "TYR")) {
	system("$gromacs/$gmx confrms -nice 0 -one -f1 Res_CG/res_cg_${loop_main}.pdb -f2 Sidechains/TYR.pdb -name -o fit_${loop_main}/fit_${loop_main}.pdb << EOD
	1
	1
	EOD");
	}
	if ($woc eq "CA") {
	print CAF $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
}

open (FIT, "< fit_${loop_main}/fit_${loop_main}.pdb") || die;

while (<FIT>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if (($woa eq ATOM)&&(($woc ne CA)&&($woc ne DA))) {
	print CAF $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
}
close (FIT);
close (CAF);
close (RES);

open (RSL, "< fit_${loop_main}/residue_${loop_main}.pdb") || die;
open RSN, ">>", 'protein_sidechains.pdb' or die;

while (<RSL>) {
           $len=length($_);
           $str1=(substr($_,0,22));
	   $str2=(substr($_,26,($len-26)));
	   chomp;
	   ($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
	   if ($woa eq "ATOM") {
            if (length($loop_main)==4) {
             print RSN $str1,"",$loop_main,"",$str2;
             } elsif (length($loop_main)==3) {
             	print RSN $str1," ",$loop_main,"",$str2;
             	} elsif (length($loop_main)==2) {
                    print RSN $str1,"  ",$loop_main,"",$str2;
            		} elsif (length($loop_main)==1) {
             		print RSN $str1,"   ",$loop_main,"",$str2;
           		 } else {print RSN $_,"\n"};
           	}
	}
close(RSL);
close(RSN);

$loop_main++;
}

open SID, "<", 'protein_sidechains.pdb' or die;
open CHA, ">", 'protein_chain.pdb' or die;

$oldres=-999;
$nchain=0;
$chain=();
@refchain='A' .. 'Z'; 
  
while (<SID>) {
           $len=length($_);
           $str1=(substr($_,0,21));
	   $str2=(substr($_,27,($len-26)));
	   chomp;
	   ($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
           $res=$woe;
	   if (($res==$chain_num[$nchain]) && ($res ne $oldres) && ($woa eq "ATOM")) {
             $oldres=$res;
             $chain=$refchain[$nchain];         
             $nchain++;
             }
           if (length($res)==4) {
             print CHA $str1,"",$chain,"",$res,"",$str2;
             } elsif (length($res)==3) {
             	print CHA $str1,"",$chain," ",$res,"",$str2;
             	} elsif (length($res)==2) {
                    print CHA $str1,"",$chain,"  ",$res,"",$str2;
            		} elsif (length($res)==1) {
             		print CHA $str1,"",$chain,"   ",$res,"",$str2;
           		 } else {print CHA $_,"\n"};
	   $oldres=$res;
           }
close(CHA);
close(SID);

open (PCH, "protein_chain.pdb") || die "Cannot open protein_chain.pdb\n";
open CHA, ">", 'chainA.pdb' or die;
open CHB, ">", 'chainB.pdb' or die;
open CHC, ">", 'chainC.pdb' or die;
open CHD, ">", 'chainD.pdb' or die;
open CHE, ">", 'chainE.pdb' or die;
open CHF, ">", 'chainF.pdb' or die;
open CHG, ">", 'chainG.pdb' or die;
open CHH, ">", 'chainH.pdb' or die;
open CHI, ">", 'chainI.pdb' or die;
open CHJ, ">", 'chainJ.pdb' or die;
open CHK, ">", 'chainK.pdb' or die;
open CHL, ">", 'chainL.pdb' or die;
open CHM, ">", 'chainM.pdb' or die;
open CHN, ">", 'chainN.pdb' or die;
open CHO, ">", 'chainO.pdb' or die;
open CHP, ">", 'chainP.pdb' or die;

while (<PCH>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woe=~ m/^A/){
	print CHA $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^B/){
	print CHB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^C/){
	print CHC $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^D/){
	print CHD $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^E/){
	print CHE $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^F/){
	print CHF $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^G/){
	print CHG $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^H/){
	print CHH $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^I/){
	print CHI $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^J/){
	print CHJ $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^K/){
	print CHK $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^L/){
	print CHL $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^M/){
	print CHM $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^N/){
	print CHN $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^O/){
	print CHO $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
	if ($woe=~ m/^P/){
	print CHP $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" ;
	}
}
close (PCH);
close (CHA);
close (CHB);
close (CHC);
close (CHD);
close (CHE);
close (CHF);
close (CHG);
close (CHH);
close (CHI);
close (CHJ);
close (CHK);
close (CHL);
close (CHM);
close (CHN);
close (CHO);
close (CHP);

system("$pulchra -s chainA.pdb");
system("$pulchra -s chainB.pdb");
system("$pulchra -s chainC.pdb");
system("$pulchra -s chainD.pdb");
system("$pulchra -s chainE.pdb");
system("$pulchra -s chainF.pdb");
system("$pulchra -s chainG.pdb");
system("$pulchra -s chainH.pdb");
system("$pulchra -s chainI.pdb");
system("$pulchra -s chainJ.pdb");
system("$pulchra -s chainK.pdb");
system("$pulchra -s chainL.pdb");
system("$pulchra -s chainM.pdb");
system("$pulchra -s chainN.pdb");
system("$pulchra -s chainO.pdb");
system("$pulchra -s chainP.pdb");

system("cat chain*.rebuilt.pdb | grep ATOM > protein_pulchra.pdb");

open SID, "<", 'protein_pulchra.pdb' or die;
open CHA, ">", 'pulchra_chain.pdb' or die;

$oldres=-999;
$nchain=0;
$chain=();
@refchain='A' .. 'Z';

while (<SID>) {
           $len=length($_);
           $str1=(substr($_,0,21));
	   $str2=(substr($_,27,($len-26)));
	   chomp;
	   ($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
           $res=$woe;
	   if (($res==$chain_num[$nchain]) && ($res ne $oldres) && ($woa eq "ATOM")) {
             $oldres=$res;
             $chain=$refchain[$nchain];         
             $nchain++;
             }
           if (length($res)==4) {
             print CHA $str1,"",$chain,"",$res,"",$str2;
             } elsif (length($res)==3) {
             	print CHA $str1,"",$chain," ",$res,"",$str2;
             	} elsif (length($res)==2) {
                    print CHA $str1,"",$chain,"  ",$res,"",$str2;
            		} elsif (length($res)==1) {
             		print CHA $str1,"",$chain,"   ",$res,"",$str2;
           		 } else {print CHA $_,"\n"};
	   $oldres=$res;
           	}
close(CHA);
close(SID);

open OPT, ">", 'optimize.py' or die;
print OPT "\nfrom modeller import *";
print OPT "\nfrom modeller.scripts import complete_pdb";
print OPT "\nfrom modeller.optimizers import conjugate_gradients, molecular_dynamics, actions";
print OPT "\nenv = environ()";
print OPT "\nenv.io.atom_files_directory = ['../atom_files']";
print OPT "\nenv.edat.dynamic_sphere = True";
print OPT "\nenv.libs.topology.read(file='\$(LIB)/top_heav.lib')";
print OPT "\nenv.libs.parameters.read(file='\$(LIB)/par.lib')";
print OPT "\ncode = 'pulchra_chain.pdb'";
print OPT "\nmdl = complete_pdb(env, code)";
print OPT "\nmdl.write(file=code+'.ini')";
print OPT "\natmsel = selection(mdl)";
print OPT "\nmdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)";
print OPT "\nmdl.restraints.write(file=code+'.rsr')";
print OPT "\nmpdf = atmsel.energy()";
print OPT "\ncg = conjugate_gradients(output='REPORT')";
print OPT "\ntrcfil = file(code+'.D00000001', 'w')";
print OPT "\ncg.optimize(atmsel, max_iterations=5)";
print OPT "\nmpdf = atmsel.energy()";
print OPT "\nmdl.write(file='model.pdb')";
close(OPT);

system("$modeller optimize.py");
}

system("cat box.pdb model.pdb > protein_box.pdb");

if ($merge eq "yes"){
if ($ff eq "GROMOS96-43a1"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos43a1 -water spc -merge all");
} elsif ($ff eq "GROMOS96-43a2"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos43a2 -water spc -merge all");
} elsif ($ff eq "GROMOS96-45a3"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos45a3 -water spc -merge all");
} elsif ($ff eq "GROMOS96-53a5"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos53a5 -water spc -merge all");
} elsif ($ff eq "GROMOS96-53a6"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos53a6 -water spc -merge all");
} elsif ($ff =~ m/^OPLS/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff oplsaa -water tip4p -merge all");
} elsif ($ff =~ m/^CHARMM27/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff charmm27 -water tip3p -merge all");
} elsif ($ff =~ m/^CHARMM36/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff charmm36 -water tip3p -merge all");
} elsif ($ff =~ m/^AMBER/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff amber99 -water spce -merge all");
	}
}elsif ($merge eq "no"){
if ($ff eq "GROMOS96-43a1"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos43a1 -water spc");
} elsif ($ff eq "GROMOS96-43a2"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos43a2 -water spc");
} elsif ($ff eq "GROMOS96-45a3"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos45a3 -water spc");
} elsif ($ff eq "GROMOS96-53a5"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos53a5 -water spc");
} elsif ($ff eq "GROMOS96-53a6"){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff gromos53a6 -water spc");
} elsif ($ff =~ m/^OPLS/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff oplsaa -water tip4p");
} elsif ($ff =~ m/^CHARMM27/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff charmm27 -water tip3p");
} elsif ($ff =~ m/^CHARMM36/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff charmm36 -water tip3p");
} elsif ($ff =~ m/^AMBER/){
system("$gromacs/$gmx pdb2gmx -f protein_box.pdb -o p2g.pdb -p protein.top -ff amber99 -water spce");
	}
}

system("$gromacs/$gmx grompp -maxwarn 10 -f em.mdp -o prot_em -c p2g.pdb -p protein.top");
system("$gromacs/$gmx mdrun -deffnm prot_em");
system("$gromacs/$gmx editconf -resnr 1 -f prot_em.gro -o atomistic-protein.pdb");

open (SID, "atomistic-protein.pdb")  || die "Energy minimisation of AT protein failed\n";

$oldres=-999;
$nchain=0;
$chain=();
@refchain='A' .. 'Z'; 
  
open CHA, ">", 'atomistic-protein-chain.pdb' or die;

while (<SID>) {
           $len=length($_);
           $str1=(substr($_,0,21));
	   $str2=(substr($_,27,($len-26)));
	   chomp;
	   ($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
           $res=$woe;
	   if (($res==$chain_num[$nchain]) && ($res ne $oldres) && ($woa eq "ATOM")) {
             $oldres=$res;
             $chain=$refchain[$nchain];         
             $nchain++;
             }
           if (length($res)==4) {
             print CHA $str1,"",$chain,"",$res,"",$str2;
             } elsif (length($res)==3) {
             	print CHA $str1,"",$chain," ",$res,"",$str2;
             	} elsif (length($res)==2) {
                    print CHA $str1,"",$chain,"  ",$res,"",$str2;
            		} elsif (length($res)==1) {
             		print CHA $str1,"",$chain,"   ",$res,"",$str2;
           		 } else {print CHA $_,"\n"};
	   $oldres=$res;
           	}
close(CHA);
close(SID);

open CHA, "<", 'atomistic-protein-chain.pdb' or die;
close(CHA);


###################################################  
#                                                 #
#           COMBINE PROTEIN WITH LIPID            #
#                                                 #
###################################################

###################################################
#                                                 #
#              LIPID DELETION START               #
#                                                 #
###################################################

open (PRO, "protein.top")  || die "Something went wrong\n";
open(TOP, "> topol.top") || die "No Output Possible\n";
my @pro = <PRO>;
if ($ff eq "GROMOS96-43a1"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-43a2"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-45a3"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-53a5"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-53a6"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-gmx53a6.itp\"\n/g for @pro;
} elsif ($ff eq "OPLSAA"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/oplsaa-lipids.itp\"\n\#include \"itp\/lipid-oplsaa.itp\"\n/g for @pro;
}elsif ($ff eq "OPLSUA"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/oplsua-lipids.itp\"\n\#include \"itp\/lipid-oplsua.itp\"\n/g for @pro;
}elsif ($ff eq "CHARMM27"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/charmm27-lipids.itp\"\n/g for @pro;
}elsif ($ff eq "CHARMM36"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/charmm36-lipids.itp\"\n/g for @pro;
}elsif ($ff eq "AMBER"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/lipid-amber.itp\"\n/g for @pro;
}
print TOP @pro;
close(PRO);
if ($lipidnum > 0) {

###################################################
#                                                 #
#                ALIGNMENT METHOD                 #
#                                                 #
###################################################

if ($method eq "align"){
if ($deletion eq "gmxsolvate"){
system("$gromacs/$gmx solvate -cp atomistic-protein.pdb -cs atomistic-lipid.pdb -o prot+lipid+box.pdb");
}


elsif (($deletion eq "g_membed")||($deletion eq "alchembed")){
if ($deletion eq "alchembed"){

if ($merge eq "no"){

if ($ff eq "GROMOS96-43a1"){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff gromos43a1  -water spc -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff gromos43a1 -water spc -merge all");
} elsif ($ff eq "GROMOS96-43a2"){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff gromos43a2 -water spc -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff gromos43a2 -water spc -merge all");
} elsif ($ff eq "GROMOS96-45a3"){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff gromos53a3 -water spc -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff gromos45a3 -water spc -merge all");
} elsif ($ff eq "GROMOS96-53a5"){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff gromos53a5 -water spc -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff gromos53a5 -water spc -merge all");
} elsif ($ff eq "GROMOS96-53a6"){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff gromos53a6 -water spc -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff gromos53a6 -water spc -merge all");
} elsif ($ff =~ m/^OPLS/){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff oplsaa -water tip4p -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff oplsaa -water tip4p -merge all");
} elsif ($ff =~ m/^CHARMM27/){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff charmm27 -water tip3p -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff charmm27 -water tip3p -merge all");
} elsif ($ff =~ m/^CHARMM36/){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff charmm36 -water tip3p -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff charmm36 -water tip3p -merge all");
} elsif ($ff =~ m/^AMBER/){
system("$gromacs/$gmx pdb2gmx -f atomistic-protein-chain.pdb -o chains.pdb -chainsep id -ignh -ff amber99 -water spce -p alchembed.top");
system("$gromacs/$gmx pdb2gmx -f chains.pdb -o protein.pdb -p protein.top -i posre.itp -ignh -ff amber99 -water spce -merge all");
	}
system("cat protein.pdb atomistic-lipid.pdb | grep ATOM > prot+lipid.pdb");
}

if ($merge eq "yes"){
system("cat atomistic-protein.pdb atomistic-lipid.pdb | grep ATOM > prot+lipid.pdb");
}

system("cat box.pdb prot+lipid.pdb > prot-lipid-box.pdb");

open (PRO, "protein.top")  || die "Something went wrong\n";
open(ALC, "> alc-membed.top") || die "No Output Possible\n";
my @pro = <PRO>;
if ($ff eq "GROMOS96-43a1"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-43a2"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-45a3"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-53a5"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-berger.itp\"\n/g for @pro;
} elsif ($ff eq "GROMOS96-53a6"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/gromos-lipids.itp\"\n\#include \"itp\/lipid-gmx53a6.itp\"\n/g for @pro;
} elsif ($ff eq "OPLSAA"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/oplsaa-lipids.itp\"\n\#include \"itp\/lipid-oplsaa.itp\"\n/g for @pro;
}elsif ($ff eq "OPLSUA"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/oplsua-lipids.itp\"\n\#include \"itp\/lipid-oplsua.itp\"\n/g for @pro;
}elsif ($ff eq "CHARMM27"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/charmm27-lipids.itp\"\n/g for @pro;
}elsif ($ff eq "CHARMM36"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/charmm36-lipids.itp\"\n/g for @pro;
}elsif ($ff eq "AMBER"){
s/forcefield.itp\"/forcefield.itp\"\n\#include \"itp\/lipid-amber.itp\"\n/g for @pro;
}
print ALC @pro;
close ALC;
}

elsif ($deletion eq "g_membed"){
system("cp topol.top alc-membed.top");
}

$popc = 0;
$pope = 0;
$pvpg = 0;
$pvpe = 0;
$popg = 0;
$pops = 0;
$popa = 0;
$dopc = 0;
$dppc = 0;
$dppe = 0;
$dppg = 0;
$ppcs = 0;
$dspc = 0;
$dspe = 0;
$dspg = 0;
$dmpc = 0;
$dmpe = 0;
$dmpg = 0;
$dhpc = 0;
$pip3 = 0;
$pip2 = 0;
$bog = 0;
$dpc = 0;
$chol = 0;
$ddm = 0;
$lmpg = 0;
$card = 0;
$lmg = 0;
$sqd = 0;
$dgd = 0;
$mag = 0;
$dag = 0;
$sds = 0;

open (INP, "atomistic-lipid.pdb") || die "Cannot open atomistic-lipid.pdb\n";
$atom = 0;
while (<INP>){
	chomp;
	$woc=substr($_,13,3);
	$wod=substr($_,17,3);
	$woc=~ s/^\s+|\s+$//g;	
	$wod=~ s/^\s+|\s+$//g;
	
        if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq C1)&&($wod eq POP))){                  
        $popc++;          
        }
        if (($ff =~ m/^CHARMM/)&&(($woc eq NC3)&&($wod eq POP))){   
        $popc++;                
        }
        if (($ff =~ m/^CHARMM/)&&(($woc eq NH3)&&($wod eq POP))){   
        $pope++;                
        }
	if (($ff =~ m/^CHARMM/)&&(($woc eq C1H)&&($wod eq POP))){   
        $popg++;                
        }
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq H1)&&($wod eq POP))){
	$pope++;
	}
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq H0)&&($wod eq POP))){
	$popg++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq H1)&&($wod eq PVP))){
	$pvpe++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq H0)&&($wod eq PVP))){
	$pvpg++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq N5)&&($wod eq POP))){
	$pops++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq P1)&&($wod eq POP))){
	$popa++;
	}
	if (($ff =~ m/^CHARMM/)&&(($woc eq NCO)&&($wod eq POP))){
        $pops++;
        }
	if ((($woc eq N1)&&($ff eq AMBER))&&($wod eq DOP)){
	$dopc++;
	}	
        if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq C1)&&($wod eq DOP))){
	$dopc++;
	}
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq C1)&&($wod eq DPP))){
	$dppc++;
	}
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq H1)&&($wod eq DPP))){
	$dppe++;
	}		
	if (($ff =~ m/^GROMOS/)&&($woc eq H0)&&($wod eq DPP)){
	$dppg++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq N14)&&($wod eq PPC)){
	$ppcs++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq C1)&&($wod eq DSP)){
	$dspc++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq H1)&&($wod eq DSP)){
	$dspe++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq H0)&&($wod eq DSP)){
	$dspg++;
	}
	if (($woc eq C1)&&($wod eq DMP)){
	$dmpc++;
	}
	if (($woc eq H1)&&($wod eq DMP)){
	$dmpe++;
	}
	if (($woc eq H0)&&($wod eq DMP)){
	$dmpg++;
	}
	if (($woc eq C1)&&($wod eq DHP)){
	$dhpc++;
	}
	if ((($woc eq HBE)||($woc eq P5))&&($wod eq PIP)){
	$pip2++;
	}
	if (($woc eq PBE)&&($wod eq PIP)){
	$pip3++;
	}
	if (($woc eq C6)&&($wod eq BOG)){
	$bog++;
	}	
	if (($woc eq C10)&&($wod eq DDM)){
	$ddm++;
	}
	if (($woc eq P9)&&($wod eq LMP)){
	$lmpg++;
	}	
	if (($woc eq N1)&&($wod eq DPC)){
	$dpc++;
	}
        if (($woc eq N4)&&($wod eq DPC)){
        $dpc++;
        }
	if (($woc eq C10)&&($wod eq CHO)){
	$chol++;
	}
	if (($woc eq P8)&&($wod eq CAR)){
	$card++;
	}
	if (($woc eq O1)&&($wod eq LMG)){
	$lmg++;
	}
	if (($woc eq O1)&&($wod eq SQD)){
	$sqd++;
	}
	if (($woc eq O1)&&($wod eq DGD)){
	$dgd++;
	}	
	if (($woc eq H1A)&&($wod eq DAG)){
	$dag++;
	}
	if (($woc eq H1A)&&($wod eq MAG)){
	$mag++;
	}
	if (($woc eq SAO)&&($wod eq SDS)){
	$sds++;
	}
	if (($ff =~ m/^AMBER/)&&($woc eq O1)&&($wod eq DPP)){
	$dppg++;
	}
}
close (INP);

open(OUT, ">> alc-membed.top") || die "No Output Possible\n";

unless ($popc == 0){ 
print OUT "\nPOPC $popc";
}
unless ($pope == 0){
print OUT "\nPOPE $pope";
}
unless ($popg == 0){ 
print OUT "\nPOPG $popg";
}
unless ($pops == 0){
print OUT "\nPOPS $pops";
}
unless ($popa == 0){ 
print OUT "\nPOPA $popa";
}
unless ($dopc == 0){
print OUT "\nDOPC $dopc";
}
unless ($dppc == 0){ 
print OUT "\nDPPC $dppc";
}
unless ($dppe == 0){ 
print OUT "\nDPPE $dppe";
}
unless ($dppg == 0){ 
print OUT "\nDPPG $dppg";
}
unless ($ppcs == 0){ 
print OUT "\nPPCS $ppcs";
}
unless ($dspc == 0){
print OUT "\nDSPC $dspc";
}
unless ($dspe == 0){
print OUT "\nDSPE $dspe";
}
unless ($dspg == 0){
print OUT "\nDSPG $dspg";
}
unless ($dmpc == 0){
print OUT "\nDMPC $dmpc";
}
unless ($dmpe == 0){
print OUT "\nDMPE $dmpe";
}
unless ($dmpg == 0){
print OUT "\nDMPG $dmpg";
}
unless ($pip2 == 0){
print OUT "\nPIP2 $pip2";
}
unless ($pip3 == 0){
print OUT "\nPIP3 $pip3";
}
unless ($dhpc == 0){
print OUT "\nDHPC $dhpc";
}
unless ($bog == 0){ 
print OUT "\nBOG $bog\n";
}
unless ($ddm == 0){ 
print OUT "\nDDM $ddm\n";
}
unless ($lmpg == 0){ 
print OUT "\nLMPG $lmpg\n";
}
unless ($dpc == 0){ 
print OUT "\nDPC $dpc\n";
}
unless ($chol == 0){ 
print OUT "\nCHOL $chol\n";
}
unless ($card == 0){ 
print OUT "\nCARD $card\n";
}
unless ($pvpe == 0){
print OUT "\nPVPE $pvpe";
}
unless ($pvpg == 0){ 
print OUT "\nPVPG $pvpg";
}
unless ($lmg == 0){ 
print OUT "\nLMG $lmg\n";
}
unless ($sqd == 0){ 
print OUT "\nSQD $sqd";
}
unless ($dgd == 0){
print OUT "\nDGD $dgd";
}
unless ($dag == 0){
print OUT "\nDAG $dag";
}
unless ($mag == 0){
print OUT "\nMAG $mag";
}
unless ($sds == 0){
print OUT "\nSDS $sds";
}
close (OUT);

system("$gromacs/$gmx make_ndx -f prot-lipid-box.pdb -o alc-membed.ndx  << EOD
del 2-30
!1|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
name 2 Lipid
q
EOD");
}

if ($deletion eq "g_membed"){
system("cat atomistic-protein.pdb atomistic-lipid.pdb | grep ATOM > prot+lipid.pdb");
system("cat box.pdb prot+lipid.pdb > prot-lipid-box.pdb");

open GMDP, ">", 'membed.mdp' or die;
print GMDP "\nintegrator      = md";
print GMDP "\ntinit           = 0.0";
print GMDP "\ndt              = 0.002";
print GMDP "\nnsteps          = 1000";
print GMDP "\nenergygrps      = Protein";
print GMDP "\nfreezegrps      = Protein";
print GMDP "\nfreezedim       = Y Y Y";
print GMDP "\nenergygrp_excl  = Protein Protein";
close(GMDP);

open GMDP, ">", 'membed.dat' or die;
print GMDP "\nnxy = 1000\nnz = 0\nxyinit = 0.150000\nxyend = 1.000000\nzinit = 1.000000\nzend = 1.000000\nrad = 0.220000\nndiff = 0\nmaxwarn = 0\npieces = 1\nasymmetry = no\n";
close(GMDP);

system("$gromacs/$gmx grompp -f membed.mdp -p alc-membed.top -c prot-lipid-box.pdb -o membed.tpr -maxwarn 5");
system("$gromacs/$gmx mdrun -membed membed.dat -s membed.tpr -p alc-membed.top <<EOD
1
2
EOD");
}

# Renkinjutsu o suru #
if ($deletion eq "alchembed"){

open ALC, ">", 'alc.mdp' || die $!;
$alc_mdp = <<END;
integrator      = md
define		    = -DPOSRES
tinit           = 0.0
dt              = 0.001
nsteps          = 1000
nstxout         = 0
nstvout         = 0
nstfout         = 0
nstlog          = 50
nstenergy       = 50
nstxout-compressed       = 50
coulombtype     = PME
rlist                = 1.35
rcoulomb             = 1.2
vdw-type        = cutoff
Tcoupl          = Berendsen 
tc-grps         = SYSTEM
tau-t           = 1
ref-t           = 323.
gen_vel         = yes  
gen-temp        = 323
gen-seed        = -1
rvdw_switch          = 0.9
rvdw                 = 1.2
cutoff-scheme        = verlet
coulomb-modifier     = potential-shift-verlet
vdw-modifier         = potential-shift-verlet
verlet-buffer-tolerance  = 0.005
constraints     = none
constraint-algorithm = Lincs
free_energy	= yes
init_lambda	= 0.00
delta_lambda	= 1e-3
sc-alpha	= 0.1
sc-power	= 1.0
sc-r-power	= 6.0
couple-moltype	= Protein_chain_A
couple-lambda0	= none
couple-lambda1	= vdw
couple-intramol = yes
END
print ALC $alc_mdp;
close(ALC);

system("$gromacs/$gmx grompp -maxwarn 5 -f alc.mdp -p alc-membed.top -c prot-lipid-box.pdb -n alc-membed.ndx -o alc.tpr");
system("$gromacs/$gmx mdrun -v -deffnm alc -c prot+lipid+alc.pdb");
system("$gromacs/$gmx trjconv -f prot+lipid+alc.pdb -n alc-membed.ndx -o atomistic-lipid-hole.pdb -s alc.tpr -pbc res  << EOD
2
EOD");

open (ALC, "< atomistic-lipid-hole.pdb") || die "Alchembed Failed\n";
close ALC;

system("cat atomistic-protein.pdb atomistic-lipid-hole.pdb | grep ATOM > prot-lipid.pdb");
system("cat box.pdb prot-lipid.pdb > prot+lipid+box.pdb");

}

elsif ($deletion eq "pymol"){

system("cat atomistic-protein.pdb atomistic-lipid.pdb | grep ATOM > prot-lipid.pdb");

open (PYMOL, "> lipid.pml") || die "Cannot open lipid.pml\n";
print PYMOL "set retain_order,1\n";
print PYMOL "set pdb_retain_ids,1\n";
print PYMOL "load prot-lipid.pdb\n";
print PYMOL "select deletion, bm. polymer\n";
print PYMOL "remove br. deletion\n";
print PYMOL "select complete-delete, polymer\n";
print PYMOL "remove complete-delete\n";
print PYMOL "save atomistic-lipids-hole.pdb\n";
close (PYMOL);

system("$pymol -c lipid.pml");

system("cat atomistic-protein.pdb atomistic-lipids-hole.pdb | grep ATOM > prot+lipid.pdb");
system("cat box.pdb prot+lipid.pdb > prot+lipid+box.pdb");
	}

elsif ($deletion eq "none"){
system("cat atomistic-protein.pdb atomistic-lipid.pdb | grep ATOM > prot+lipid.pdb");
system("cat box.pdb prot+lipid.pdb > prot+lipid+box.pdb");
	}

}

###################################################
#                                                 #
#              ALIGNMENT METHOD END               #
#                                                 #
###################################################

elsif ($method eq "full"){
system("cat atomistic-protein.pdb atomistic-lipid.pdb | grep ATOM > prot+lipid.pdb");
system("cat box.pdb prot+lipid.pdb > prot+lipid+box.pdb");
}

###################################################
#                                                 #
#            LIPID DELETION COMPLETE              #
#                                                 #
###################################################

###################################################
#                                                 #
#        ENERGY MINIMISE PROTEIN & LIPID          #
#                                                 #
###################################################

$popc = 0;
$pope = 0;
$popg = 0;
$pvpg = 0;
$pvpe = 0;
$pops = 0;
$popa = 0;
$dopc = 0;
$dppc = 0;
$dppe = 0;
$dppg = 0;
$ppcs = 0;
$dspc = 0;
$dspe = 0;
$dspg = 0;
$dmpc = 0;
$dmpe = 0;
$dmpg = 0;
$dhpc = 0;
$pip3 = 0;
$pip2 = 0;
$bog = 0;
$chol = 0;
$ddm = 0;
$lmpg = 0;
$dpc = 0;
$card = 0;
$lmg = 0;
$sqd = 0;
$dgd = 0;
$mag = 0;
$dag = 0;
$sds = 0;

open (INP, "prot+lipid+box.pdb") || die "Cannot open prot+lipid+box.pdb\n";
$atom = 0;
while (<INP>){
	chomp;
	$woc=substr($_,13,3);
	$wod=substr($_,17,3);
	$woc=~ s/^\s+|\s+$//g;	
	$wod=~ s/^\s+|\s+$//g;
        if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq C1)&&($wod eq POP))){                  
        $popc++;          
        }
        if (($ff =~ m/^CHARMM/)&&(($woc eq NC3)&&($wod eq POP))){   
        $popc++;                
        }
        if (($ff =~ m/^CHARMM/)&&(($woc eq NH3)&&($wod eq POP))){   
        $pope++;                
        }
	if (($ff =~ m/^CHARMM/)&&(($woc eq C1H)&&($wod eq POP))){   
        $popg++;                
        }
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq H1)&&($wod eq POP))){
	$pope++;
	}
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq H0)&&($wod eq POP))){
	$popg++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq H1)&&($wod eq PVP))){
	$pvpe++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq H0)&&($wod eq PVP))){
	$pvpg++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq P1)&&($wod eq POP))){
	$popa++;
	}
	if (($ff =~ m/^GROMOS/)&&(($woc eq N5)&&($wod eq POP))){
	$pops++;
	}
	if (($ff =~ m/^CHARMM/)&&(($woc eq NCO)&&($wod eq POP))){
        $pops++;
        }
	if ((($woc eq N1)&&($ff eq AMBER))&&($wod eq DOP)){
	$dopc++;
	}	
        if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq C1)&&($wod eq DOP))){
	$dopc++;
	}
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq C1)&&($wod eq DPP))){
	$dppc++;
	}
	if (($ff =~ m/^GROMOS/||$ff =~ m/^OPLS/)&&(($woc eq H1)&&($wod eq DPP))){
	$dppe++;
	}		
	if (($ff =~ m/^GROMOS/)&&($woc eq H0)&&($wod eq DPP)){
	$dppg++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq N14)&&($wod eq PPC)){
	$ppcs++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq C1)&&($wod eq DSP)){
	$dspc++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq H1)&&($wod eq DSP)){
	$dspe++;
	}
	if (($ff =~ m/^GROMOS/)&&($woc eq H0)&&($wod eq DSP)){
	$dspg++;
	}
	if (($woc eq C1)&&($wod eq DMP)){
	$dmpc++;
	}
	if (($woc eq H1)&&($wod eq DMP)){
	$dmpe++;
	}
	if (($woc eq H0)&&($wod eq DMP)){
	$dmpg++;
	}
	if ((($woc eq HBE)||($woc eq P5))&&($wod eq PIP)){
	$pip2++;
	}
	if (($woc eq PBE)&&($wod eq PIP)){
	$pip3++;
	}
	if (($woc eq C6)&&($wod eq BOG)){
	$bog++;
	}	
	if (($woc eq C10)&&($wod eq DDM)){
	$ddm++;
	}
	if (($woc eq P9)&&($wod eq LMP)){
	$lmpg++;
	}
	if (($woc eq N1)&&($wod eq DPC)){
	$dpc++;
	}
        if (($woc eq N4)&&($wod eq DPC)){
        $dpc++;
        }
	if (($woc eq C10)&&($wod eq CHO)){
	$chol++;
	}
	if (($woc eq P8)&&($wod eq CAR)){
	$card++;
	}
	if (($woc eq O1)&&($wod eq LMG)){
	$lmg++;
	}
	if (($woc eq O1)&&($wod eq SQD)){
	$sqd++;
	}
	if (($woc eq O1)&&($wod eq DGD)){
	$dgd++;
	}	
	if (($woc eq H1A)&&($wod eq DAG)){
	$dag++;
	}
	if (($woc eq H1A)&&($wod eq MAG)){
	$mag++;
	}
	if (($woc eq SAO)&&($wod eq SDS)){
	$sds++;
	}
	if (($ff =~ m/^AMBER/)&&($woc eq O1)&&($wod eq DPP)){
	$dppg++;
	}
}
close (INP);
unless ($popc == 0){ 
print TOP "\nPOPC $popc";
}
unless ($pope == 0){
print TOP "\nPOPE $pope";
}
unless ($popg == 0){ 
print TOP "\nPOPG $popg";
}
unless ($pops == 0){
print TOP "\nPOPS $pops";
}
unless ($popa == 0){
print TOP "\nPOPA $popa";
}
unless ($dopc == 0){
print TOP "\nDOPC $dopc";
}
unless ($dppc == 0){ 
print TOP "\nDPPC $dppc";
}
unless ($dppe == 0){ 
print TOP "\nDPPE $dppe";
}
unless ($dppg == 0){ 
print TOP "\nDPPG $dppg";
}
unless ($ppcs == 0){ 
print TOP "\nPPCS $ppcs";
}
unless ($dspc == 0){
print TOP "\nDSPC $dspc";
}
unless ($dspe == 0){
print TOP "\nDSPE $dspe";
}
unless ($dspg == 0){
print TOP "\nDSPG $dspg";
}
unless ($dmpc == 0){
print TOP "\nDMPC $dmpc";
}
unless ($dmpe == 0){
print TOP "\nDMPE $dmpe";
}
unless ($dmpg == 0){
print TOP "\nDMPG $dmpg";
}
unless ($pip2 == 0){
print TOP "\nPIP2 $pip2";
}
unless ($pip3 == 0){
print TOP "\nPIP3 $pip3";
}
unless ($dhpc == 0){
print TOP "\nDHPC $dhpc";
}
unless ($bog == 0){ 
print TOP "\nBOG $bog\n";
}
unless ($ddm == 0){ 
print TOP "\nDDM $ddm\n";
}
unless ($lmpg == 0){ 
print TOP "\nLMPG $lmpg\n";
}
unless ($dpc == 0){ 
print TOP "\nDPC $dpc\n";
}
unless ($chol == 0){ 
print TOP "\nCHOL $chol\n";
}
unless ($card == 0){ 
print TOP "\nCARD $card\n";
}
unless ($pvpe == 0){
print TOP "\nPVPE $pvpe";
}
unless ($pvpg == 0){ 
print TOP "\nPVPG $pvpg";
}
unless ($lmg == 0){ 
print TOP "\nLMG $lmg\n";
}
unless ($sqd == 0){ 
print TOP "\nSQD $sqd";
}
unless ($dgd == 0){
print TOP "\nDGD $dgd";
}
unless ($dag == 0){
print TOP "\nDAG $dag";
}
unless ($mag == 0){
print TOP "\nMAG $mag";
}
unless ($sds == 0){
print TOP "\nSDS $sds";
}
close (TOP);
close (LIP);

if ($deletion eq "none"){
system("$gromacs/$gmx grompp -maxwarn 10 -f em-cg.mdp -o prot_lipid_em -c prot+lipid+box.pdb -p topol.top");
}
else{
system("$gromacs/$gmx grompp -maxwarn 10 -f em.mdp -o prot_lipid_em -c prot+lipid+box.pdb -p topol.top");
}
system("$gromacs/$gmx mdrun -deffnm prot_lipid_em");
system("$gromacs/$gmx trjconv -f prot_lipid_em.gro -o atomistic-protein-lipid.pdb -s prot_lipid_em.tpr -pbc res <<EOD
0
EOD");

open (ATM, "atomistic-protein-lipid.pdb") || die "Energy minimisation of AT protein-lipid complex failed\n";
close(ATM);

`cp atomistic-protein-lipid.pdb temp.pdb`;
	}
else{
`cp atomistic-protein.pdb temp.pdb`;
	}
}
else{
`cp atomistic-lipid.pdb temp.pdb`;
open (LIP, "lipid.top")  || die "Something went wrong\n";
open(TOP, "> topol.top") || die "No Output Possible\n";
my @lip = <LIP>;
if ($ff =~ m/^GROMOS/){
s/\[ system \]/\n\#include "gromos53a6.ff\/spc.itp\"\n\#include "gromos53a6.ff\/ions.itp\"\n\[ system \]\n/g for @lip;
}
if ($ff =~ m/^OPLS/){
s/\[ system \]/\n\#include "oplsaa.ff\/tip4p.itp\"\n\#include "oplsaa.ff\/ions.itp\"\n\[ system \]\n/g for @lip;
}
if ($ff =~ m/^CHARMM27/){
s/\[ system \]/\n\#include "charmm27.ff\/tip3p.itp\"\n\#include "charmm27.ff\/ions.itp\"\n\[ system \]\n/g for @lip;
}
if ($ff =~ m/^CHARMM36/){
s/\[ system \]/\n\#include "charmm36.ff\/tip3p.itp\"\n\#include "charmm36.ff\/ions.itp\"\n\[ system \]\n/g for @lip;
}
if ($ff =~ m/^AMBER/){
s/\[ system \]/\n\#include "amber99.ff\/spce.itp\"\n\#include "amber99.ff\/ions.itp\"\n\[ system \]\n/g for @lip;
}
print TOP @lip;
close (TOP);
close (LIP);
}

###################################################
#                                                 #
#                  ADD WATER                      #
#                                                 #
###################################################

unless ($solvate eq "no"){

### GENBOX ADDITION ###

if ($delete_waters eq "gmxsolvate"){
open VDW, ">", 'vdwradii.dat' or die;
print VDW "\nAAA  C     0.15";
print VDW "\n???  F     0.12";
print VDW "\n???  H     0.04";
print VDW "\n???  N     0.110";
print VDW "\n???  O     0.105";
print VDW "\n???  S     0.16";
print VDW "\nCHO  C2    0.15";
print VDW "\nPOP  CA1   0.5";
print VDW "\nPOP  CA2   0.5";
print VDW "\nPVP  CA1   0.5";
print VDW "\nPOP  C51   0.5";
print VDW "\nPOP  C52   0.5";
print VDW "\nPVP  CA2   0.5";
print VDW "\nDSP  CA1   0.5";
print VDW "\nDSP  CA2   0.5";
print VDW "\nDSP  CA3   0.5";
print VDW "\nDSP  CA4   0.5";
print VDW "\nPOP  C54   0.5";
print VDW "\nPOP  C55   0.5";
print VDW "\nPOP  H18R  0.7";	#### charmm POPC+POPS atoms
print VDW "\nPOP  H18S  0.7";	#### charmm POPC+POPS atoms
print VDW "\nPOP  H18T  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H17R  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H17S  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H16R  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H16S  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H15R  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H15S  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H16X  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H16Y  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H16Z  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H15X  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H15Y  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H14X  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H14Y  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H13X  0.7";   #### charmm POPC+POPS atoms
print VDW "\nPOP  H13Y  0.7";   #### charmm POPC+POPS atoms
print VDW "\nCHO  C27   0.5";
print VDW "\nCHO  C28   0.5";
print VDW "\nCHO  C29   0.5";
print VDW "\nDOP  CA1   0.5";
print VDW "\nDOP  CA2   0.5";
print VDW "\nDOP  CA3   0.5";
print VDW "\nDOP  CA4   0.5";
print VDW "\nDPP  C30   0.5";
print VDW "\nDPP  C31   0.5";
print VDW "\nDPP  C49   0.5";
print VDW "\nDPP  C50   0.5";
print VDW "\nDAG  C30   0.5";
print VDW "\nDAG  C31   0.5";
print VDW "\nDAG  C49   0.5";
print VDW "\nDAG  C50   0.5";
print VDW "\nMAG  C30   0.5";
print VDW "\nMAG  C31   0.5";
print VDW "\nDMP  C28   0.5";
print VDW "\nDMP  C29   0.5";
print VDW "\nDMP  C45   0.5";
print VDW "\nDMP  C46   0.5";
print VDW "\nPIP  CAA   0.5";
print VDW "\nPIP  CAB   0.5";
print VDW "\nPIP  CCO   0.5";
print VDW "\nPIP  CCN   0.5";
print VDW "\nDDM  C24   0.5";
print VDW "\nDDM  C23   0.5";
print VDW "\nBOG  C8B   0.5";
print VDW "\nBOG  C14   0.5";
print VDW "\nCAR  CA1   0.5";
print VDW "\nPPC  C32   0.5";
print VDW "\nPPC  C50   0.5";
print VDW "\nCAR  CA2   0.5\n";
close(VDW);
}

if ($ff =~ m/^GROMOS/){
system("$gromacs/$gmx solvate -cs spc216.gro -cp temp.pdb -o prot+lipid+sol.pdb -scale 1");
}
if ($ff =~ m/^OPLS/){
system("$gromacs/$gmx solvate -cs tip4p.gro -cp temp.pdb -o prot+lipid+sol.pdb -scale 1");
}
if ($ff =~ m/^CHARMM/){
system("$gromacs/$gmx solvate -cs spc216.gro -cp temp.pdb -o prot+lipid+sol.pdb -scale 1");
}
if ($ff =~ m/^AMBER/){
system("$gromacs/$gmx solvate -cs spc216.gro -cp temp.pdb -o prot+lipid+sol.pdb -scale 1");
}

system("rm vdwradii.dat");

### IF GENBOX ADDITION FAILS DUE TO MEMORY ISSUES ###

unless(-e "prot+lipid+sol.pdb"){
$delete_waters = "pymol";
if ($ff =~ m/^GROMOS/){
system("$gromacs/$gmx solvate -cs spc216.gro -cp temp.pdb -o solvated.pdb");
}
if ($ff =~ m/^OPLS/){
system("$gromacs/$gmx solvate -cs tip4p.gro -cp temp.pdb -o solvated.pdb");
}
if ($ff =~ m/^CHARMM/){
system("$gromacs/$gmx solvate -cs spc216.gro -cp temp.pdb -o solvated.pdb");
}
if ($ff =~ m/^AMBER/){
system("$gromacs/$gmx solvate -cs spc216.gro -cp temp.pdb -o solvated.pdb");
}
}

if ($delete_waters eq "pymol") {

open (PYMOL, "> delete_waters.pml") || die "Cannot open delete_waters.pml\n";
print PYMOL "set retain_order,1\n";
print PYMOL "set pdb_retain_ids,1\n";
print PYMOL "load solvated.pdb\n";
print PYMOL "select headgroup, bm. (name p*,o*,n*) and (resn pop*,pvp*,dop*,dpp*,ppc*,dsp*,dmp*,dhp*,pip*,bog*,cho*,ddm*,lmp*,dpc*,car*,lmg*,sqd,dgd,dag,sds) around 6 and resname sol
\n";
print PYMOL "remove bm. (name C* and (not name C1,C2,C3,C4,C5,C6)) and (resn pop*,pvp*,dop*,dpp*,ppc*,dsp*,dmp*,dhp*,pip*,bog*,cho*,ddm*,lmp*,dpc*,car*,lmg*,sqd,dgd,dag,sds) around 6 and resn sol and not headgroup\n";
print PYMOL "save sol_kill.pdb\n";
close (PYMOL);

system("$pymol -c delete_waters.pml");

system("cat box.pdb sol_kill.pdb > sol_kill_box.pdb");

system("$gromacs/$gmx editconf -resnr 1 -f sol_kill_box.pdb -o prot+lipid+sol.pdb");
}

elsif ($delete_waters eq "killwater"){
$solvated = "solvated.pdb";
$outfile = "sol_kill.pdb";
$outfile2 = "prot+lipid+sol.pdb";

$numWaterRemoved = 0;			

$protMinX = 1000;
$protMinY = 1000;
$protMinZ = 1000;
$protMaxX = 0;
$protMaxY = 0;
$protMaxZ = 0;

open(IN, $solvated);
	until(eof IN){
		$linebuffer = <IN>;
		@tokens = split ' ', $linebuffer;
		if((($tokens[0] eq "ATOM")&&($ff ne "AMBER")&&(($tokens[2] eq "C20")||($tokens[2] eq "C25")||($tokens[2] eq "C10")))||($tokens[0] eq "ATOM")&&($ff eq "AMBER")&&($tokens[2] eq "C30")){
			$atomZ = $tokens[7];
				if($atomZ < $protMinZ){
					$protMinZ = $atomZ;
				}
				if($atomZ > $protMaxZ){
					$protMaxZ = $atomZ;
				}
			}
		}


open(OUT, ">./$outfile");

$protMinX += 10;
$protMinY += 10;
$protMaxX -= 10;
$protMaxY -= 10;
			
close(IN);
open(IN, $solvated);

#now kill the waters
until(eof IN){
	$linebuffer = <IN>;
	@tokens = split ' ', $linebuffer;
	
	$atomX = $tokens[5];
	$atomY = $tokens[6];
	$atomZ = $tokens[7];

	if(($tokens[0] eq "ATOM") && ($tokens[3] eq "SOL") && 
	   ($tokens[2] eq "OW")  && ($atomZ > $protMinZ) && ($atomZ < $protMaxZ)){
		##ITS a water molecule in the kill z range, but is it in the protein box?
		if(($atomX > $protMinX) && ($atomX < $protMaxX) && ($atomY > $protMinY) && ($atomY < $protMaxY)){
			print "Water #$tokens[4] in protein and saved\n";
			print OUT $linebuffer;
		}
		else{
			$linebuffer = <IN>;
			$linebuffer = <IN>;
			if ($ff =~ m/^OPLS/){
			$linebuffer = <IN>;
			}			
			$numWaterRemoved++;
		}
	}
	else{
		print OUT $linebuffer;
	}
}

print "$numWaterRemoved Waters were removed\n";

close(IN);
close(OUT);

system("$gromacs/$gmx editconf -resnr 1 -f $outfile -o $outfile2");
	}

open (SOL, "prot+lipid+sol.pdb") || die "Cannot open prot+lipid+sol.pdb\n";

$sol = 0;

while (<SOL>){
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
	if ($woc eq OW){
	$sol++;
	}
}
close (SOL);

open(TOP, "topol.top") || die "Something went wrong\n";
open(SOLTOP, "> sol.top") || die "No Output Possible\n";
	my @top = <TOP>;
	if ($ff =~ m/^GROMOS/){
	s/\#include \"spc.itp\"/\#ifdef FLEX_SPC\n\#include \"gromos53a6.ff\/flexspc.itp\"\n\#else\n#include \"gromos53a6.ff\/spc.itp"\n\#endif/g for @top;
	}
	if ($ff =~ m/^OPLS/){
	s/\#include \"spc.itp\"/\#include \"oplsaa.ff\/tip4p.itp\"/g for @top;
	}
	if ($ff =~ m/^CHARMM/){
	s/\#include \"spc.itp\"/\#include \"charmm27.ff\/tip3p.itp\"/g for @top;
	}
	print SOLTOP @top;
        printf (SOLTOP "\nSOL %d\n", $sol );
close (TOP);
close (SOLTOP);

open(TOP, "topol.top") || die "Something went wrong\n";
open(OUT, "> topol1.top") || die "No Output Possible\n";
	my @top = <TOP>;
	if ($ff =~ m/^GROMOS/){
	s/\#include \"spc.itp\"/\#ifdef FLEX_SPC\n\#include \"gromos53a6.ff\/flexspc.itp\"\n\#else\n#include \"gromos53a6.ff\/spc.itp"\n\#endif/g for @top;
	}
	if ($ff =~ m/^OPLS/){
	s/\#include \"spc.itp\"/\#include \"oplsaa.ff\/tip4p.itp\"/g for @top;
	}
	if ($ff =~ m/^CHARMM27/){
	s/\#include \"spc.itp\"/\#include \"charmm27.ff\/tip3p.itp\"/g for @top;
	}
	if ($ff =~ m/^CHARMM36/){
	s/\#include \"spc.itp\"/\#include \"charmm36.ff\/tip3p.itp\"/g for @top;
	}
	print OUT @top;
close (TOP);

###################################################
#                                                 #
#                   ADD IONS                      #
#                                                 #
####################################################

`$gromacs/$gmx grompp -maxwarn 10 -f em.mdp -o em -c prot+lipid+sol.pdb -p sol.top`;

system("$gromacs/$gmx make_ndx -f prot+lipid+sol.pdb -o solvent.ndx << EOD
del 0-25
r SOL
q
EOD");

system("$gromacs/$gmx genion -s em.tpr -o prot+lipid+sol+ion.pdb -neutral -conc 0.15 -n solvent.ndx << EOD
0
EOD");

my $na = 0;
my $sol = 0;
my $cl = 0;

open (INP, "prot+lipid+sol+ion.pdb") || die "Cannot open prot+lipid+box.pdb\n";
$atom = 0;
while (<INP>){
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
	if (($woc eq NA)){
	$na++;
	}
	if (($woc eq CL)){
	$cl++;
	}
	if (($woc eq OW)){
	$sol++;
	}
}
close (INP);

print OUT "\nSOL $sol";
print OUT "\nNA $na";
print OUT "\nCL $cl";

close (OUT);

open(EMV, "em.mdp") || die "Something went wrong\n";
open(EMS, "> em-system.mdp") || die "No Output Possible\n";

while (<EMV>) {
if ($_ !~ /epsilon_r/) {
	print EMS $_;
	}
}
print EMS "\ndefine                   = -DFLEX_SPC\n";
close (EMS);
close (EMV);

system("$gromacs/$gmx grompp -maxwarn 10 -f em-system.mdp -o system_em -c prot+lipid+sol+ion.pdb -p topol1.top");
system("$gromacs/$gmx mdrun -deffnm system_em");
system("$gromacs/$gmx trjconv -f system_em.gro -o atomistic-system.pdb -s system_em.tpr -pbc res <<EOD
0
EOD");

if (($lipidnum > 0)&&($protnum > 0)) {
system("$gromacs/$gmx make_ndx -f atomistic-system.pdb -o atomistic-system.ndx <<EOD
del 0
del 1-30
rSOL|aNA|aCL
name 1 SOL_ION
0|1
2|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
!2|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
name 4 Lipid
del 2
del 2
q
EOD");
`cp $cg2atdir/mdp/pr-gmx5.mdp pr.mdp`;
`cp $cg2atdir/mdp/100ns-gmx5.mdp 100ns.mdp`;
} elsif (($lipidnum > 0)&&($protnum == 0)){
system("$gromacs/$gmx make_ndx -f atomistic-system.pdb -o atomistic-system.ndx <<EOD
del 0-30
rSOL|aNA|aCL
name 0 SOL_ION
0|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
!0|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDPC*|rDSP*|rCAR*|rPVP*|rLMG*|rLMP*|rSQD*|rDGD*|rDAG*|rSDS*|rPPC*|rMAG*
name 2 Lipid
del 1
q
EOD");
`cp $cg2atdir/mdp/pr-lipid-gmx5.mdp pr.mdp`;
`cp $cg2atdir/mdp/100ns-lipid-gmx5.mdp 100ns.mdp`;
} elsif (($lipidnum == 0)&&($protnum > 0)){
system("$gromacs/$gmx make_ndx -f atomistic-system.pdb -o atomistic-system.ndx <<EOD
del 0
del 1-30
rSOL|aNA|aCL
name 1 SOL_ION
q
EOD");
`cp $cg2atdir/mdp/pr-protein-gmx5.mdp pr.mdp`;
`cp $cg2atdir/mdp/100ns-protein-gmx5.mdp 100ns.mdp`;
}
`mkdir PR/`;

mkdir('itp_files');
`cp *.itp atomistic-system.ndx itp_files`;
`mv topol1.top topol`;

}

######### Cleaning up ##########
`rm -rf ddm.pdb dspe.pdb order-* cg.top em.tpr martini* ff_v1.4_x.itp chol.pdb Res_CG/ fit_*/ Sidechains/ *.ndx *v*.itp z *.log wat.pdb order.pdb box.pdb ion_atoms.pdb lipid.pdb p2g.pdb em_inbox* *_em* *bilayer* '#'* *_chain* NA* CL* pip* pop* dpp* bog* ppc* temp* sol* *.py  dhp* dmp* *label* mdout.mdp model.pdb chain* dop* pro*.pdb inbox.pdb coord.xvg charge.txt sqd* dgd* lmg* pvp* *old* ener.edr dag* sds*`;
`mv topol topol.top`;
`mv itp_files/* .`;
`rm -rf itp_files #* *membed* *.cpt nohydrogens.pdb order-center* em.tpr ddm* ener* chol*`;
`rm -rf cg.top traj.trr nohydrogens.pdb lipid.top membed* ddm.top dpc.pdb card.pdb protein.top protein-cg.itp ff_v1.4_x.itp martini_v2.2.itp alc*`;
mkdir('mdp_files');
`mv *.mdp mdp_files/`;
#################################


