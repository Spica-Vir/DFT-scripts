#!/bin/bash
###################################################
# Formate transformation: geometry optimization   #
# output file (.gui, or fort.34) of CRYSTAL-17 to #
# POSCAR file of VASP. Place this file in work    #
# dictionary.                                     #
# Note: the .xyz file is problematic in VESTA.    #
# Spica.H.Zhou,ICL,spica.h.zhou@gmail.com         #
# 20:09:57,04/05/21                               #
###################################################
# Note: Use atomic positions reported in .gui or fort.34, or .f34
FILENAME=$1 # output files of geometry optimization
OUTPUTNAME=$2
WORKPATH=`pwd`
# rm $OUTPUTNAME
# 
# label=`grep 'OPT END - CONVERGED' $WORKPATH/$FILENAME.out`
# if [[ ! -n "$label"  ]]; then
# 	echo "Error! Geometry optimization not converged!"
# 	exit
# fi
# lattice data
avect=(`awk 'NR==2 {print $1,$2,$3}' $WORKPATH/$FILENAME`)
bvect=(`awk 'NR==3 {print $1,$2,$3}' $WORKPATH/$FILENAME`)
cvect=(`awk 'NR==4 {print $1,$2,$3}' $WORKPATH/$FILENAME`)
# Generate POSCAR
echo "Transformated from CRYSTAL-17" > $WORKPATH/$OUTPUTNAME
echo "1.00" >> $WORKPATH/$OUTPUTNAME
echo "    ${avect[@]}" >> $WORKPATH/$OUTPUTNAME
echo "    ${bvect[@]}" >> $WORKPATH/$OUTPUTNAME
echo "    ${cvect[@]}" >> $WORKPATH/$OUTPUTNAME
echo "Cartesian" >> $WORKPATH/$OUTPUTNAME
# skip symmetric operators
nop=`awk 'NR==5 {printf "%d",$1}' $WORKPATH/$FILENAME`
lineop=`echo "scale=0;$nop*4+5" | bc`
# atomic positions
	# generate periodic table (that's crazy) up to 83, radioelement after that
pridtb=('X' 'H' 'He' 'Li' 'Be' 'B' 'C' 'N' 'O' 'F' 'Ne' 'Na' 'Mg' 'Al' 'Si' 'P'\
	'S' 'Cl' 'Ar' 'K' 'Ca' 'Sc' 'Ti' 'V' 'Cr' 'Mn' 'Fe' 'Co' 'Ni' 'Cu' 'Zn' \
	'Ga' 'Ge' 'As' 'Se' 'Br' 'Kr' 'Rb' 'Sr' 'Y' 'Zr' 'Nb' 'Mo' 'Tc' 'Ru' 'Rh' \
	'Pd' 'Ag' 'Cd' 'In' 'Sn' 'Sb' 'Te' 'I' 'Xe' 'Cs' 'Ba' 'La' 'Ce' 'Pr' 'Nd' \
	'Pm' 'Sm' 'Eu' 'Gd' 'Tb' 'Dy' 'Ho' 'Er' 'Tm' 'Yb' 'Lu' 'Hf' 'Ta' 'W' 'Re' \
	'Os' 'Ir' 'Pt' 'Au' 'Hg' 'Tl' 'Pb' 'Bi')
#
natom=`awk -v a="$lineop" 'NR==a+1 {printf "%d",$1}' $WORKPATH/$FILENAME`
bglineatom=`echo "scale=0;$lineop+2" | bc`
edlineatom=`echo "scale=0;$natom+$lineop+1" | bc`
#
label=()
num=()
nkind=0
for (( i = $bglineatom; i <= $edlineatom; i++ )); do
	atompos=(`awk -v a="$i" 'NR==a {printf "%10.6f%10.6f%10.6f",$2,$3,$4}' \
		$WORKPATH/$FILENAME`)
	echo "    ${atompos[@]}" >> $WORKPATH/$OUTPUTNAME
	atomlabel=`awk -v a="$i" 'NR==a {printf "%d",$1}' $WORKPATH/$FILENAME`
	# too large or conventional atomic number
	if [[ $atomlabel > 83 ]] && [[ $atomlabel < 101 ]]; then
		echo "Error: Atomic number too large at line $i , enlarge the periodic \
		table!"
		exit
	fi
	while (( $atomlabel > 100 )); do
		atomlabel=`echo "scale=0;$atomlabel-100" | bc`
	done
	# classify and count atoms
	if [[ $nkind == 0 ]]; then
		label[$nkind]=$atomlabel
		num[$nkind]=1
		nkind=`echo "scale=0;$nkind+1" | bc`
	else
		isnew=1
		for (( j = 0; j < ${#label[@]}; j++ )); do
			if [[ ${label[$j]} == $atomlabel ]]; then
				isnew=0
				num[$j]=`echo "scale=0;${num[$j]}+1" | bc`
			else
				continue
			fi
		done
		if [[ $isnew == 1 ]]; then
				label[$nkind]=$atomlabel
				num[$nkind]=1
				nkind=`echo "scale=0;$nkind+1" | bc`
		fi
	fi
done
# enter atom names and numbers
for i in ${label[@]}; do
	name=`echo "$name    ${pridtb[$i]}"`
done
for i in ${num[@]}; do
	num_str=`echo "$num_str    $i"`
done
sed -i "/Cartesian/i\\$name" $WORKPATH/$OUTPUTNAME
sed -i "/Cartesian/i\\$num_str" $WORKPATH/$OUTPUTNAME
