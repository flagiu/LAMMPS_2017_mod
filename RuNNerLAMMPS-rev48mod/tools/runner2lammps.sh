#!/bin/bash

####################################################################################################
#
# File converter RuNNer DATA -> LAMMPS DATA
#
# author:   Andreas Singraber
# date:     2013-02-28
# email:    andreas.singraber@univie.ac.at
#
# Copyright (c) 2013 Andreas Singraber
#
####################################################################################################

awk 'BEGIN {
	UC_LENGTH = 0.529177208613983;
    #UC_LENGTH = 1.0;
	ilat = 0;
	i = 0;
	j = 0;
    na = 0;
	ne = 0;
    iatom = 0;
    lelem[  1] = "H" ;
    lelem[  2] = "He";
    lelem[  3] = "Li";
    lelem[  4] = "Be";
    lelem[  5] = "B" ;
    lelem[  6] = "C" ;
    lelem[  7] = "N" ;
    lelem[  8] = "O" ;
    lelem[  9] = "F" ;
    lelem[ 10] = "Ne";
    lelem[ 11] = "Na";
    lelem[ 12] = "Mg";
    lelem[ 13] = "Al";
    lelem[ 14] = "Si";
    lelem[ 15] = "P" ;
    lelem[ 16] = "S" ;
    lelem[ 17] = "Cl";
    lelem[ 18] = "Ar";
    lelem[ 19] = "K" ;
    lelem[ 20] = "Ca";
    lelem[ 21] = "Sc";
    lelem[ 22] = "Ti";
    lelem[ 23] = "V" ;
    lelem[ 24] = "Cr";
    lelem[ 25] = "Mn";
    lelem[ 26] = "Fe";
    lelem[ 27] = "Co";
    lelem[ 28] = "Ni";
    lelem[ 29] = "Cu";
    lelem[ 30] = "Zn";
    lelem[ 31] = "Ga";
    lelem[ 32] = "Ge";
    lelem[ 33] = "As";
    lelem[ 34] = "Se";
    lelem[ 35] = "Br";
    lelem[ 36] = "Kr";
    lelem[ 37] = "Rb";
    lelem[ 38] = "Sr";
    lelem[ 39] = "Y" ;
    lelem[ 40] = "Zr";
    lelem[ 41] = "Nb";
    lelem[ 42] = "Mo";
    lelem[ 43] = "Tc";
    lelem[ 44] = "Ru";
    lelem[ 45] = "Rh";
    lelem[ 46] = "Pd";
    lelem[ 47] = "Ag";
    lelem[ 48] = "Cd";
    lelem[ 49] = "In";
    lelem[ 50] = "Sn";
    lelem[ 51] = "Sb";
    lelem[ 52] = "Te";
    lelem[ 53] = "I" ;
    lelem[ 54] = "Xe";
    lelem[ 55] = "Cs";
    lelem[ 56] = "Ba";
    lelem[ 57] = "La";
    lelem[ 58] = "Ce";
    lelem[ 59] = "Pr";
    lelem[ 60] = "Nd";
    lelem[ 61] = "Pm";
    lelem[ 62] = "Sm";
    lelem[ 63] = "Eu";
    lelem[ 64] = "Gd";
    lelem[ 65] = "Tb";
    lelem[ 66] = "Dy";
    lelem[ 67] = "Ho";
    lelem[ 68] = "Er";
    lelem[ 69] = "Tm";
    lelem[ 70] = "Yb";
    lelem[ 71] = "Lu";
    lelem[ 72] = "Hf";
    lelem[ 73] = "Ta";
    lelem[ 74] = "W" ;
    lelem[ 75] = "Re";
    lelem[ 76] = "Os";
    lelem[ 77] = "Ir";
    lelem[ 78] = "Pt";
    lelem[ 79] = "Au";
    lelem[ 80] = "Hg";
    lelem[ 81] = "Tl";
    lelem[ 82] = "Pb";
    lelem[ 83] = "Bi";
    lelem[ 84] = "Po";
    lelem[ 85] = "At";
    lelem[ 86] = "Rn";
    lelem[ 87] = "Fr";
    lelem[ 88] = "Ra";
    lelem[ 89] = "Ac";
    lelem[ 90] = "Th";
    lelem[ 91] = "Pa";
    lelem[ 92] = "U" ;
    lelem[ 93] = "Np";
    lelem[ 94] = "Pu";
    lelem[ 95] = "Am";
    lelem[ 96] = "Cm";
    lelem[ 97] = "Bk";
    lelem[ 98] = "Cf";
    lelem[ 99] = "Es";
    lelem[100] = "Fm";
    lelem[101] = "Md";
    lelem[102] = "No";
}
{
	if($1 == "lattice") {
		ilat++;
		olat[ilat,1] = $2 * UC_LENGTH;
        olat[ilat,2] = $3 * UC_LENGTH;
        olat[ilat,3] = $4 * UC_LENGTH;
	}
	if($1 == "atom") {
		iatom++;
		x[iatom] = $2 * UC_LENGTH;
        y[iatom] = $3 * UC_LENGTH;
        z[iatom] = $4 * UC_LENGTH;
        es[iatom] = $5;
	}
	if($1 == "comment") {
		$1 = "";
		comment = $0
	}
}
END {
	na = iatom;

    # DETERMINE NUMBER OF ELEMENTS
	ne = 1;
	ename[1] = es[1];
	for(i=2; i<=na; i++) {
		j=1;
		while(j <= ne) {
			if(es[i] == ename[j]) {
				break;
			}
			j++;
		}
		if(j == ne+1) {
			ename[j] = es[i];
			ne++;
		}
	}

    # GET NUCLEAR CHARGE OF ALL ELEMENTS
    for(i=1; i<=ne; i++) {
        for(j=1; j<=102; j++) {
            if(ename[i] == lelem[j]) {
                nucch[i] = j;
                break;
            }
        }
        et[i] = i;
    }

    # SORT ELEMENTS IN ORDER OF INCREASING NUCLEAR CHARGE
    i = 1;
    if(ne > 1) {
        do {
            if(nucch[i] > nucch[i+1]) {
                tmp = et[i];
                et[i] = et[i+1];
                et[i+1] = tmp;
                tmp = nucch[i];
                nucch[i] = nucch[i+1];
                nucch[i+1] = tmp;
                tmp = ename[i];
                ename[i] = ename[i+1];
                ename[i+1] = tmp;
                i = 0;
            }
            i++;
        } while(i<ne);
    }

    # ASSIGN LAMMPS ELEMENT INDEX TO ATOMS
	for(i=1; i<=na; i++) {
		for (j=1; j<=ne; j++) {
			if(es[i] == ename[j]) {
                ei[i] = j;
                nape[j]++;
            }
		}
	}

    # APPLY LAMMPS LATTICE VECTOR CONVENTIONS
    A_x = olat[1,1]; A_y = olat[1,2]; A_z = olat[1,3];
    B_x = olat[2,1]; B_y = olat[2,2]; B_z = olat[2,3];
    C_x = olat[3,1]; C_y = olat[3,2]; C_z = olat[3,3];
    norm_A = sqrt(A_x * A_x + A_y * A_y + A_z * A_z);
    norm_B = sqrt(B_x * B_x + B_y * B_y + B_z * B_z);
    norm_C = sqrt(C_x * C_x + C_y * C_y + C_z * C_z);
    a_x = norm_A;
    a_y = 0.0;
    a_z = 0.0;
    b_x = 1.0 / norm_A * (B_x * A_x + B_y * A_y + B_z * A_z);
    b_y = sqrt(norm_B * norm_B - b_x * b_x);
    b_z = 0.0;
    c_x = 1.0 / norm_A * (C_x * A_x + C_y * A_y + C_z * A_z);
    c_y = 1.0 / b_y * (B_x * C_x + B_y * C_y + B_z * C_z - b_x * c_x);
    c_z = sqrt(norm_C * norm_C - c_x * c_x - c_y * c_y);
    nlat[1,1] = a_x; nlat[1,2] = a_y; nlat[1,3] = a_z; 
    nlat[2,1] = b_x; nlat[2,2] = b_y; nlat[2,3] = b_z; 
    nlat[3,1] = c_x; nlat[3,2] = c_y; nlat[3,3] = c_z; 

    # WRITE OUTPUT
	printf "Generated by runner2lammps.sh, original comment line: %s\n", comment;
	printf "\n";
	printf "%5d atoms\n", na;
	printf "%3d atom types\n", ne;
	printf "\n";
	printf "%24.16E %24.16E xlo xhi\n", 0.0, nlat[1,1];
	printf "%24.16E %24.16E ylo yhi\n", 0.0, nlat[2,2];
	printf "%24.16E %24.16E zlo zhi\n", 0.0, nlat[3,3];
	if(nlat[2,1] != 0 || nlat[3,1] != 0 || nlat[3,2] !=0) {
        printf "%24.16E %24.16E %24.16E xy xz yz\n", nlat[2,1], nlat[3,1], nlat[3,2];
    }
	printf "\n";
	printf "Atoms\n";
	printf "\n";
	for(i=1; i<=na; i++) {
	    printf "%5d %3d %24.16E %24.16E %24.16E\n", i, ei[i], x[i], y[i], z[i];
	}
}' 
