/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    BZY_whole

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields_anode.H"
    #include "createFields_cathode.H"
    #include "createFields_electrolyte.H"
    #include "createFields_overall.H"
    #include "readSolverControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
	int	   step_count = 0;
	int    counter2 = 0;
	double C_H_equi = 1870.8; //0.0225*83147;
	double C_V_equi = 4677.02; //0.05625*83147;
	double C_H_trap = 1247.2; //0.015*83147;
	double C_eA_equi = 97342; //0.015*83147;
	double C_eC_equi = 97342; //0.015*83147;
	double grid_ele_L = gridnum_0+gridnum_1+gridnum_2;
	double grid_ele_R = gridnum_0+gridnum_1+gridnum_2+gridnum_3+gridnum_4+gridnum_5;
	double gridnum_ele = gridnum_3+gridnum_4+gridnum_5;
	double gridnum_a = gridnum_0+gridnum_1;
	double gridnum_c = gridnum_7+gridnum_8;
	double permittivity;
	double eqnResidualPhi;
	double eqnResidualH;
	double eqnResidualV;
	double ratio_surface;

//***************************************************************************************************************//	
if (runTime.timeName() == "0")
{
	//------------------------anode initialization-----------------------------
	forAll(anodemesh.C(),counter)
	{
		//concentration 
		C_H2ads.internalField()[counter] = k1A_plus/k1A_minus*C_H2.internalField()[counter];
		C_H_A.internalField()[counter] = Foam::sqrt ( k2A_plus/k2A_minus*C_H2ads.internalField()[counter]*C_max ) ;
		//initialize the reaction rates
		V1_A.internalField()[counter] = 0;
		V2_A.internalField()[counter] = 0;
	}
	//------------------------cathode initialization-------------------------------
	forAll(cathodemesh.C(),counter)
	{
		//concentration
		C_O2ads.internalField()[counter] = k2C_plus/k2C_minus*C_O2.internalField()[counter];
		C_O.internalField()[counter] = Foam::sqrt ( k3C_plus/k3C_minus*C_O2ads.internalField()[counter]*C_max ) ;
		C_H2Oad.internalField()[counter] = k6C_minus/k6C_plus*C_H2O.internalField()[counter];
		C_OH.internalField()[counter] = Foam::sqrt (k4C_plus*k5C_minus/k5C_plus/k4C_minus*C_O.internalField()[counter]*C_H2Oad.internalField()[counter]) ;
		C_H_C.internalField()[counter] = k4C_minus/k4C_plus*C_max*C_OH.internalField()[counter]/C_O.internalField()[counter];
		//initialize the reaction rates
		V2_C.internalField()[counter] = 0;
		V3_C.internalField()[counter] = 0;
		V4_C.internalField()[counter] = 0;
		V5_C.internalField()[counter] = 0;
		V6_C.internalField()[counter] = 0;
	}

	//------------------------charge initialization-------------------------------
    #include "charge_calculation_init.H"
}
//*****************************************************************************************************************//
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	//#include "Poisson.H"
    //#include "Electron_transport.H"
    //#include "charge_calculation.H"

	//do
	//{
	if (Efield.internalField()[grid_ele_R+gridnum_6+gridnum_7+100].y() <-0.02 || Efield.internalField()[grid_ele_R+gridnum_6+gridnum_7+100].y() >0.02 || Efield.internalField()[100].y() <-0.02 || Efield.internalField()[100].y() >0.02)
	{
		#include "Electron_transport.H"
		#include "charge_calculation.H"
		#include "Poisson.H"
	}
	else {
		#include "Electrode_transport.H"
		#include "Electrolyte_transport.H"
		#include "charge_calculation.H"
		#include "Poisson.H"
	}
	//}
	//while (Efield.internalField()[grid_ele_R+gridnum_5+100].y() <-1e-4 || Efield.internalField()[grid_ele_R+gridnum_5+100].y() >1e-4 || Efield.internalField()[100].y() <-1e-4 || Efield.internalField()[100].y() >1e-4);
	runTime.write();
	step_count++;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
