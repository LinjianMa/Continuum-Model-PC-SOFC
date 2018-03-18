//*****************************************************************************************************************//
	//-------------------------Solving the Poisson Equations-----------------------------
	Phi.storePrevIter();
	solve
		(	fvm::laplacian(Phi) == -charge/perm_0	);
	Phi.relax();
	Efield = -fvc::grad(Phi);                   //Electric field output
	//---------------------------Map Phi to Phi_ele---------------------------
	solve
	(	fvm::ddt(Phi_ele)	);
	solve
	(	fvm::ddt(Efield_ele)	);
	solve
	(	fvm::ddt(Efield_bulk)	);             //Electric field in bulk electrolyte
	forAll(electrolytemesh.C(),counter)
	{
		Phi_ele.internalField()[counter] = Phi.internalField()[counter+grid_ele_L];
		Efield_ele.internalField()[counter] = Efield.internalField()[counter+grid_ele_L].y();            //Electric field output
	}
	Efield_bulk.internalField()[0] = Efield_ele.internalField()[gridnum_ele/2];
	Efield_bulk.internalField()[gridnum_ele-1] = Efield_ele.internalField()[gridnum_ele/2];
	//---------------------------Map Phi to Phi_A---------------------------
	solve
	(	fvm::ddt(Phi_A)	);
	solve
	(	fvm::ddt(Efield_A)	);
	solve
	(	fvm::ddt(Efield_bulk_A)	);             //Electric field in Anode
	forAll(anodemesh.C(),counter)
	{
		Phi_A.internalField()[counter] = Phi.internalField()[counter];
		Efield_A.internalField()[counter] = Efield.internalField()[counter].y();            //Electric field output
	}
	Efield_bulk_A.internalField()[0] = Efield_A.internalField()[gridnum_a/2];
	Efield_bulk_A.internalField()[gridnum_a-1] = Efield_A.internalField()[gridnum_a/2];
	//---------------------------Map Phi to Phi_C---------------------------
	solve
	(	fvm::ddt(Phi_C)	);
	solve
	(	fvm::ddt(Efield_C)	);
	solve
	(	fvm::ddt(Efield_bulk_C)	);             //Electric field in Cathode
	forAll(cathodemesh.C(),counter)
	{
		Phi_C.internalField()[counter] = Phi.internalField()[counter+grid_ele_R+gridnum_6];
		Efield_C.internalField()[counter] = Efield.internalField()[counter+grid_ele_R+gridnum_6].y();            //Electric field output
	}
	Efield_bulk_C.internalField()[0] = Efield_C.internalField()[gridnum_c/2];
	Efield_bulk_C.internalField()[gridnum_c-1] = Efield_C.internalField()[gridnum_c/2];
	//-------------------------Map Phi to overpotential_A and Overpotential_C---------------
	solve
	(	fvm::ddt(overpotential_A)	);
	solve
	(	fvm::ddt(overpotential_C)	);
	overpotential_A.internalField()[gridnum_a-1] = Phi.internalField()[grid_ele_L+gridnum_3-1] - Phi.internalField()[gridnum_0];
	//if (overpotential_A.internalField()[gridnum_0-1] > 0.49 )
	//{
	//	overpotential_A.internalField()[gridnum_0-1] = 0.49;
	//}
	overpotential_C.internalField()[0] = Phi.internalField()[grid_ele_R+gridnum_6+gridnum_7-1] - Phi.internalField()[grid_ele_R-gridnum_5];
//*****************************************************************************************************************//