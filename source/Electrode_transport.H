//*****************************************************************************************************************//
	//----------------------Solve the anode and cathode transport equations-------------------
	forAll(anodemesh.C(),counter)
	{
		//Anode
		V1_A.internalField()[counter] = 1e-6*(k1A_plus*C_H2.internalField()[counter] - k1A_minus*C_H2ads.internalField()[counter]);
		V2_A.internalField()[counter] = 1e-6*(k2A_plus*C_H2ads.internalField()[counter] - k2A_minus*C_H_A.internalField()[counter]*C_H_A.internalField()[counter]/C_max);
	}
	forAll(cathodemesh.C(),counter)
	{
		//Cathode
		V2_C.internalField()[counter] = 1e-6*(k2C_plus*C_O2.internalField()[counter] - k2C_minus*C_O2ads.internalField()[counter]);
		V3_C.internalField()[counter] = 1e-6*(k3C_plus*C_O2ads.internalField()[counter] - k3C_minus*C_O.internalField()[counter]*C_O.internalField()[counter]/C_max);
		V4_C.internalField()[counter] = 1e-6*(k4C_plus*C_H_C.internalField()[counter]*C_O.internalField()[counter]/C_max - k4C_minus*C_OH.internalField()[counter]);
		V5_C.internalField()[counter] = 1e-6*(k5C_plus*C_H_C.internalField()[counter]*C_OH.internalField()[counter]/C_max - k5C_minus*C_H2Oad.internalField()[counter]);
		V6_C.internalField()[counter] = 1e-6*(k6C_plus*C_H2Oad.internalField()[counter] - k6C_minus*C_H2O.internalField()[counter]);	
	}
	//-----------------------------Anode----------------------------------
	C_H2.storePrevIter();
	solve
		(	ratio_H2*fvm::ddt(C_H2) - fvm::laplacian(D_H2,C_H2) == -V1_A	);
	C_H2.relax();
	C_H2ads.storePrevIter();
	solve
		(	ratio_H2ads*fvm::ddt(C_H2ads) - fvm::laplacian(D_H2ads,C_H2ads) == V1_A-V2_A	);
	C_H2ads.relax();
	C_H_A.storePrevIter();
	solve
		(	ratio_H_A*fvm::ddt(C_H_A) - fvm::laplacian(D_H_A,C_H_A) == 2*V2_A	);
	C_H_A.relax();
	Flux_H_A = -D_H_A*fvc::grad(C_H_A);
	solve
		(	fvm::laplacian(D_H_ele,current_boundary_A)	);  //update anode current density
	//-----------------------------Anode End------------------------------
	//-----------------------------Cathode--------------------------------
	C_H_C.storePrevIter();
	solve
		(	ratio_H_C*fvm::ddt(C_H_C) - fvm::laplacian(D_H_C,C_H_C) == -V4_C-V5_C	);
	C_H_C.relax();
	C_O2.storePrevIter();
	solve
		(	ratio_O2*fvm::ddt(C_O2) - fvm::laplacian(D_O2,C_O2) == -V2_C	);
	C_O2.relax();
	C_O2ads.storePrevIter();
	solve
		(	ratio_O2ads*fvm::ddt(C_O2ads) - fvm::laplacian(D_O2ads,C_O2ads) == V2_C-V3_C	);
	C_O2ads.relax();
	C_O.storePrevIter();
	solve
		(	ratio_O*fvm::ddt(C_O) - fvm::laplacian(D_O,C_O) == 2*V3_C-V4_C	);
	C_O.relax();
	C_OH.storePrevIter();
	solve
		(	ratio_OH*fvm::ddt(C_OH) - fvm::laplacian(D_OH,C_OH) == V4_C-V5_C	);
	C_OH.relax();
	C_H2Oad.storePrevIter();
	solve
		(	ratio_H2Oad*fvm::ddt(C_H2Oad) - fvm::laplacian(D_H2Oad,C_H2Oad) == V5_C-V6_C	);
	C_H2Oad.relax();
	C_H2O.storePrevIter();
	solve
		(	ratio_H2O*fvm::ddt(C_H2O) - fvm::laplacian(D_H2O,C_H2O) == V6_C	);
	C_H2O.relax();
	Flux_H_C = -D_H_C*fvc::grad(C_H_C);
	Flux_O = -D_O*fvc::grad(C_O);
	Flux_OH = -D_OH*fvc::grad(C_OH);
	solve
		(	fvm::laplacian(D_H_ele,current_boundary_C)	);  //update cathode current density
	/*solve
		(	fvm::ddt(current_boundary_C1)	);
	solve
		(	fvm::ddt(current_boundary_C2)	);
	solve
		(	fvm::ddt(current_boundary_C3)	);*/
	//-----------------------------Cathode End------------------------------
	//---------Link anode and cathode current density with electrolyte I---------
	solve
		(	fvm::ddt(I_A)	);
	solve
		(	fvm::ddt(I_C)	);
	//solve
	//	(	fvm::ddt(ratio_H_anode)	);
	//solve
	//	(	fvm::ddt(ratio_H_cathode)	);
	I_A.internalField()[0] = current_boundary_A.internalField()[gridnum_0-1];
	I_C.internalField()[gridnum_ele-1] = current_boundary_C.internalField()[0];
	/*
	if (I_A.internalField()[0] < 5000 && I_A.internalField()[0] > -5000)
	{
		ratio_H_anode.internalField()[0] = 1;
		ratio_H_cathode.internalField()[gridnum_ele-1] = 1;
	}
	else 
	{	
		ratio_H_anode.internalField()[0] = 1e-1;
		ratio_H_cathode.internalField()[gridnum_ele-1] = 1e-1;
	}
	*/
//*****************************************************************************************************************//