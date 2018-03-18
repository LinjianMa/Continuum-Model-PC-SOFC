//*****************************************************************************************************************//% 
	//------------------------------Transport of electron in anode---------------------------
	forAll(anodemesh.C(),counter)
	{
		grad_phi_A.internalField()[counter].y() = D_e_A.value()*96485/8.31/873*Efield_A.internalField()[counter];           
	}
	dphi_e_A = linearInterpolate(grad_phi_A) & anodemesh.Sf();
	C_e_A.storePrevIter();
	solve
		(	ratio_e_A*fvm::ddt(C_e_A) - fvm::laplacian(D_e_A,C_e_A) - fvm::div(dphi_e_A, C_e_A)	);
	C_e_A.relax();
	Current_e_A = 96485*(D_e_A*fvc::grad(C_e_A)+grad_phi_A*C_e_A);
	//------------------------------Transport of electron in cathode---------------------------
	forAll(cathodemesh.C(),counter)
	{
		grad_phi_C.internalField()[counter].y() = D_e_C.value()*96485/8.31/873*Efield_C.internalField()[counter];           
	}
	dphi_e_C = linearInterpolate(grad_phi_C) & cathodemesh.Sf();
	C_e_C.storePrevIter();
	solve
		(	ratio_e_C*fvm::ddt(C_e_C) - fvm::laplacian(D_e_C,C_e_C) - fvm::div(dphi_e_C, C_e_C)	);
	C_e_C.relax();
	Current_e_C = 96485*(D_e_C*fvc::grad(C_e_C)+grad_phi_C*C_e_C);
//*****************************************************************************************************************//