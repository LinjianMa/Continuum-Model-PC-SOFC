//*****************************************************************************************************************//
	//------------------------------Transport of proton---------------------------
	forAll(electrolytemesh.C(),counter)
	{
		grad_phi_ele.internalField()[counter].y() = -D_H_ele.value()*96485/8.31/873*Efield_ele.internalField()[counter];           
	}
	//grad_phi_ele = D_H_ele*96485/8.31/873*fvc::grad(Phi_ele);
	dphi_H = linearInterpolate(grad_phi_ele) & electrolytemesh.Sf();
	C_H_ele.storePrevIter();
    fvScalarMatrix HEqn
		(	ratio_H_ele*fvm::ddt(C_H_ele) - fvm::laplacian(D_H_ele,C_H_ele) - fvm::div(dphi_H, C_H_ele)	);
	C_H_ele.relax();
	Current_H_ele = 96485*(D_H_ele*fvc::grad(C_H_ele)+grad_phi_ele*C_H_ele);
    eqnResidualH = HEqn.solve().initialResidual(); 
	//------------------------------Transport of vacancy--------------------------
	/*
	forAll(electrolytemesh.C(),counter)
	{
		grad_phi_ele.internalField()[counter].y() = -2*D_V_ele.value()*96485/8.31/873*Efield_ele.internalField()[counter];         
	}
    //grad_phi_ele = 2*D_V_ele*96485/8.31/873*fvc::grad(Phi_ele);
	dphi_V = linearInterpolate(grad_phi_ele) & electrolytemesh.Sf();
	C_V_ele.storePrevIter();
    fvScalarMatrix VEqn
		(	ratio_V_ele*fvm::ddt(C_V_ele) - fvm::laplacian(D_V_ele,C_V_ele) - fvm::div(dphi_V, C_V_ele)	);
	C_V_ele.relax();
	Current_V_ele = 2*96485*(D_V_ele*fvc::grad(C_H_ele)+grad_phi_ele*C_V_ele);
    eqnResidualV = VEqn.solve().initialResidual(); 
	*/
	//--------------------------update surface concentrations--------------------------
	solve
		(	fvm::ddt(C_V_Asuf)	);
	solve
		(	fvm::ddt(C_H_Asuf)	);
	solve
		(	fvm::ddt(C_V_Csuf)	);
	solve
		(	fvm::ddt(C_H_Csuf)	);
	ratio_surface = Foam::exp(-Segeng_V/0.075307)*C_V_ele.internalField()[0]/(83147-C_V_ele.internalField()[0]);
	C_V_Asuf.internalField()[gridnum_a-1] = ratio_surface*83147/(1+ratio_surface);
	ratio_surface = Foam::exp(-Segeng_H/0.075307)*C_H_ele.internalField()[0]/(83147-C_H_ele.internalField()[0]);
	C_H_Asuf.internalField()[gridnum_a-1] = ratio_surface*83147/(1+ratio_surface);
	ratio_surface = Foam::exp(-Segeng_V/0.075307)*C_V_ele.internalField()[gridnum_ele-1]/(83147-C_V_ele.internalField()[gridnum_ele-1]);
	C_V_Csuf.internalField()[0] = ratio_surface*83147/(1+ratio_surface);
	ratio_surface = Foam::exp(-Segeng_H/0.075307)*C_H_ele.internalField()[gridnum_ele-1]/(83147-C_H_ele.internalField()[gridnum_ele-1]);
	C_H_Csuf.internalField()[0] = ratio_surface*83147/(1+ratio_surface);
//*****************************************************************************************************************//