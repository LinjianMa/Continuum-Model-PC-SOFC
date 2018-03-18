//*****************************************************************************************************************//
	//-----------------------------Charge density calculations---------------------------------
	for (counter2=(gridnum_0);counter2<gridnum_0+gridnum_1;counter2++)
	{
		permittivity = perm_r_A.value() + (-perm_r_A.value()+perm_r.value())*(counter2-gridnum_0)/(gridnum_1);
		charge.internalField()[counter2] = -96485*(C_e_A.internalField()[counter2]-C_eA_equi)/perm_r_A.value();
		//charge.internalField()[counter2] = charge.internalField()[counter2]+1000000*96485/8.31/873*Efield.internalField()[100].y()/permittivity;
	}
	//Surface charge for Pt_BZY
	for (counter2=0;counter2<gridnum_3;counter2++)
	{
		//charge.internalField()[counter2+gridnum_0+gridnum_1] = 96485*(2*(C_V_ele.internalField()[0]-C_V_equi)+(C_H_ele.internalField()[0]-C_H_equi))/perm_r.value();
		charge.internalField()[counter2+gridnum_0+gridnum_1+gridnum_2] = 96485*(2*(C_V_ele.internalField()[counter2]-C_V_equi)+(C_H_ele.internalField()[counter2]-C_H_equi))/perm_r.value();
	}
	//Electrolyte
	for (counter2=0;counter2<gridnum_4;counter2++)
	{
		//charge.internalField()[counter2+grid_ele_L] = 96485*(2*(C_V_ele.internalField()[counter2]-C_V_equi)+(C_H_ele.internalField()[counter2]-C_H_equi))/perm_r.value();
		charge.internalField()[counter2+grid_ele_L+gridnum_3] = 96485*(2*(C_V_ele.internalField()[counter2+gridnum_3]-C_V_equi)+(C_H_ele.internalField()[counter2+gridnum_3]-C_H_equi))/perm_r.value();
	}	
	//Surface charge for Ag_BZY
	for (counter2=0;counter2<gridnum_5;counter2++)
	{
		//charge.internalField()[counter2+grid_ele_L+gridnum_3] = 96485*(2*(C_V_ele.internalField()[gridnum_ele-1]-C_V_equi)+(C_H_ele.internalField()[gridnum_ele-1]-C_H_equi))/perm_r.value();
		charge.internalField()[counter2+grid_ele_L+gridnum_3+gridnum_4] = 96485*(2*(C_V_ele.internalField()[counter2+gridnum_3+gridnum_4]-C_V_equi)+(C_H_ele.internalField()[counter2+gridnum_3+gridnum_4]-C_H_equi))/perm_r.value();
	}
	//Surface charge for Ag
	for (counter2=0;counter2<gridnum_7;counter2++)
	{
		permittivity = perm_r.value() + (-perm_r.value()+perm_r_C.value())*(counter2)/(gridnum_7);
		charge.internalField()[counter2+grid_ele_R+gridnum_6] = -96485*(C_e_C.internalField()[counter2]-C_eC_equi)/perm_r_C.value();
		//charge.internalField()[counter2+grid_ele_R+gridnum_4] = charge.internalField()[counter2+grid_ele_R+gridnum_4]-1000000*96485/8.31/873*Efield.internalField()[grid_ele_R+gridnum_4+gridnum_5+100].y()/permittivity;
	}
//*****************************************************************************************************************//