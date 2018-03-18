//*****************************************************************************************************************//
	//Surface charge for Pt
	for (counter2=(grid_ele_L-gridnum_1);counter2<grid_ele_L;counter2++)
	{
		charge.internalField()[counter2] = -aaa;
	}
	//Surface charge for Pt_BZY
	for (counter2=0;counter2<gridnum_2;counter2++)
	{
		charge.internalField()[counter2+grid_ele_L] = aaa;
	}
	//Surface charge for Ag_BZY
	for (counter2=gridnum_2+gridnum_3;counter2<gridnum_ele;counter2++)
	{
		charge.internalField()[counter2+grid_ele_L] = -aaa;
	}
	//Surface charge for Ag
	for (counter2=0;counter2<gridnum_5;counter2++)
	{
		charge.internalField()[counter2+grid_ele_R] = aaa;
	}
//*****************************************************************************************************************//