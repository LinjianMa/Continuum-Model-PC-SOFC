	Info<< "Creating meshes\n" << endl;

// **********************************************************************
//  anode mesh
// **********************************************************************
	Info<< "\tCreating anode mesh\n" << endl;
	fvMesh anodemesh
	(
		IOobject
		(
			"anode",
			runTime.timeName(),
			runTime,
			IOobject::MUST_READ
		)
	);

// **********************************************************************
//  cathode mesh
// **********************************************************************
	Info<< "\tCreating cathode mesh\n" << endl;
	fvMesh cathodemesh
	(
		IOobject
		(
			"cathode",
			runTime.timeName(),
			runTime,
			IOobject::MUST_READ
		)
	);

// **********************************************************************
//  electrolyte mesh
// **********************************************************************
	Info<< "\tCreating electrolyte mesh\n" << endl;
	fvMesh electrolytemesh
	(
		IOobject
		(
			"electrolyte",
			runTime.timeName(),
			runTime,
			IOobject::MUST_READ
		)
	);

// **********************************************************************
//  overall mesh
// **********************************************************************
	Info<< "\tCreating overall mesh\n" << endl;
	fvMesh overallmesh
	(
		IOobject
		(
			"overall",
			runTime.timeName(),
			runTime,
			IOobject::MUST_READ
		)
	);