Info<< "Reading transportProperties\n" << endl;

IOdictionary electrolyte_transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        electrolytemesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity D\n" << endl;

//electrolyte parameters
dimensionedScalar D_H_ele
(
    electrolyte_transportProperties.lookup("D_H_ele")
);

dimensionedScalar D_V_ele
(
    electrolyte_transportProperties.lookup("D_V_ele")
);

/***************************************************************/

Info<< "Reading field C\n" << endl;

//volScalarField parameters
volScalarField Phi_ele
(
    IOobject
    (
        "Phi_ele",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField C_H_ele 
(
    IOobject
    (
        "C_H_ele",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField C_V_ele 
(
    IOobject
    (
        "C_V_ele",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);


volScalarField ratio_H_anode
(
    IOobject
    (
        "ratio_H_anode",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);
volScalarField ratio_H_cathode
(
    IOobject
    (
        "ratio_H_cathode",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField charge_HA
(
    IOobject
    (
        "charge_HA",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField charge_VA
(
    IOobject
    (
        "charge_VA",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField charge_HC
(
    IOobject
    (
        "charge_HC",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField charge_VC
(
    IOobject
    (
        "charge_VC",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField charge_sumA
(
    IOobject
    (
        "charge_sumA",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField charge_sumC
(
    IOobject
    (
        "charge_sumC",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField I_A
(
    IOobject
    (
        "I_A",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);

volScalarField I_C
(
    IOobject
    (
        "I_C",
        runTime.timeName(),
        electrolytemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    electrolytemesh
);
//end

//volVectorField variables
    volVectorField Current_H_ele
    (
        IOobject
        (
            "Current_H_ele",
            runTime.timeName(),
            electrolytemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_H_ele*fvc::grad(Phi_ele)
    );
    volVectorField Current_V_ele
    (
        IOobject
        (
            "Current_V_ele",
            runTime.timeName(),
            electrolytemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_V_ele*fvc::grad(Phi_ele)
    );

    volVectorField grad_phi_ele
    (
        IOobject
        (
            "grad_phi_ele",
            runTime.timeName(),
            electrolytemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_H_ele*fvc::grad(Phi_ele)
    );

	volScalarField Efield_ele
	(
		IOobject
		(
			"Efield_ele",
			runTime.timeName(),
			electrolytemesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		electrolytemesh
	);
	volScalarField Efield_bulk
	(
		IOobject
		(
			"Efield_bulk",
			runTime.timeName(),
			electrolytemesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		electrolytemesh
	);

	surfaceScalarField dphi_V
	(
		IOobject
		(
			"dphi_V",
			runTime.timeName(),
			electrolytemesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	linearInterpolate(grad_phi_ele) & electrolytemesh.Sf()
	);
	surfaceScalarField dphi_H
	(
		IOobject
		(
			"dphi_H",
			runTime.timeName(),
			electrolytemesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	linearInterpolate(grad_phi_ele) & electrolytemesh.Sf()
	);
