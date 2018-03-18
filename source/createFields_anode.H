Info<< "Reading transportProperties\n" << endl;

IOdictionary anode_transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        anodemesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity D\n" << endl;

dimensionedScalar D_H2
(
    anode_transportProperties.lookup("D_H2")
);

dimensionedScalar D_H2ads
(
    anode_transportProperties.lookup("D_H2ads")
);

dimensionedScalar D_H_A
(
    anode_transportProperties.lookup("D_H_A")
);

dimensionedScalar D_e_A
(
    anode_transportProperties.lookup("D_e_A")
);

dimensionedScalar perm_r_A
(
    anode_transportProperties.lookup("perm_r_A")
);
/***************************************************************/

Info<< "Reading field C\n" << endl;

volScalarField C_H2
(
    IOobject
    (
        "C_H2",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField C_H2ads
(
    IOobject
    (
        "C_H2ads",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField C_H_A
(
    IOobject
    (
        "C_H_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField C_H_Asuf
(
    IOobject
    (
        "C_H_Asuf",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField C_V_Asuf
(
    IOobject
    (
        "C_V_Asuf",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField C_e_A
(
    IOobject
    (
        "C_e_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField Phi_A
(
    IOobject
    (
        "Phi_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volVectorField Current_e_A
(
    IOobject
    (
        "Current_e_A",
        runTime.timeName(),
        anodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
   -D_e_A*fvc::grad(Phi_A)
);

	volScalarField Efield_A
	(
		IOobject
		(
			"Efield_A",
			runTime.timeName(),
			anodemesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		anodemesh
	);
	volScalarField Efield_bulk_A
	(
		IOobject
		(
			"Efield_bulk_A",
			runTime.timeName(),
			anodemesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		anodemesh
	);

volVectorField grad_phi_A
(
    IOobject
    (
        "grad_phi_A",
        runTime.timeName(),
        anodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
   -D_e_A*fvc::grad(Phi_A)
);

surfaceScalarField dphi_e_A
(
	IOobject
	(
		"dphi_e_A",
		runTime.timeName(),
		anodemesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
linearInterpolate(grad_phi_A) & anodemesh.Sf()
);

Info<< "Reading source term V\n" << endl;
volScalarField V1_A
(
    IOobject
    (
        "V1_A",
        runTime.timeName(),
        anodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField V2_A
(
    IOobject
    (
        "V2_A",
        runTime.timeName(),
        anodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

Info<< "Reading/calculating flux \n" << endl;
    volVectorField Flux_H_A
    (
        IOobject
        (
            "Flux_H_A",
            runTime.timeName(),
            anodemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_H_A*fvc::grad(C_H_A)
    );

/**********************k and overpotential*****************************************/

volScalarField Cmax_A
(
    IOobject
    (
        "Cmax_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField Area_A
(
    IOobject
    (
        "Area_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField k1AT_plus
(
    IOobject
    (
        "k1AT_plus",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField k1AT_minus
(
    IOobject
    (
        "k1AT_minus",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField k2AT_plus
(
    IOobject
    (
        "k2AT_plus",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField k2AT_minus
(
    IOobject
    (
        "k2AT_minus",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField k3AT_plus
(
    IOobject
    (
        "k3AT_plus",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField k3AT_minus
(
    IOobject
    (
        "k3AT_minus",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField overpotential_A
(
    IOobject
    (
        "overpotential_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField current_boundary_A
(
    IOobject
    (
        "current_boundary_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField alpha_A
(
    IOobject
    (
        "alpha_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);

volScalarField ratio_A
(
    IOobject
    (
        "ratio_A",
        runTime.timeName(),
        anodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    anodemesh
);