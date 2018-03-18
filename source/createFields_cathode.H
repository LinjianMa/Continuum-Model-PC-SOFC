Info<< "Reading transportProperties\n" << endl;

IOdictionary cathode_transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        cathodemesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity D\n" << endl;

dimensionedScalar D_H_C
(
    cathode_transportProperties.lookup("D_H_C")
);

dimensionedScalar D_O2
(
    cathode_transportProperties.lookup("D_O2")
);

dimensionedScalar D_O2ads
(
    cathode_transportProperties.lookup("D_O2ads")
);

dimensionedScalar D_O
(
    cathode_transportProperties.lookup("D_O")
);

dimensionedScalar D_OH
(
    cathode_transportProperties.lookup("D_OH")
);

dimensionedScalar D_H2Oad
(
    cathode_transportProperties.lookup("D_H2Oad")
);

dimensionedScalar D_H2O
(
    cathode_transportProperties.lookup("D_H2O")
);

dimensionedScalar D_e_C
(
    cathode_transportProperties.lookup("D_e_C")
);

dimensionedScalar perm_r_C
(
    cathode_transportProperties.lookup("perm_r_C")
);

/***************************************************************/

Info<< "Reading field C\n" << endl;

volScalarField C_H_C
(
    IOobject
    (
        "C_H_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_O2
(
    IOobject
    (
        "C_O2",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_O2ads
(
    IOobject
    (
        "C_O2ads",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_O
(
    IOobject
    (
        "C_O",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_OH
(
    IOobject
    (
        "C_OH",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_H2Oad
(
    IOobject
    (
        "C_H2Oad",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_H2O
(
    IOobject
    (
        "C_H2O",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

Info<< "Reading source term V\n" << endl;

volScalarField V2_C
(
    IOobject
    (
        "V2_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField V3_C
(
    IOobject
    (
        "V3_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField V4_C
(
    IOobject
    (
        "V4_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField V5_C
(
    IOobject
    (
        "V5_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField V6_C
(
    IOobject
    (
        "V6_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

Info<< "Reading/calculating flux \n" << endl;

    volVectorField Flux_H_C
    (
        IOobject
        (
            "Flux_H_C",
            runTime.timeName(),
            cathodemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_H_C*fvc::grad(C_H_C)
    );

    volVectorField Flux_O
    (
        IOobject
        (
            "Flux_O",
            runTime.timeName(),
            cathodemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_O*fvc::grad(C_O)
    );

    volVectorField Flux_OH
    (
        IOobject
        (
            "Flux_OH",
            runTime.timeName(),
            cathodemesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_OH*fvc::grad(C_OH)
    );

/**********************k and overpotential*****************************************/
volScalarField C_H_Csuf
(
    IOobject
    (
        "C_H_Csuf",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField C_V_Csuf
(
    IOobject
    (
        "C_V_Csuf",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField Cmax_C
(
    IOobject
    (
        "Cmax_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField Area_C
(
    IOobject
    (
        "Area_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k2CT_plus
(
    IOobject
    (
        "k2CT_plus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k2CT_minus
(
    IOobject
    (
        "k2CT_minus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k3CT_plus
(
    IOobject
    (
        "k3CT_plus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k3CT_minus
(
    IOobject
    (
        "k3CT_minus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k4CT_plus
(
    IOobject
    (
        "k4CT_plus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k4CT_minus
(
    IOobject
    (
        "k4CT_minus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k5CT_plus
(
    IOobject
    (
        "k5CT_plus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k5CT_minus
(
    IOobject
    (
        "k5CT_minus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k6CT_plus
(
    IOobject
    (
        "k6CT_plus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k6CT_minus
(
    IOobject
    (
        "k6CT_minus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k7CT_plus
(
    IOobject
    (
        "k7CT_plus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField k7CT_minus
(
    IOobject
    (
        "k7CT_minus",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField overpotential_C
(
    IOobject
    (
        "overpotential_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField current_boundary_C
(
    IOobject
    (
        "current_boundary_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField current_boundary_C1
(
    IOobject
    (
        "current_boundary_C1",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField current_boundary_C2
(
    IOobject
    (
        "current_boundary_C2",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField current_boundary_C3
(
    IOobject
    (
        "current_boundary_C3",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField alpha_C
(
    IOobject
    (
        "alpha_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField ratio_C
(
    IOobject
    (
        "ratio_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

/////////////////////////////////////////////////////
volScalarField C_e_C
(
    IOobject
    (
        "C_e_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volScalarField Phi_C
(
    IOobject
    (
        "Phi_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    cathodemesh
);

volVectorField Current_e_C
(
    IOobject
    (
        "Current_e_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
   -D_e_C*fvc::grad(Phi_C)
);

	volScalarField Efield_C
	(
		IOobject
		(
			"Efield_C",
			runTime.timeName(),
			cathodemesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		cathodemesh
	);
	volScalarField Efield_bulk_C
	(
		IOobject
		(
			"Efield_bulk_C",
			runTime.timeName(),
			cathodemesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		cathodemesh
	);

volVectorField grad_phi_C
(
    IOobject
    (
        "grad_phi_C",
        runTime.timeName(),
        cathodemesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
   -D_e_C*fvc::grad(Phi_C)
);

surfaceScalarField dphi_e_C
(
	IOobject
	(
		"dphi_e_C",
		runTime.timeName(),
		cathodemesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
linearInterpolate(grad_phi_C) & cathodemesh.Sf()
);