Info<< "Reading transportProperties\n" << endl;

IOdictionary overall_transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        overallmesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading permittivity \n" << endl;

dimensionedScalar perm_0
(
    overall_transportProperties.lookup("perm_0")
);

dimensionedScalar perm_r
(
    overall_transportProperties.lookup("perm_r")
);


/***************************************************************/

Info<< "Reading field C\n" << endl;

//electrolyte parameters
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        overallmesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    overallmesh
);

volScalarField charge
(
    IOobject
    (
        "charge",
        runTime.timeName(),
        overallmesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    overallmesh
);

volVectorField Efield
(
    IOobject
    (
        "Efield",
        runTime.timeName(),
        overallmesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
	-fvc::grad(Phi)
);