 Info<< "Reading field phin\n" << endl;

    volScalarField phin
    (
        IOobject
        (
            "phin",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field phib\n" << endl;

    volScalarField phib
    (
        IOobject
        (
            "phib",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


 Info<< "Reading field phis\n" << endl;

    volScalarField phis
    (
        IOobject
        (
            "phis",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1-phin-phib
    );

    Info<< "Reading field c\n" << endl;

    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



/////////////////////////////////////////////////////

    label pdRefCell = 0;
    scalar pdRefValue = 0.0;

    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


   Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );  


    volVectorField Utotal
    (
        IOobject
        (
            "Utotal",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );  



    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd 
    );


 volScalarField pc
    (
        IOobject
        ( 
            "pc",
            phin.time().timeName(),
            phin.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    setRefCell
    (
        pd,
        mesh.solutionDict().subDict("PIMPLE"),
        pdRefCell,
        pdRefValue
    );


    scalar pRefValue = 0.0;

    if (pd.needReference())
    {
        pRefValue = readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("pRefValue")
        );

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pdRefCell)
        );
    }


Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

    Info<< "Reading field divU\n" << endl;

    volScalarField divU
    (
        IOobject
        (
            "divU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(U)
    );

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////



    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Lambda
    (
        transportProperties.lookup("Lambda")
    );

    dimensionedScalar Gamma1
    (
        transportProperties.lookup("Gamma1")
    );

    dimensionedScalar Gamma2
    (
        transportProperties.lookup("Gamma2")
    );

    dimensionedScalar Res
    (
        transportProperties.lookup("Res")
    );

    dimensionedScalar Reb
    (
        transportProperties.lookup("Reb")
    );

    dimensionedScalar Rep
    (
        transportProperties.lookup("Rep")
    );

    dimensionedScalar Ds
    (
        transportProperties.lookup("Ds")
    );

    dimensionedScalar mu0
    (
        transportProperties.lookup("mu0")
    );

    dimensionedScalar Kc
    (
        transportProperties.lookup("Kc")
    );


   Info<< "Reading mu \n" << endl;

    dimensionedScalar Lambda1
    (
        transportProperties.lookup("Lambda1")
    );

   Info<< "Reading C1 \n" << endl;

    dimensionedScalar C1
    (
        transportProperties.lookup("C1")
    );

    dimensionedScalar C2
    (
        transportProperties.lookup("C2")
    );

    dimensionedScalar CB
    (
        transportProperties.lookup("CB")
    );

    dimensionedScalar K1
    (
        transportProperties.lookup("K1")
    );

    dimensionedScalar K2
    (
        transportProperties.lookup("K2")
    );

// Need to be specified 
    dimensionedScalar Np
    (
        transportProperties.lookup("Np")
    );

    dimensionedScalar Nb
    (
        transportProperties.lookup("Nb")
    );

    dimensionedScalar chibn
    (
        transportProperties.lookup("chibn")
    );

    dimensionedScalar chisn
    (
        transportProperties.lookup("chisn")
    );

    dimensionedScalar chisb
    (
        transportProperties.lookup("chisb")
    );

    dimensionedScalar deltaphin
    (
        transportProperties.lookup("deltaphin")
    );

    dimensionedScalar scaling
    (
        transportProperties.lookup("scaling")
    );

////////////////////////////////

    Info<< "Creating field blockTn\n" << endl;
    volVector2Field blockTn
    (
        IOobject
        (
            "blockTn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector2("zero", dimless, vector2::zero)
    );

    Info<< "Creating field blockT2n\n" << endl;
    volVector2Field blockT2n
    (
        IOobject
        (
            "blockT2n",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector2("zero", dimless, vector2::zero)
    );

    Info<< "Creating field blockTb\n" << endl;
    volVector2Field blockTb
    (
        IOobject
        (
            "blockTb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector2("zero", dimless, vector2::zero)
    );

    Info<< "Creating field blockT2b\n" << endl;
    volVector2Field blockT2b
    (
        IOobject
        (
            "blockT2b",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector2("zero", dimless, vector2::zero)
    );

    volScalarField psin
	(
        IOobject
        (
            "psin",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
    );

    volScalarField psib
	(
        IOobject
        (
            "psib",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
    );



/*
    Info<< "Calculating field usource2\n" << endl;
    volVectorField usource2
    (
        IOobject
        (
            "usource2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
 -2 *psi*(fvc::grad(phin))
    );


    Info<< "Calculating field usource2\n" << endl;
    volVectorField usourcenew
    (
        IOobject
        (
            "usourcenew",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
  (fvc::grad(phin))
    );

    Info<< "Calculating field volgn\n" << endl;
    volScalarField volgn
    (
        IOobject
        (
            "volgn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        scaling*mu*c/(K_c+c)
    );


  Info<< "Calculating field temp2\n" << endl;
    volScalarField temp2
    (
        IOobject
        (
            "temp2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
   Lambda*Gamma_2*( (1/N) + (phin.oldTime()/(scalar(1)-phin.oldTime())) - scalar(2)*chi*phin.oldTime())
    );

  Info<< "Calculating field temp\n" << endl;
    volScalarField temp
    (
        IOobject
        (
            "temp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
   Lambda*Gamma_2*( (1/N) + (phin/(scalar(1)-phin)) - scalar(2)*chi*phin)
    );

    Info<< "Calculating field modpsi\n" << endl;
    volScalarField modpsi
    (
        IOobject
        (
            "modpsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	psi
    );



   Info<< "Calculating field fourth\n" << endl;
    volScalarField fourth
    (
        IOobject
        (
            "fourth",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
  fvc::laplacian(Lambda *phin, psi) 
    );

  Info<< "Reading nuA \n" << endl;
    volScalarField nuA
    (
        IOobject
        (
            "nuA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	nuP * max(phin) + nuW * max(scalar(1) - phin ) 
    );

*/
