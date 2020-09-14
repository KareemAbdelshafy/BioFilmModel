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

    dimensionedScalar nuW
    (
        transportProperties.lookup("nuW")
    );

    dimensionedScalar nuP
    (
        transportProperties.lookup("nuP")
    );


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar Lambda
    (
        transportProperties.lookup("Lambda")
    );

    Info<< "Reading Gamma_1 \n" << endl;

    dimensionedScalar Gamma_1
    (
        transportProperties.lookup("Gamma_1")
    );

   Info<< "Reading Gamma_2 \n" << endl;

    dimensionedScalar Gamma_2
    (
        transportProperties.lookup("Gamma_2")
    );

   Info<< "Reading chi \n" << endl;

    dimensionedScalar chi
    (
        transportProperties.lookup("chi")
    );

   Info<< "Reading N \n" << endl;

    dimensionedScalar N
    (
        transportProperties.lookup("N")
    );

    Info<< "Reading D_s \n" << endl;

    dimensionedScalar D_s
    (
        transportProperties.lookup("D_s")
    );

    Info<< "Reading A \n" << endl;

    dimensionedScalar A
    (
        transportProperties.lookup("A")
    );

   Info<< "Reading mu \n" << endl;

    dimensionedScalar mu
    (
        transportProperties.lookup("mu")
    );

   Info<< "Reading K_1 \n" << endl;

    dimensionedScalar K_1
    (
        transportProperties.lookup("K_1")
    );

  Info<< "Reading K_c \n" << endl;

    dimensionedScalar K_c
    (
        transportProperties.lookup("K_c")
    );

  Info<< "Reading deltaphin \n" << endl;

    dimensionedScalar deltaphin
    (
        transportProperties.lookup("deltaphin")
    );

 Info<< "Reading scaling \n" << endl;

    dimensionedScalar scaling
    (
        transportProperties.lookup("scaling")
    );

////////////////////////////////
/*
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
    Info<< "Creating field blockT\n" << endl;
    volVector2Field blockT
    (
        IOobject
        (
            "blockT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector2("zero", dimless, vector2::zero)
    );

    Info<< "Creating field blockT2\n" << endl;
    volVector2Field blockT2
    (
        IOobject
        (
            "blockT2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector2("zero", dimless, vector2::zero)
    );

    volScalarField psi
	(
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
    );

    volVectorField gradpsi
    (
        IOobject
        (
            "gradpsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
  fvc::grad(psi) 
    );

    volVectorField gradphin
    (
        IOobject
        (
            "gradphin",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
  fvc::grad(phin) 
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
