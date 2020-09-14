

	volVectorField gradphin = fvc::grad(phin);
	surfaceVectorField gradphinf = 2* Gamma1 * fvc::interpolate(gradphin);
	volScalarField modpsin = fvc::div(gradphinf & mesh.Sf());

	volVectorField gradphib = fvc::grad(phib);
	surfaceVectorField gradphibf = 2* Gamma1 * fvc::interpolate(gradphib);
	volScalarField modpsib = fvc::div(gradphibf & mesh.Sf());	
///////////////////////////////////////////////////////
	volScalarField nu2 = (1/Rep) * phin + (1/Reb) * phib + (1/Res) * (scalar(1) - phin - phib ) ;

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
	gMax(phin)
    );

	nuA =   gMax(nu2);
	surfaceScalarField nu2f = fvc::interpolate(nu2);
	surfaceScalarField nuAf = fvc::interpolate(nuA);

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvc::laplacian(nu2f, U)
      - fvm::laplacian(nuAf, U)
      + fvc::laplacian(nuAf, U)
      - (fvc::grad(U) & fvc::grad(nu2f))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (momentumPredictor)
    {
  residuleU=      solve
        (
            UEqn
         ==
            fvc::reconstruct
            (  
                (
	          2 * fvc::interpolate(modpsin)*fvc::snGrad(phin) +
	          2 * fvc::interpolate(modpsib)*fvc::snGrad(phib) +
	          fvc::interpolate(modpsib)*fvc::snGrad(phin) + fvc::interpolate(modpsin)*fvc::snGrad(phib)
                  //- ghf*fvc::snGrad(rho)
                  - fvc::snGrad(pd)
                )*mesh.magSf()
            )
        ).initialResidual();
    }


