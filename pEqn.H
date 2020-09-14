{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA,  U, phi)
    );

    adjustPhi(phiU, U, pd);

    phi = phiU +
        (
	 2 * fvc::interpolate(modpsin)*fvc::snGrad(phin) +
	 2 * fvc::interpolate(modpsib)*fvc::snGrad(phib) +
	 fvc::interpolate(modpsib)*fvc::snGrad(phin) + fvc::interpolate(modpsin)*fvc::snGrad(phib)
         // - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();

residulePold = residuleP;

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);
	//pdEqn.relax();

        if ( ~ (flag ==1 && nonOrth == nNonOrthCorr) )
        {
   residuleP =  pdEqn.solve(mesh.solutionDict().solver(pd.name())).initialResidual();
        }
        else
        {
            pdEqn.solve(mesh.solutionDict().solver(pd.name() + "Final"));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pdEqn.flux();
        }
    }

  
	  U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
