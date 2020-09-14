	scalar dt = runTime.deltaT().value();

	volScalarField temp = (1/Np)*(Foam::log(phin+deltaphin) + 1);  // note the *10
	volScalarField dfdphin =  Gamma2 * (temp + chibn * phib + chisn * phis);
	// where psin = -Gamma1 * fvc::laplacian(phin)

	temp = (1/Nb)*(Foam::log(phib+deltaphin) + 1);
	volScalarField dfdphib =  Gamma2 * (temp + chibn * phin + chisb * phis);
	// where psib = -Gamma1 * fvc::laplacian(phib)

	temp = (Foam::log(phis+deltaphin) + 1);
	volScalarField dfdphis = -Gamma1 * fvc::laplacian(phis) + Gamma2 * (temp + chisn * phin + chisb * phib);

	volScalarField fn1 =  dfdphin - 0.5 *(psib +dfdphib + dfdphis);
	volScalarField fn2 = dfdphis - (psib +dfdphib);

	volScalarField fb1 =  dfdphib - 0.5 *(psin +dfdphin + dfdphis);
	volScalarField fb2 = dfdphis - (psin +dfdphin);

	volScalarField gn = mu0 * phib * c / (Kc + c); 
	volScalarField gb =( (C1*c/(K1 + c)) - CB );
/*

	fvc::ddt(phin) == fvc::laplacian(Lambda*phin*(1-phin), fn1+ psin) + fvc::laplacian(Lambda*phin*0.5*(phib-phis), fn2) + gn;

	fvc::ddt(phib) == fvc::laplacian(Lambda*phib*(1-phib), fb1) + fvc::laplacian(Lambda*phib*0.5*(phin-phis), fb2) + fvc::(gb,phib);

	//////////////////////

	fvm::Sp(scalar(1.5/dt), phin) - (2*phin - 0.5* phin.oldTime())/(dt)  == fvc::laplacian(Lambda*phin*(1-phin), psin) +
	 fvc::laplacian(Lambda*phin*(1-phin), fn1) + fvc::laplacian(Lambda*phin*0.5*(phib-phis), fn2) + gn;

	fvm::Sp(scalar(1.5/dt), phib) - (2*phib - 0.5* phib.oldTime())/(dt)  == fvc::laplacian(Lambda*phib*(1-phib), psib +fb1) + fvc::laplacian(Lambda*phib*0.5*(phin-phis), fb2) + fvc::(gb,phib);

*/



  

