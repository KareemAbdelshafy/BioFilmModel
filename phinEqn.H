
#	include "phinvariables.H"

	 Info<< "solving phin\n" << endl;

surfaceScalarField phiU = fvc::interpolate(U) & mesh.Sf();

	//Info<< "detlTa \n" << runTime.deltaT().value() << endl;
	//Foam::dimensioned<double> tempMax = Foam::gMax(temp);-Gamma1 * fvc::laplacian(phib)
	Foam::dimensioned<double> phinMax = Foam::gMax(phin);

	//Foam::dimensioned<double> tempMax2 = Foam::gMax(temp.oldTime());
	//Foam::dimensioned<double> phinMax2 = Foam::gMax(phin.oldTime());

	//Info << "tempMin" << tempMax <<endl ;


///---------------- psin Equation --------------- ////

fvScalarMatrix psinEq ( -fvm::laplacian(Lambda*phin*(1-phin), psin) );  

 
fvScalarMatrix psin_phin ( fvm::Sp(scalar(1/dt), phin) + fvm::div(phiU,phin) );  


volScalarField psin_S =  (phin)/(dt) + gn +fvc::laplacian(Lambda*phin*(1-phin), fn1) + fvc::laplacian(Lambda*phin*0.5*(phib-phis), fn2) ;
  
///----------------- phin Equation --------------- ////

fvScalarMatrix phinEq ( - fvm::laplacian(Gamma1  , phin) );   // note I removed the -ve sign

fvScalarMatrix phin_psin (- fvm::Sp(scalar(1) ,psin) );  

/////////////// 

            BlockLduMatrix<vector2> blockM(mesh);
            BlockLduMatrix<vector2> blockM2(mesh);
               
            //- Transfer the coupled interface list for processor/cyclic/etc.
            // boundaries
            blockM.interfaces() = blockTn.boundaryField().blockInterfaces();
            blockM2.interfaces() = blockT2n.boundaryField().blockInterfaces();

            // Grab block diagonal and set it to zero
            Field<tensor2>& d = blockM.diag().asSquare();
            Field<tensor2>& d2 = blockM2.diag().asSquare();
            d = tensor2::zero;
            d2 = tensor2::zero;

            // Grab linear off-diagonal and set it to zero
            Field<tensor2>& l = blockM.lower().asSquare();
            Field<tensor2>& l2 = blockM2.lower().asSquare();
            Field<tensor2>& u = blockM.upper().asSquare();
            Field<tensor2>& u2 = blockM2.upper().asSquare();
            u = tensor2::zero;
            u2 = tensor2::zero;
            l = tensor2::zero;
            l2 = tensor2::zero;

            vector2Field& blockX = blockTn.internalField();
            vector2Field& blockX2 = blockT2n.internalField();
            vector2Field blockB(mesh.nCells(), vector2::zero);
            vector2Field blockB2(mesh.nCells(), vector2::zero);


            //- Inset equations into block Matrix
            blockMatrixTools::insertEquation(0, psinEq, blockM, blockX, blockB);
            blockMatrixTools::insertEquation(1, phinEq, blockM, blockX, blockB);
            blockMatrixTools::insertEquation(0, psin_phin, blockM2, blockX2, blockB2);
            blockMatrixTools::insertEquation(1, phin_psin, blockM2, blockX2, blockB2);

//*
	forAll(u, j )
	{
	u[j](1, 0) = u2[j](1, 1);
	u[j](0, 1) = u2[j](0, 0);
	}

	forAll(l, j )
	{
	l[j](1, 0) = l2[j](1, 1);
	l[j](0, 1) = l2[j](0, 0);
	}
//*/
	forAll(d, i)
	{
	d[i](1,0) = d2[i](1,1);
	d[i](0,1) = d2[i](0,0);
	blockB[i][0] = psin_S[i] *mesh.V()[i];
	blockB[i][1] = 0;
	}


#	include "cyclicBC.H"

            BlockSolverPerformance<vector2> solverPerf =
                BlockLduSolver<vector2>::New
                (
                    blockTn.name(),
                    blockM,
                    mesh.solutionDict().solver(blockTn.name())
                )->solve(blockX, blockB);

	 Foam::VectorN<double,2>& eqnResidualP  = solverPerf.initialResidual();

	 PP = eqnResidualP[0] + eqnResidualP[1];

	solverPerf.print();

	blockMatrixTools::blockRetrieve(0, psin.internalField(), blockX);
	blockMatrixTools::blockRetrieve(1, phin.internalField(), blockX);


	phin.correctBoundaryConditions();
	psin.correctBoundaryConditions();

phinMax = Foam::gMax(phin);
Info << "max phic 1 = "<< phinMax << endl;

//surfaceScalarField phinS = fvc::interpolate(phin);

phis = 1- phin - phib;

#	include "phibEqn.H"

