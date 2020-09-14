
//#	include "phinvariables.H"

	 Info<< "solving phib\n" << endl;

	Foam::dimensioned<double> phibMax = Foam::gMax(phib);

///---------------- psib Equation --------------- ////

fvScalarMatrix psibEq ( -fvm::laplacian(Lambda*phib*(1-phib), psib) );  
 
fvScalarMatrix psib_phib ( fvm::Sp(scalar(1/dt), phib) - fvm::Sp(gb, phib) + fvm::div(phiU,phib) );  

volScalarField psib_S =  (phib)/(dt) +fvc::laplacian(Lambda*phib*(1-phib), fb1) + fvc::laplacian(Lambda*phib*0.5*(phin-phis), fb2) ;
  
///----------------- phib Equation --------------- ////

fvScalarMatrix phibEq (- fvm::laplacian(Gamma1  , phib) );   

fvScalarMatrix phib_psib ( - fvm::Sp(scalar(1) ,psib) );  // note I removed the -ve sign

/////////////// 

         //   BlockLduMatrix<vector2> blockM(mesh);
         //   BlockLduMatrix<vector2> blockM2(mesh);
               
            //- Transfer the coupled interface list for processor/cyclic/etc.
            // boundaries
            blockM.interfaces() = blockTb.boundaryField().blockInterfaces();
            blockM2.interfaces() = blockT2b.boundaryField().blockInterfaces();

            // Grab block diagonal and set it to zero
          //  Field<tensor2>& d = blockM.diag().asSquare();
          //  Field<tensor2>& d2 = blockM2.diag().asSquare();
            d = tensor2::zero;
            d2 = tensor2::zero;

            // Grab linear off-diagonal and set it to zero
          //  Field<tensor2>& l = blockM.lower().asSquare();
          //  Field<tensor2>& l2 = blockM2.lower().asSquare();
          //  Field<tensor2>& u = blockM.upper().asSquare();
         //   Field<tensor2>& u2 = blockM2.upper().asSquare();
            u = tensor2::zero;
            u2 = tensor2::zero;
            l = tensor2::zero;
            l2 = tensor2::zero;

         /* vector2Field& */ blockX = blockTb.internalField();
          /* vector2Field& */ blockX2 = blockT2b.internalField();
          // vector2Field&  blockB(mesh.nCells(), vector2::zero);
          // vector2Field& blockB2(mesh.nCells(), vector2::zero);


            //- Inset equations into block Matrix
            blockMatrixTools::insertEquation(0, psibEq, blockM, blockX, blockB);
            blockMatrixTools::insertEquation(1, phibEq, blockM, blockX, blockB);
            blockMatrixTools::insertEquation(0, psib_phib, blockM2, blockX2, blockB2);
            blockMatrixTools::insertEquation(1, phib_psib, blockM2, blockX2, blockB2);

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
	blockB[i][0] = psib_S[i] *mesh.V()[i];
	blockB[i][1] = 0;
	}

	
#	include "cyclicBC2.H"

phibMax = Foam::gMax(phib);
Info << "max phib 2 = "<< phibMax << endl;

            BlockSolverPerformance<vector2> solverPerfb =
                BlockLduSolver<vector2>::New
                (
                    blockTb.name(),
                    blockM,
                    mesh.solutionDict().solver(blockTb.name())
                )->solve(blockX, blockB);

	 Foam::VectorN<double,2>& eqnResidualPb  = solverPerfb.initialResidual();

	 PP = eqnResidualP[0] + eqnResidualP[1];

	solverPerfb.print();

	blockMatrixTools::blockRetrieve(0, psib.internalField(), blockX);
	blockMatrixTools::blockRetrieve(1, phib.internalField(), blockX);


	phib.correctBoundaryConditions();
	psib.correctBoundaryConditions();

phibMax = Foam::gMax(phib);
Info << "max phib 1 ="<< phibMax << endl;

//surfaceScalarField phinS = fvc::interpolate(phin);

phis = 1- phin - phib;

volScalarField volgc = C2 * phib / (K2 + c);

	 Info<< "solving c\n" << endl;

	   solve
            (
                fvm::ddt(phis,c) ==
                fvm::laplacian(Ds*phis,c)
	      - fvm::Sp(volgc,c)
            );

phibMax = Foam::gMax(c);
Info << "max c 1 = "<< phibMax << endl;

