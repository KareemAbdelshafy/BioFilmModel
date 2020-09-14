//fvScalarMatrix phinTot (phinEq + phin_psi);

 
blockM.interfaces() = blockTb.boundaryField().blockInterfaces();
blockM2.interfaces() = blockT2b.boundaryField().blockInterfaces();

	forAll(mesh.boundaryMesh(),patchI)
{
	if (blockM.interfaces().set(patchI))
	{
		Field<tensor2>& coupledLower = blockM.coupleLower()[patchI].asSquare();
		Field<tensor2>& coupledUpper = blockM.coupleUpper()[patchI].asSquare();
           	//coupledLower = tensor2::zero;
           	//coupledUpper = tensor2::zero;
 
		const scalarField& TLower = psibEq.internalCoeffs()[patchI];
		const scalarField& TUpper = psibEq.boundaryCoeffs()[patchI];
		const scalarField& TsLower = phibEq.internalCoeffs()[patchI];
		const scalarField& TsUpper = phibEq.boundaryCoeffs()[patchI];

//		blockMatrixTools::blockInsert(0, TLower, coupledLower);
//		blockMatrixTools::blockInsert(1, TsLower, coupledLower);
//		blockMatrixTools::blockInsert(0, TUpper, coupledUpper);
//		blockMatrixTools::blockInsert(1, TsUpper, coupledUpper);


		Field<tensor2>& coupledLower2 = blockM2.coupleLower()[patchI].asSquare();
		Field<tensor2>& coupledUpper2 = blockM2.coupleUpper()[patchI].asSquare();
//           	coupledLower2 = tensor2::zero;
//           	coupledUpper2 = tensor2::zero;

		const scalarField& TLower2 = psib_phib.internalCoeffs()[patchI];
		const scalarField& TUpper2 = psib_phib.boundaryCoeffs()[patchI];
		const scalarField& TsLower2 = phib_psib.internalCoeffs()[patchI];
		const scalarField& TsUpper2 = phib_psib.boundaryCoeffs()[patchI];

//		blockMatrixTools::blockInsert(0, TLower2, coupledLower2);
//		blockMatrixTools::blockInsert(1, TsLower2, coupledLower2);
//		blockMatrixTools::blockInsert(0, TUpper2, coupledUpper2);
//		blockMatrixTools::blockInsert(1, TsUpper2, coupledUpper2);

	forAll(coupledLower, j )
	{
	coupledLower[j](1, 0) = coupledLower2[j](1, 1);
	coupledLower[j](0, 1) = coupledLower2[j](0, 0);
	}

	forAll(coupledUpper, j )
	{
	coupledUpper[j](1, 0) = coupledUpper2[j](1, 1);
	coupledUpper[j](0, 1) = coupledUpper2[j](0, 0);
	}

//	Info << "scalarField coupledLower " << coupledLower << endl;
//	Info << "scalarField coupledUpper " << coupledUpper << endl;

	}
}


/*

	forAll(mesh.boundaryMesh(),patchII)
{
	if (blockM2.interfaces().set(patchII))
	{
		Field<vector2>& coupledLower2 = blockM2.coupleLower()[patchII].asLinear();
		Field<vector2>& coupledUpper2 = blockM2.coupleUpper()[patchII].asLinear();

		const scalarField& TLower2 = psi_phin.internalCoeffs()[patchII];
		const scalarField& TUpper2 = psi_phin.boundaryCoeffs()[patchII];
		const scalarField& TsLower2 = phin_psi.internalCoeffs()[patchII];
		const scalarField& TsUpper2 = phin_psi.boundaryCoeffs()[patchII];

		blockMatrixTools::blockInsert(0, TLower2, coupledLower2);
		blockMatrixTools::blockInsert(1, TsLower2, coupledLower2);
		blockMatrixTools::blockInsert(0, TUpper2, coupledUpper2);
		blockMatrixTools::blockInsert(1, TsUpper2, coupledUpper2);



	Info << "scalarField coupledLower " << coupledLower << endl;
	Info << "scalarField coupledUpper " << coupledUpper << endl;

	}
}

*/


