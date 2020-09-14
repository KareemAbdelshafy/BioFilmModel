/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    bioFilm3

Description
    my first try to solve biofilm model without navier-stokes equations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fieldTypes.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "fvBlockMatrix.H"
#include "blockMatrixTools.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readPIMPLEControls.H"
#   include "initContinuityErrs.H"

if (Foam::gMax(phib) < 1.1e-3)
phib = 10 * phib;
// Tolerance variables

	 scalar PP = 0;  scalar minphin = 0; scalar residuleP = 1; scalar residuleU = 1;
scalar gh =0; scalar rho = 1;  scalar residulePold = 1;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    while (runTime.loop())
    {

	Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPIMPLEControls.H"
#       include "readTimeControls.H"
#       include "CourantNo.H"
//#       include "setDeltaT.H"

        int oCorr = 0;  int corr = 0; int flag=0;

#           include "phinEqn.H"

#           include "UEqn.H"


            do {
            
	if (residuleP < 1e-5 ) flag=1;

#               include "pEqn.H"

            }while (flag == 0);

		flag=0;
		residuleP = 1;
		residulePold = 1;

#           include "continuityErrs.H"

            p = pd;

            if (pd.needReference())
            {
                p += dimensionedScalar
                (
                    "p",
                    p.dimensions(),
                    pRefValue - getRefCellValue(p, pdRefCell)
                );
            }
        

divU = fvc::div(U); 

Utotal= -Lambda*(1-phin)*fvc::grad(fn1) - Lambda*0.5*(phib-phis) * fvc::grad(fn2) - Lambda*(1-phin)* fvc::grad(psin) + U; 

//Info << Utotal << endl;
        runTime.write();
//}
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
